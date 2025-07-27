"""
Provides functionality for fragment removal mutation using reverse reactions.
"""
import __future__
import copy
import json
import os
import random
from typing import Any, Dict, List, Optional, Tuple, Union

from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.mutation import MutationBase
from autogrow.plugins.registry_base import plugin_managers
from autogrow.types import Compound
from autogrow.utils.logging import log_warning
import autogrow.utils.mol_object_handling as MOH

class FragmentRemoval(MutationBase):
    """
    A mutation plugin that removes fragments from ligands using reverse
    reactions (to maintain synthetic accessibility).
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the plugin category
            and a list of ArgumentVars.
        """
        return (
            "Fragment Removal Mutation",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Run the Fragment Removal Mutation plugin. Creates new molecules by removing fragments from existing molecules using reverse reactions.",
                ),
                # This plugin uses the same rxn_library_path as FragmentAddition
            ],
        )

    def validate(self, params: dict):
        """
        Validate the provided arguments.

        Args:
            params (dict): A dictionary of parameters to validate.

        Raises:
            ValueError: If rxn_library_path is not provided or is invalid.
        """
        if "rxn_library_path" not in params:
            raise ValueError("rxn_library_path must be provided for FragmentRemoval.")
        if not os.path.exists(params["rxn_library_path"]):
            internal_lib = os.path.join(
                os.path.dirname(__file__),
                "reaction_libraries",
                params["rxn_library_path"],
            )
            if os.path.exists(internal_lib):
                params["rxn_library_path"] = internal_lib
            else:
                raise ValueError(
                    "rxn_library_path is not a valid path. "
                    "Please provide a valid path to the reaction library."
                )

    def setup(self, **kwargs):
        """
        Setup the plugin with provided arguments.
        """
        if not hasattr(self, "reverse_reactions"):
            self._load_and_compile_reverse_reactions()

    def _load_and_compile_reverse_reactions(self):
        """
        Load and compile all reverse reaction strings from the reaction library.
        This is done once during setup to avoid costly recompilation.
        """
        rxn_library_path = self.params["rxn_library_path"]
        rxn_library_file = os.path.join(rxn_library_path, "rxn_library.json")
        try:
            with open(rxn_library_file, "r") as f:
                reaction_dict = json.load(f)
        except Exception as e:
            raise Exception("Failed to load rxn_library.json") from e

        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        self.reverse_reactions: List[Tuple[Any, int]] = []
        for rxn_name, rxn_info in reaction_dict.items():
            if "reverse_reaction_strings" in rxn_info and rxn_info["reverse_reaction_strings"]:
                rxn_num = rxn_info["RXN_NUM"]
                for reverse_smarts in rxn_info["reverse_reaction_strings"]:
                    try:
                        rxn = chemtoolkit.reaction_from_smarts(reverse_smarts)
                        rxn.Initialize()
                        self.reverse_reactions.append((rxn, rxn_num))
                    except Exception:
                        log_warning(
                            f"Could not parse reverse reaction SMARTS in '{rxn_name}': "
                            f"{reverse_smarts}. Skipping this reverse reaction."
                        )

    def run_mutation(
        self, cmpd: Compound
    ) -> Optional[List[Tuple[str, int, Union[str, None]]]]:
        """
        Run the fragment removal mutation on a parent molecule.

        Args:
            cmpd (Compound): The compound to mutate.

        Returns:
            Optional[List[Tuple[str, int, Union[str, None]]]]: A list containing a
            single tuple with the new fragment's SMILES, the original reaction ID,
            and None. Returns None if no valid fragment could be generated.
        """
        mol_to_mutate = self._prepare_mol(cmpd.smiles)
        if mol_to_mutate is None:
            return None

        shuffled_reactions = copy.deepcopy(self.reverse_reactions)
        random.shuffle(shuffled_reactions)

        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        for rxn, rxn_num in shuffled_reactions:
            products_tuple = rxn.RunReactants((mol_to_mutate,))

            if not products_tuple:
                continue

            product_sets = list(products_tuple)
            random.shuffle(product_sets)

            for product_set in product_sets:
                fragments = list(product_set)
                sorted_fragments = []

                if len(fragments) > 1:
                    # Calculate MCS for each fragment against the parent
                    fragment_mcs_scores = []
                    for frag in fragments:
                        # Sanitize fragment before MCS calculation for robustness
                        sane_frag = MOH.check_sanitization(copy.deepcopy(frag))
                        if sane_frag is None:
                            continue  # Skip fragments that don't sanitize

                        mcs_result = chemtoolkit.find_mcs([mol_to_mutate, sane_frag])
                        mcs_size = mcs_result.numAtoms if mcs_result is not None else 0
                        fragment_mcs_scores.append((frag, mcs_size))

                    # Sort fragments by MCS size in descending order
                    fragment_mcs_scores.sort(key=lambda x: x[1], reverse=True)

                    # Get the sorted list of fragments
                    sorted_fragments = [item[0] for item in fragment_mcs_scores]
                else:
                    # If there's only one fragment, no need to sort
                    sorted_fragments = fragments

                for fragment in sorted_fragments:
                    validated_smiles = self._validate_product(fragment, cmpd)
                    if validated_smiles is not None:
                        # Found a valid fragment most similar to original in
                        # terms of MCS, return it
                        return [(validated_smiles, rxn_num, None)]

        # No reaction produced a valid fragment
        return None

    def _prepare_mol(self, smiles: str) -> Optional[Any]:
        """
        Prepare a molecule for reaction by sanitizing it.

        Args:
            smiles (str): The SMILES string of the molecule.

        Returns:
            Optional[Any]: An RDKit molecule object if successful, else None.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        try:
            mol = chemtoolkit.mol_from_smiles(smiles, sanitize=False)
        except Exception:
            return None

        mol = MOH.check_sanitization(mol)
        if mol is None:
            return None

        mol = MOH.try_reprotanation(mol) # Reactions work better with explicit hydrogens
        return mol

    def _validate_product(self, product_mol: Any, parent_info: Compound) -> Optional[str]:
        """
        Validate a reaction product.

        Args:
            product_mol (Any): The RDKit molecule object of the product.
            parent_info (Compound): The parent compound, needed for filters.

        Returns:
            Optional[str]: The SMILES string of the validated product, or None.
        """
        product_mol = MOH.check_sanitization(product_mol)
        if product_mol is None:
            return None

        product_mol = MOH.handle_frag_check(product_mol)
        if product_mol is None:
            return None

        product_mol = MOH.check_for_unassigned_atom(product_mol)
        if product_mol is None:
            return None
        
        product_mol = MOH.try_reprotanation(product_mol)
        if product_mol is None:
            return None
            
        product_mol = MOH.try_deprotanation(product_mol)
        if product_mol is None:
            return None

        product_mol = MOH.check_sanitization(product_mol)
        if product_mol is None:
            return None

        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        product_smiles: str = chemtoolkit.mol_to_smiles(
            product_mol, isomeric_smiles=True
        )

        if product_smiles == parent_info.smiles:
            return None

        # Run through filters
        tmp_predock_cmpd = Compound(smiles=product_smiles, id="tmp")
        assert self.plugin_managers is not None, "Plugin managers not set"
        passed_filter = (
            len(self.plugin_managers.SmilesFilter.run(predock_cmpds=[tmp_predock_cmpd])) > 0
        )

        return product_smiles if passed_filter else None