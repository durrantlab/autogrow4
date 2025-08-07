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
from autogrow.plugins.mutation.utils import (
    load_reaction_library,
    prepare_mol_for_reaction,
    validate_product,
    validate_rxn_library_path,
)
from autogrow.plugins.registry_base import plugin_managers
from autogrow.types import Compound
from autogrow.utils.logging import log_warning
import autogrow.utils.mol_object_handling as MOH

# Set to True to enable image-based debugging for this plugin
DEBUG = False
MIN_FRAGMENT_HEAVY_ATOMS = 5


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
        validate_rxn_library_path(params)

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
        reaction_dict = load_reaction_library(rxn_library_path)

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
        if DEBUG:
            os.makedirs("./debug", exist_ok=True)
            from rdkit.Chem import Draw

        mol_to_mutate = prepare_mol_for_reaction(cmpd.smiles)
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
                if not fragments:
                    continue

                # Calculate MCS score for each fragment against the parent
                fragment_mcs_scores = []
                for frag in fragments:
                    sane_frag = MOH.check_sanitization(copy.deepcopy(frag))
                    if sane_frag is None:
                        continue  # Skip fragments that don't sanitize
                    mcs_result = chemtoolkit.find_mcs([mol_to_mutate, sane_frag])
                    mcs_size = mcs_result.numAtoms if mcs_result is not None else 0
                    fragment_mcs_scores.append((frag, mcs_size))

                if not fragment_mcs_scores:
                    continue  # No sanitizable fragments

                # Sort fragments by MCS size to find the largest
                fragment_mcs_scores.sort(key=lambda x: x[1], reverse=True)

                largest_mcs_fragment = fragment_mcs_scores[0][0]

                kept_fragments = [largest_mcs_fragment]

                # Filter other fragments by heavy atom count
                for frag, mcs_size in fragment_mcs_scores:
                    if frag is largest_mcs_fragment:
                        continue

                    sane_frag = MOH.check_sanitization(copy.deepcopy(frag))
                    if sane_frag is None:
                        continue

                    heavy_atom_count = chemtoolkit.lipinski_heavy_atom_count(sane_frag)
                    if heavy_atom_count >= MIN_FRAGMENT_HEAVY_ATOMS:
                        kept_fragments.append(frag)

                # Remove duplicates. It's possible the largest MCS fragment was also
                # added again due to passing the size filter.
                unique_smiles = set()
                unique_kept_fragments = []
                for frag in kept_fragments:
                    smi = chemtoolkit.mol_to_smiles(frag)
                    if smi not in unique_smiles:
                        unique_smiles.add(smi)
                        unique_kept_fragments.append(frag)

                kept_fragments = unique_kept_fragments

                random.shuffle(kept_fragments)

                for fragment in kept_fragments:
                    validated_smiles = validate_product(
                        fragment, cmpd, self.plugin_managers
                    )
                    if validated_smiles is not None:
                        # Found a valid fragment
                        if DEBUG:
                            self._save_debug_image(
                                parent_mol=mol_to_mutate,
                                fragment_mcs_scores=fragment_mcs_scores,
                                chosen_fragment=fragment,
                                parent_compound=cmpd,
                                rxn_num=rxn_num,
                            )
                        return [(validated_smiles, rxn_num, None)]

        # No reaction produced a valid fragment
        return None

    def _save_debug_image(
        self,
        parent_mol: Any,
        fragment_mcs_scores: List[Tuple[Any, int]],
        chosen_fragment: Any,
        parent_compound: Compound,
        rxn_num: int,
    ):
        """
        Saves a debug image showing the parent, all fragments, and the chosen fragment.

        Args:
            parent_mol (Any): The RDKit molecule object of the parent.
            fragment_mcs_scores (List[Tuple[Any, int]]): A list of tuples, where each
                tuple contains a fragment molecule and its MCS score with the parent.
            chosen_fragment (Any): The RDKit molecule object of the chosen fragment.
            parent_compound (Compound): The parent Compound object.
            rxn_num (int): The reaction number used.
        """
        mols_to_draw = []
        legends = []

        # Add parent molecule
        mols_to_draw.append(parent_mol)
        legends.append(f"Parent: {parent_compound.id}")

        # Add all generated fragments with their MCS scores
        for i, (frag, score) in enumerate(fragment_mcs_scores):
            mols_to_draw.append(frag)
            legends.append(f"Fragment {i+1} (MCS: {score})")

        # Add the chosen fragment separately for emphasis
        mols_to_draw.append(chosen_fragment)
        legends.append("Chosen Fragment")

        # Create the image grid
        img = Draw.MolsToGridImage(
            mols_to_draw,
            legends=legends,
            molsPerRow=4,
            subImgSize=(300, 300),
            useSVG=False,
        )

        # Save the image to a unique file
        filename = f"./debug/{parent_compound.id}_rxn{rxn_num}_{random.randint(1000, 9999)}.png"
        img.save(filename)
        log_warning(f"Saved fragment removal debug image to: {filename}")

