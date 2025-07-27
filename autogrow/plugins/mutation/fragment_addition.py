"""
Provides functionality for fragment addition mutation.

This module implements the FragmentAddition class, a mutation plugin that adds
fragments to existing molecules using specified reaction libraries. It includes
methods for loading and managing reaction libraries, identifying functional
groups, and performing chemical reactions on molecules.

Key components:
- FragmentAddition: Main class for performing fragment addition mutations.
- Reaction library management: Functions for loading and validating reaction
  libraries.
- Molecule preparation: Methods for preparing molecules for reactions.
- Reaction execution: Functions for trying single and multi-reactant reactions.
"""

import contextlib
import json
import os
from typing import Any, Dict, List, Optional, Tuple, Union
from autogrow.config.argument_vars import ArgumentVars

from autogrow.plugins.mutation import MutationBase
from autogrow.plugins.mutation.utils import (
    load_reaction_library,
    validate_product,
    validate_rxn_library_path,
)
from autogrow.types import Compound
from autogrow.utils.logging import log_debug, log_warning
import copy
import random
import glob

import autogrow.utils.mol_object_handling as MOH
from autogrow.plugins.registry_base import plugin_managers
import gzip


class FragmentAddition(MutationBase):
    """Plugin that dds fragments to existing mols using reaction libraries."""

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        Returns:
         Tuple[str, List[ArgumentVars]]: A tuple containing the plugin category
          and a list of ArgumentVars.
        """
        built_in_libs = glob.glob(
            os.path.join(os.path.dirname(__file__), "reaction_libraries") + "/*"
        )

        # Keep only dirs
        built_in_libs = [x for x in built_in_libs if os.path.isdir(x)]

        # Keep only basenames
        built_in_libs = [os.path.basename(x) for x in built_in_libs]

        rxn_library_path_help = (
            "The set of reactions and reactants for use in the mutation operation. Valid values include the path to a reaction library or "
            + (", ".join(built_in_libs))
            + "."
        )
        return (
            "Fragment Addition Mutation",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Run the Fragment Addition Mutation plugin. Creates new molecules by adding fragments to existing molecules, per user-specified reaction libraries.",
                ),
                ArgumentVars(
                    name="rxn_library_path",
                    type=str,
                    default="all_rxns",
                    help=rxn_library_path_help,
                ),
                ArgumentVars(
                    name="min_fragment_mol_weight",
                    type=float,
                    default=None,
                    help=f"The minimum molecular weight of fragments to be used in reactions. If not specified, no minimum MW filter is applied. Used with the {self.name} plugin.",
                ),
                ArgumentVars(
                    name="max_fragment_mol_weight",
                    type=float,
                    default=None,
                    help=f"The maximum molecular weight of fragments to be used in reactions. If not specified, no maximum MW filter is applied. Used with the {self.name} plugin.",
                ),
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

        min_mw = params.get("min_fragment_mol_weight")
        max_mw = params.get("max_fragment_mol_weight")

        if min_mw is not None:
            try:
                min_mw = float(min_mw)
            except (ValueError, TypeError) as e:
                raise ValueError(
                    f"min_fragment_mol_weight must be a number. Got: '{min_mw}'"
                ) from e
            if min_mw < 0:
                raise ValueError("min_fragment_mol_weight must be non-negative.")

        if max_mw is not None:
            try:
                max_mw = float(max_mw)
            except (ValueError, TypeError) as e:
                raise ValueError(
                    f"max_fragment_mol_weight must be a number. Got: '{max_mw}'"
                ) from e
            if max_mw < 0:
                raise ValueError("max_fragment_mol_weight must be non-negative.")

        if min_mw is not None and max_mw is not None and min_mw > max_mw:
            raise ValueError(
                "min_fragment_mol_weight cannot be greater than max_fragment_mol_weight."
            )

    def setup(self, **kwargs):
        """
        Setup the plugin with provided arguments.

        This method is required because one can imagine a scenario where you
        want to setup a plugin only once, then execute the run function multiple
        times.

        Args:
            **kwargs (Any): Arbitrary keyword arguments used for plugin setup.
        """
        self._load_rxn_data()

    def run_mutation(
        self, cmpd: Compound
    ) -> Optional[List[Tuple[str, int, Union[str, None]]]]:
        """
        Run the mutation on the parent molecule.

        This will take the shuffled list of reaction names
        (self.shuffled_reaction_list) and test the molecule to see if it is
        capable of being used in the reaction. If the molecule is unable to be
        used in the reaction, then we move on to the next reaction in the list.
        If none work, we return None.

        Args:
            cmpd (Compound): Compound of a molecule to be reacted.

        Returns:
            Optional[List[Tuple[str, int, Union[str, None]]]]: A list of tuples,
                each containing the reaction product SMILES, the id_number of the
                reaction as found in the reaction_dict, and the id for the
                complementary mol. (None if it was a single reactant reaction).
            Returns None if all reactions failed or input failed to convert to
                a sanitizable rdkit mol.
        """
        # Prepare the molecule
        mol_data = self._prepare_mol(cmpd.smiles)
        if mol_data is None:
            return None
        mol, mol_reprotanated, mol_deprotanated, list_react_grps_in_mol = mol_data

        # Try reactions on the molecule
        reaction_result = self._try_reactions(
            mol_reprotanated, mol_deprotanated, list_react_grps_in_mol, cmpd
        )
        if not reaction_result:
            return None
        return reaction_result

    def _load_rxn_data(self):
        """Load the reaction data, if it hasn't been previously loaded."""
        rxn_library_path = self.params["rxn_library_path"]

        if not hasattr(self, "reaction_dict"):
            # Only load if not already loaded
            self.reaction_dict = load_reaction_library(rxn_library_path)

        if not (hasattr(self, "functional_group_dict")):
            # Only load if not already loaded
            self.functional_group_dict = self._extract_functional_groups(
                self.reaction_dict
            )
        if not hasattr(self, "complementary_mol_dict"):
            # Only load if not already loaded
            self.complementary_mol_dict = self._load_complementary_mols(
                rxn_library_path  # , complementary_mols
            )

        # List of already predicted smiles (to make sure no duplicates)
        # TODO: Think more about this. Needs to be only once, but if here, it will be with every reaction.
        # Now called from setup, so only once, but think more about implementation.
        # existing_smiles: List[Compound],
        # self.existing_smiles = [
        #  x.smiles for x in existing_smiles
        # ]

        self._validate_rxn_lib()

    def _extract_functional_groups(
        self, reaction_dict: Dict[str, Any]
    ) -> Dict[str, str]:
        """
        Extract functional groups and their SMARTS from the reaction dictionary.

        Args:
            reaction_dict (Dict[str, Any]): The dictionary of reactions.

        Returns:
            Dict[str, str]: A dictionary mapping functional group names to SMARTS strings.
        """
        functional_groups = {}
        for rxn_name, rxn_details in reaction_dict.items():
            if "functional_groups" in rxn_details and "group_smarts" in rxn_details:
                for fg_name, fg_smarts in zip(
                    rxn_details["functional_groups"], rxn_details["group_smarts"]
                ):
                    if fg_name not in functional_groups:
                        functional_groups[fg_name] = fg_smarts
        return functional_groups

    def _validate_functional_group_definitions(self):
        """
        Validates that functional groups have consistent SMARTS definitions across all reactions.
        Raises:
            ValueError: If a functional group has multiple, conflicting SMARTS definitions.
        """
        func_grp_defs = {}
        for reaction_name, reaction_info in self.reaction_dict.items():
            for i, func_grp_name in enumerate(reaction_info["functional_groups"]):
                func_grp_smarts = reaction_info["group_smarts"][i]
                if func_grp_name not in func_grp_defs:
                    func_grp_defs[func_grp_name] = func_grp_smarts
                else:
                    # Check if the existing definition matches the new one
                    if func_grp_defs[func_grp_name] != func_grp_smarts:
                        error_msg = (
                            f"Multiple definitions for functional group '{func_grp_name}' "
                            f"found in rxn_library.json. This must be consistent.\n"
                            f"  Reaction '{reaction_name}' defines it as: {func_grp_smarts}\n"
                            f"  Previous definition was: {func_grp_defs[func_grp_name]}"
                        )
                        raise ValueError(error_msg)

    def _validate_rxn_lib(self):
        """
        Validate the reaction library data.

        Raises:
         AssertionError: If any of the validation checks fail.
        """
        for key, val in self.reaction_dict.items():
            # Make sure these keys exist: ['reaction_name',
            # 'example_rxn_product', 'example_rxn_reactants',
            # 'functional_groups', 'group_smarts', 'num_reactants',
            # 'reaction_string', 'RXN_NUM']
            assert all(
                x in val
                for x in [
                    "reaction_name",
                    "example_rxn_product",
                    "example_rxn_reactants",
                    "functional_groups",
                    "group_smarts",
                    "num_reactants",
                    "reaction_string",
                    "RXN_NUM",
                ]
            ), f"Missing key in reaction_dict: {key}"

            # Make sure they are all strings except functional_groups, which should be a list of strings.
            for x in val:
                if x in [
                    "functional_groups",
                    "group_smarts",
                    "example_rxn_reactants",
                    "reverse_reaction_strings",
                ]:
                    continue
                if x == "num_reactants":
                    continue
                assert isinstance(
                    val[x], str
                ), f"Non-string value in reaction_dict: {val[x]}"

            assert isinstance(
                val["functional_groups"], list
            ), f"functional_groups is not a list in reaction_dict: {key}"
            assert all(
                isinstance(x, str) for x in val["functional_groups"]
            ), f"Non-string value in functional_groups in reaction_dict: {key}"
            assert isinstance(
                val["group_smarts"], list
            ), f"group_smarts is not a list in reaction_dict: {key}"
            assert all(
                isinstance(x, str) for x in val["group_smarts"]
            ), f"Non-string value in group_smarts in reaction_dict: {key}"
            assert len(val["functional_groups"]) == len(
                val["group_smarts"]
            ), f"Mismatch between functional_groups and group_smarts for {key}"
            assert isinstance(
                val["num_reactants"], int
            ), "num_reactants is not an integer in reaction_dict"

            # key must match reaction name
            assert (
                key == val["reaction_name"]
            ), f"Key {key} does not match reaction name {val['reaction_name']}"

            for group in val["functional_groups"]:
                # Make sure all functional groups in the reaction are in the functional_group_dict
                assert (
                    group in self.functional_group_dict
                ), f"Functional group {group} not found in functional_group_dict"

                # Make sure self.functional_group_dict[group] is a string
                assert isinstance(
                    self.functional_group_dict[group], str
                ), f"Functional group {group} is not a string"

                # Also make sure the file exists
                assert (
                    group in self.complementary_mol_dict
                ), f"Complementary molecules for functional group {group} not loaded."
                assert isinstance(
                    self.complementary_mol_dict[group], list
                ), f"Complementary molecules for {group} is not a list."

        # After all other checks, validate functional group consistency
        self._validate_functional_group_definitions()

    def _load_complementary_mols(
        self, rxn_library_path: str
    ) -> Dict[str, List[List[str]]]:
        """
        Load and filter the complementary molecules for reactions.
        Based on user-controlled variables, this definition will retrieve a
        dictionary of molecules separated into classes by their functional
        groups. The molecules are filtered by molecular weight if specified.
        Args:
         rxn_library_path (str): A string defining the choice of the reaction
          library.
        Returns:
         Dict[str, List[List[str]]]: A dictionary of complementary molecules where
          keys are functional group names and values are lists of [SMILES, ID]
          pairs.
        Raises:
         Exception: If any required .smi.gz files are missing in the
          complementary_mols directory.
        """
        min_mw = self.params.get("min_fragment_mol_weight")
        max_mw = self.params.get("max_fragment_mol_weight")
        perform_mw_filter = min_mw is not None or max_mw is not None

        complementary_mols_dir = os.path.join(rxn_library_path, "complementary_mols")

        # script_dir = os.path.dirname(os.path.realpath(__file__))

        # # if complementary_mols == "":
        # if rxn_library_path == "click_chem_rxns":
        #     complementary_mols = os.path.join(
        #         script_dir,
        #         "reaction_libraries",
        #         "click_chem_rxns",
        #         "complementary_mols",
        #     )
        # elif rxn_library_path == "robust_rxns":
        #     complementary_mols = os.path.join(
        #         script_dir,
        #         "reaction_libraries",
        #         "robust_rxns",
        #         "complementary_mols",
        #     )
        # elif rxn_library_path == "all_rxns":
        #     complementary_mols = os.path.join(
        #         script_dir,
        #         "reaction_libraries",
        #         "all_rxns",
        #         "complementary_mols",
        #     )
        # else:
        #     raise Exception(
        #         "rxn_library_path is not incorporated into smiles_click_chem.py"
        #     )

        # elif os.path.isdir(complementary_mols) is False:
        #     raise Exception(
        #         "complementary_mols is not a directory. It must be a \
        #             directory with .smi files containing SMILES specified by \
        #             functional groups.These .smi files must be named the same \
        #             as the files in the complementary_mols."
        #     )

        # Make a list of all the functional groups. These will be the name of
        # the .smi folders already separated by group.
        functional_groups = self.functional_group_dict.keys()

        missing_smi_files = []
        complementary_mols_dict = {}
        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        for group in functional_groups:
            filepath = f"{complementary_mols_dir}{os.sep}{group}.smi.gz"
            if not os.path.isfile(filepath):
                missing_smi_files.append(filepath)
                log_warning(
                    f"Could not find the following .smi.gz file for complementary molecules for Mutation: {filepath}"
                )
                continue

            filtered_mols = []
            with gzip.open(filepath, "rt") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < 3:
                        continue
                    smiles, mol_id, mw = parts[0], f"CID-{parts[1]}", float(parts[2])

                    if perform_mw_filter:
                        passes_min = (min_mw is None) or (mw >= min_mw)
                        passes_max = (max_mw is None) or (mw <= max_mw)
                        if passes_min and passes_max:
                            filtered_mols.append([smiles, mol_id])
                        # If mol is None, it's skipped, which is fine.
                    else:
                        # No filtering, just add the molecule
                        filtered_mols.append([smiles, mol_id])

            if perform_mw_filter and not filtered_mols:
                log_warning(
                    f"No fragments for functional group '{group}' passed the MW filter (min: {min_mw}, max: {max_mw})."
                )

            complementary_mols_dict[group] = filtered_mols

        if missing_smi_files:
            raise Exception(
                "The following .smi.gz file for complementary molecules "
                + "for Mutation is missing: ",
                missing_smi_files,
            )

        return complementary_mols_dict

    def _prepare_mol(
        self, ligand_smiles: str
    ) -> Optional[Tuple[Any, Any, Any, List[str]]]:
        """
        Prepare the molecule for reaction.

        Args:
            ligand_smiles (str): SMILES string of a molecule to be reacted.

        Returns:
            Optional[Tuple[Chem.Mol, Chem.Mol, Chem.Mol, List[str]]]: A tuple
                containing the original molecule, reprotanated molecule,
                deprotanated molecule, and a list of functional groups within
                the molecule. Returns None if preparation fails.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        try:
            mol = chemtoolkit.mol_from_smiles(
                ligand_smiles, sanitize=False
            )  # This is the input molecule which serves as the parent molecule
        except Exception:
            # mol object failed to initialize
            return None

        # try sanitizing, which is necessary later
        mol = MOH.check_sanitization(mol)
        if mol is None:
            return None

        # Is important for some functional groups while being deprotanated are
        # useful for other reaction
        mol_reprotanated = copy.deepcopy(mol)
        mol_reprotanated = MOH.try_reprotanation(mol_reprotanated)
        if mol_reprotanated is None:
            return None

        mol_deprotanated = copy.deepcopy(mol)
        mol_deprotanated = MOH.try_deprotanation(mol_deprotanated)
        if mol_deprotanated is None:
            return None

        # Determine which functional groups are within a ligand
        list_subs_within_mol = self._identify_functional_grps(
            mol_deprotanated, mol_reprotanated
        )
        if len(list_subs_within_mol) == 0:
            log_warning(f"Detected no reactive functional groups: {ligand_smiles}")
            return None

        return mol, mol_reprotanated, mol_deprotanated, list_subs_within_mol

    def _try_reactions(
        self,
        mol_reprotanated: Any,
        mol_deprotanated: Any,
        list_subs_within_mol: List[str],
        parent_info: Compound,
    ) -> Optional[List[Tuple[str, int, Optional[str]]]]:
        """
        Try reactions on the molecule.

        Args:
            mol_reprotanated (Chem.Mol): Reprotanated molecule.
            mol_deprotanated (Chem.Mol): Deprotanated molecule.
            list_subs_within_mol (List[str]): List of functional groups within
                the molecule.

        Returns:
            Optional[List[Tuple[str, int, Optional[str]]]]: A list of tuples, each
                containing the reaction product SMILES, reaction ID number, and
                complementary molecule name (if applicable). Returns None if all
                reactions fail.
        """
        # Randomize the order of the list of reactions
        shuffled_reaction_list = self._shuffle_dict_keys(self.reaction_dict)

        for reaction_name in shuffled_reaction_list:
            a_reaction_dict = self.reaction_dict[reaction_name]

            result = self._try_reaction(
                a_reaction_dict,
                mol_deprotanated,
                mol_reprotanated,
                list_subs_within_mol,
                parent_info,
            )
            if result is not None:
                return result

        return None

    def _identify_functional_grps(
        self, mol_deprotanated: Any, mol_reprotanated: Any
    ) -> List[str]:
        """
        Identify functional groups present in a molecule.

        This function will take a molecule and find which functional groups it
        has. This will save time for picking reactions, particularly as reaction
        lists become larger.

        Args:
            mol_deprotanated (Chem.Mol): An rdkit molecule which has been
                sanitized and deprotanated.
            mol_reprotanated (Chem.Mol): An rdkit molecule which has been
                sanitized and fully protanated.

        Returns:
            List[str]: A list of the names of every functional group found
                within the molecule. These will be used later to filter for
                reactions.
        """
        list_subs_within_mol = []
        functional_group_dict = self.functional_group_dict

        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        for key in list(functional_group_dict.keys()):
            substructure = chemtoolkit.mol_from_smarts(functional_group_dict[key])
            if substructure is None:
                raise ValueError(
                    f"Invalid SMARTS string for functional group '{key}': "
                    f"'{functional_group_dict[key]}'. Please check your "
                    f"rxn_library.json file."
                )
            if mol_reprotanated.HasSubstructMatch(substructure):
                list_subs_within_mol.append(key)
            elif mol_deprotanated.HasSubstructMatch(substructure):
                list_subs_within_mol.append(key)
        return list_subs_within_mol

    def _shuffle_dict_keys(self, dictionary: Dict[Any, Any]) -> List[Any]:
        """
        Get a random ordered list of all the keys from a dictionary.

        Args:
            dictionary (Dict[Any, Any]): Any dictionary.

        Returns:
            List[Any]: A randomly ordered list containing all the keys from the
                dictionary.
        """
        keys = list(dictionary.keys())  # List of keys
        random.shuffle(keys)
        return keys

    def _try_reaction(
        self,
        a_reaction_dict: Dict[str, Any],
        mol_deprotanated: Any,
        mol_reprotanated: Any,
        list_subs_within_mol: List[str],
        parent_info: Compound,
    ) -> Optional[List[Tuple[str, int, Optional[str]]]]:
        """
        Try to perform the reaction specified in a_reaction_dict on the mol.

        Args:
            a_reaction_dict (Dict[str, Any]): The reaction to try.
            mol_deprotanated (Chem.Mol): Deprotanated molecule.
            mol_reprotanated (Chem.Mol): Reprotanated molecule.
            list_subs_within_mol (List[str]): List of functional groups within
                the molecule.

        Returns:
            Optional[List[Tuple[str, int, Optional[str]]]]: A list of tuples,
                each containing the reaction product SMILES, reaction ID number, and
                complementary molecule name (if applicable). Returns None if
                reaction fails.
        """
        fun_groups_in_rxn = a_reaction_dict["functional_groups"]
        contains_group = None
        for i in range(len(fun_groups_in_rxn)):
            if fun_groups_in_rxn[i] in list_subs_within_mol:
                contains_group = i
                # The number i which contains_group is now equal to will
                # be used to remember the placement of the molecule later
                # in the reaction.
                break

        if contains_group is None:
            # Reaction doesn't contain a functional group found in the
            # reactant molecule. So lets move on to the next reaction
            return None

        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        # Determine whether to react using the protanated or
        # deprotanated form of the ligand
        substructure = chemtoolkit.mol_from_smarts(
            self.functional_group_dict[fun_groups_in_rxn[contains_group]]
        )

        if mol_deprotanated.HasSubstructMatch(substructure) is True:
            mol_to_use = copy.deepcopy(mol_deprotanated)
        else:
            mol_to_use = copy.deepcopy(mol_reprotanated)
        substructure = None

        rxn = chemtoolkit.reaction_from_smarts(str(a_reaction_dict["reaction_string"]))
        rxn.Initialize()

        # if the reaction requires only a single reactant we will attempt
        # to run the reaction
        if a_reaction_dict["num_reactants"] == 1:
            return self._try_single_reactant_reaction(rxn, mol_to_use, a_reaction_dict, parent_info)
        else:
            return self._try_multi_reactant_reaction(
                rxn, mol_to_use, a_reaction_dict, contains_group, parent_info
            )

    def _try_single_reactant_reaction(
        self, rxn: Any, mol_to_use: Any, a_reaction_dict: Dict[str, Any], parent_info: Compound
    ) -> Optional[List[Tuple[str, int, Optional[str]]]]:
        """
        Try a single reactant reaction.

        Args:
            rxn (AllChem.ChemicalReaction): The reaction object.
            mol_to_use (Chem.Mol): The molecule to react.
            a_reaction_dict (Dict[str, Any]): The reaction dictionary.

        Returns:
            Optional[List[Tuple[str, int, Optional[str]]]]: A list with a single
                tuple containing the reaction product SMILES, reaction ID number,
                and None (as there's no complementary molecule). Returns None if
                reaction fails.
        """
        with contextlib.suppress(Exception):
            # if reaction works keep it
            reaction_products_list = [x[0] for x in rxn.RunReactants((mol_to_use,))]

            # randomly shuffle the lists of products so that we don't
            # bias a single product type. ie ClickChem Reactions
            # 5_Alkyne_and_Azide produces two products: a 1,5 isomer
            # and a 1,4 isomer; This will shuffle the list and try
            # each option
            random.shuffle(reaction_products_list)

            if reaction_products_list:
                for reaction_product in reaction_products_list:
                    # Filter and check the product is valid
                    reaction_product_smiles = validate_product(
                        reaction_product, parent_info, self.plugin_managers
                    )
                    if reaction_product_smiles is not None:
                        reaction_id_number = a_reaction_dict["RXN_NUM"]
                        return [(reaction_product_smiles, reaction_id_number, None)]
        return None

    def _try_multi_reactant_reaction(
        self,
        rxn: Any,
        mol_to_use: Any,
        a_reaction_dict: Dict[str, Any],
        contains_group: int,
        parent_info: Compound,
    ) -> Optional[List[Tuple[str, int, Optional[str]]]]:
        """
        Try a multi-reactant reaction.

        Args:
            rxn (AllChem.ChemicalReaction): The reaction object.
            mol_to_use (Chem.Mol): The molecule to react.
            a_reaction_dict (Dict[str, Any]): The reaction dictionary.
            contains_group (int): Index of the functional group in the reaction
                that is in the molecule.

        Returns:
            Optional[List[Tuple[str, int, Optional[str]]]]: A list of tuples,
                each containing the reaction product SMILES, reaction ID number,
                and complementary molecule name(s). Returns None if reaction
                fails.
        """
        mutants_per_batch = self.params.get("mutants_per_batch", 1)
        products = []
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        fun_groups_in_rxn = a_reaction_dict["functional_groups"]

        for _ in range(mutants_per_batch):
            list_reactant_mols = []
            comp_mol_id = []
            reactants_found_for_this_mutant = True
            for i in range(len(fun_groups_in_rxn)):
                if i == contains_group:
                    # This is where the molecule goes
                    list_reactant_mols.append(mol_to_use)
                else:
                    # for reactants which need to be taken from the
                    # complementary dictionary. Find the reactants functional
                    # group
                    functional_group_name = str(a_reaction_dict["functional_groups"][i])

                    # Determine the substructure
                    substructure = chemtoolkit.mol_from_smarts(
                        self.functional_group_dict[fun_groups_in_rxn[i]]
                    )
                    found_comp_mol = False

                    # Let's give up to 100 tries to find a complementary
                    # molecule
                    for _ in range(100):
                        # Find that group in the complementary dictionary.
                        # comp_molecule = ["cccc", "ZINC123"]
                        comp_molecule = self._get_random_complementary_mol(
                            functional_group_name
                        )
                        if comp_molecule is None:
                            # No fragments available for this group, so we can't form this mutant
                            break

                        # zinc_database name
                        zinc_database_comp_mol_name = comp_molecule[1]
                        # SMILES string of complementary molecule
                        comp_smiles = comp_molecule[0]
                        # Check if this is a sanitizable molecule
                        comp_mol = chemtoolkit.mol_from_smiles(
                            comp_smiles, sanitize=False
                        )
                        # Try sanitizing, which is necessary later
                        comp_mol = MOH.check_sanitization(comp_mol)
                        if comp_mol is None:
                            continue
                        # Try with deprotanated molecule to recognize for the
                        # reaction
                        comp_mol_deprotanated = MOH.try_deprotanation(
                            copy.deepcopy(comp_mol)
                        )
                        if (
                            comp_mol_deprotanated is not None
                            and comp_mol_deprotanated.HasSubstructMatch(substructure)
                        ):
                            list_reactant_mols.append(comp_mol_deprotanated)
                            comp_mol_id.append(zinc_database_comp_mol_name)
                            found_comp_mol = True
                            break
                        # Try with reprotanated molecule
                        comp_mol_reprotanated = MOH.try_reprotanation(
                            copy.deepcopy(comp_mol)
                        )
                        if (
                            comp_mol_reprotanated is not None
                            and comp_mol_reprotanated.HasSubstructMatch(substructure)
                        ):
                            list_reactant_mols.append(comp_mol_reprotanated)
                            comp_mol_id.append(zinc_database_comp_mol_name)
                            found_comp_mol = True
                            break
                    if not found_comp_mol:
                        # Could not find a complementary molecule
                        reactants_found_for_this_mutant = False
                        break
            if not reactants_found_for_this_mutant:
                continue  # to next mutant in batch
            # Convert list to tuple
            tuple_reactant_mols = tuple(list_reactant_mols)
            # Try to run reaction
            try:
                # If reaction works, keep it
                reaction_products_list = [
                    x[0] for x in rxn.RunReactants(tuple_reactant_mols)
                ]
                # Randomly shuffle the list of products to avoid bias
                random.shuffle(reaction_products_list)
            except Exception:
                continue  # to next mutant in batch
            if reaction_products_list:
                for reaction_product in reaction_products_list:
                    # Filter and check if the product is valid
                    reaction_product_smiles = validate_product(
                        reaction_product, parent_info, self.plugin_managers
                    )
                    if reaction_product_smiles is not None:
                        reaction_id_number = a_reaction_dict["RXN_NUM"]
                        if len(comp_mol_id) == 1:
                            zinc_database_comp_mol_names = comp_mol_id[0]
                        else:
                            zinc_database_comp_mol_names = "+".join(comp_mol_id)
                        products.append(
                            (
                                reaction_product_smiles,
                                reaction_id_number,
                                zinc_database_comp_mol_names,
                            )
                        )
                        break  # Move to the next mutant in the batch
        return products if products else None

    def _get_random_complementary_mol(
        self, functional_group: str
    ) -> Optional[List[str]]:
        """
        Get a random complementary molecule for a given functional group.
        Args:
         functional_group (str): The functional group of the needed
          complementary molecule for the reaction.
        Returns:
         Optional[List[str]]: A list containing the SMILES string and name of the
          randomly chosen complementary molecule, or None if no fragments
          are available for the given functional group.
        """
        fragment_list = self.complementary_mol_dict.get(functional_group)
        if not fragment_list:
            log_warning(f"No available fragments for functional group: {functional_group}")
            return None
        return random.choice(fragment_list)