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
from autogrow.config.argparser import ArgumentVars

# from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.mutation import MutationBase
from autogrow.utils.logging import log_debug, log_warning
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
import copy
import random
import glob

import autogrow.utils.mol_object_handling as MOH

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class FragmentAddition(MutationBase):
    """
    A mutation plugin that adds fragments to existing molecules using specified
    reaction libraries.
    """

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
            raise ValueError("rxn_library_path must be provided.")

        # Make sure the rxn_library_path is a valid path
        if os.path.exists(params["rxn_library_path"]) is False:
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

        This method is required because one can imagine a scenario where you
        want to setup a plugin only once, then execute the run function multiple
        times.

        Args:
            **kwargs (Any): Arbitrary keyword arguments used for plugin setup.
        """
        self._load_rxn_data()

    def run_mutation(
        self, parent_smiles: str
    ) -> Optional[Tuple[str, int, Union[str, None]]]:
        """
        Run the mutation on the parent molecule.

        This will take the shuffled list of reaction names
        (self.shuffled_reaction_list) and test the molecule
        to see if it is capable of being used in the reaction. If the molecule
        is unable to be used in the reaction, then we move on to the next
        reaction in the list. If none work, we return None.

        Args:
            parent_smiles (str): SMILES string of a molecule to be reacted.

        Returns:
            Optional[Tuple[str, int, Union[str, None]]]: A tuple containing the
                reaction product SMILES, the id_number of the reaction as found
                in the reaction_dict, and the id for the complementary mol
                (None if it was a single reactant reaction).
            Returns None if all reactions failed or input failed to convert to
                a sanitizable rdkit mol.
        """
        # Prepare the molecule
        mol_data = self._prepare_mol(parent_smiles)
        if mol_data is None:
            return None
        mol, mol_reprotanated, mol_deprotanated, list_subs_within_mol = mol_data

        # Try reactions on the molecule
        reaction_result = self._try_reactions(
            mol_reprotanated, mol_deprotanated, list_subs_within_mol
        )

        if reaction_result is None:
            return None

        (
            reaction_product_smiles,
            reaction_id_number,
            zinc_database_comp_mol_names,
        ) = reaction_result

        return (
            reaction_product_smiles,
            reaction_id_number,
            zinc_database_comp_mol_names,
        )

    def _load_rxn_data(self):
        """
        Load the reaction data, if it hasn't been previously loaded.
        """
        rxn_library_path = self.params["rxn_library_path"]

        if not hasattr(self, "reaction_dict"):
            # Only load if not already loaded
            self.reaction_dict = self._load_rxn_lib(rxn_library_path)

        if not (hasattr(self, "functional_group_dict")):
            # Only load if not already loaded
            self.functional_group_dict = self._load_functional_grps(rxn_library_path)

        if not hasattr(self, "complementary_mol_dict"):
            # Only load if not already loaded
            self.complementary_mol_dict = self._load_complementary_mols(
                rxn_library_path  # , complementary_mols
            )

        # List of already predicted smiles (to make sure no duplicates)
        # TODO: Think more about this. Needs to be only once, but if here, it will be with every reaction.
        # Now called from setup, so only once, but think more about implementation.
        # existing_smiles: List[PreDockedCompound],
        # self.existing_smiles = [
        #     x.smiles for x in existing_smiles
        # ]

        self._validate_rxn_lib()

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
                if x == "functional_groups":
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
                assert os.path.exists(
                    self.complementary_mol_dict[group]
                ), f"Complementary mol file {self.complementary_mol_dict[group]} not found"

        # TODO: Could be more validation (for example of smi files in complementary_mol_dict)

    def _load_rxn_lib(self, rxn_library_path: str) -> Dict[str, Dict[str, Any]]:
        """
        Load the chemical reactions for SmartClickChem.

        This is where all the chemical reactions for SmartClickChem are
        retrieved.

        The reactions are written as SMARTS-reaction strings. This dictionary
        uses the reaction name as the key and the Reaction Smarts as the value.

        Args:
            rxn_library_path (str): A string defining the choice of the reaction
                library.

        Returns:
            Dict[str, Dict[str, Any]]: A dictionary containing all the reactions
                and all the information required to run the reaction.

        Raises:
            Exception: If the rxn_library file cannot be imported or if the
                rxn_library_path is incorrectly formatted.
        """
        return self._load_reformatted_rxn_dict(
            rxn_library_path,
            "rxn_library.json",
            "rxn_library_file json file not able to be imported.",
            " Check that the rxn_library_path is formatted correctly",
        )

    def _load_functional_grps(self, rxn_library_path: str) -> Dict[str, str]:
        """
        Load the functional groups required for the respective reactions.

        This retrieves a dictionary of all functional groups required for the
        respective reactions. This dictionary will be used to identify possible
        reactions.

        Note: If your functional groups involve stereochemistry notations such
        as '\', please replace with '\\' (all functional groups should be
        formatted as SMARTS)

        Args:
            rxn_library_path (str): A string defining the choice of the reaction library.

        Returns:
            Dict[str, str]: A dictionary containing all SMARTS for identifying
                the functional groups.

        Raises:
            Exception: If the function_group_library json file cannot be
                imported or if the rxn_library_path is incorrectly formatted.
        """
        return self._load_reformatted_rxn_dict(
            rxn_library_path,
            "functional_groups.json",
            "function_group_library json file not able to be imported. ",
            "Check that the rxn_library_path is formatted correctly",
        )

    def _load_reformatted_rxn_dict(self, rxn_library_path, arg1, arg2, arg3):
        rxn_library_file = os.path.join(rxn_library_path, arg1)
        try:
            with open(rxn_library_file, "r") as rxn_file:
                reaction_dict_raw = json.load(rxn_file)
        except Exception as e:
            raise Exception((arg2 + arg3)) from e
        return self._reformat_rxn_dict(reaction_dict_raw)

    def _load_complementary_mols(
        self, rxn_library_path: str  # , complementary_mols: str
    ) -> Dict[str, str]:
        """
        Load the complementary molecules for reactions.

        Based on user-controlled variables, this definition will retrieve a
        dictionary of molecules separated into classes by their functional
        groups. The sorting of a .smi file into this should be handled in the
        user parameter testing when autogrow is initially started.

        Args:
            rxn_library_path (str): A string defining the choice of the reaction
                library.

        Returns:
            Dict[str, str]: A dictionary of complementary molecules where keys
                are functional group names and values are paths to corresponding
                .smi files.

        Raises:
            Exception: If any required .smi files are missing in the
                complementary_mols directory.
        """
        complementary_mols = os.path.join(rxn_library_path, "complementary_mols")

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
        for group in functional_groups:
            filepath = f"{complementary_mols}{os.sep}{group}.smi"

            if os.path.isfile(filepath) is True:
                complementary_mols_dict[group] = filepath

            else:
                missing_smi_files.append(filepath)
                print(
                    f"Could not find the following .smi file for complementary  molecules for Mutation: {filepath}"
                )

        if missing_smi_files:
            raise Exception(
                "The following .smi file for complementary molecules "
                + "for Mutation is missing: ",
                missing_smi_files,
            )

        return complementary_mols_dict

    def _reformat_rxn_dict(self, old_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        Convert json dictionary items to appropriate Python data types.

        Json dictionaries import as type unicode. This script converts all the
        keys and items to strings, with a few specific exceptions. It takes both
        the functional group dictionary and the reaction library.

        The reaction library is a dictionary of dictionaries and has a few
        exceptions which are not intended to be strings, i.e., the num_reactants
        which converts to integer and functional_groups which convert to a list
        of strings.

        The functional_group_dictionary is simply a dictionary with all items
        and keys needing to be strings.

        Args:
            old_dict (Dict[str, Any]): A dictionary of the reaction library or
                functional groups. This is what is imported from the .json file.

        Returns:
            Dict[str, Any]: A dictionary of the reaction library or functional
                groups where the unicode type items have been replaced with the
                proper Python data types.
        """
        new_dict = {}
        for rxn_key, rxn_dic_old in old_dict.items():
            key_str = str(rxn_key)

            # For reaction libraries
            if type(rxn_dic_old) == dict:
                new_sub_dict = {}
                for key in rxn_dic_old.keys():
                    sub_key_str = str(key)
                    item = rxn_dic_old[key]

                    if sub_key_str == "functional_groups":
                        new_list = []
                        for i in item:
                            i_str = str(i)
                            new_list.append(i_str)

                        item = new_list
                    elif sub_key_str == "num_reactants":
                        item = int(item)
                    else:
                        item = str(item)

                    new_sub_dict[sub_key_str] = item
                new_dict[key_str] = new_sub_dict

            else:
                item = old_dict[rxn_key]
                new_dict[key_str] = str(item)

        return new_dict

    def _prepare_mol(
        self, ligand_smiles: str
    ) -> Optional[Tuple[Chem.Mol, Chem.Mol, Chem.Mol, List[str]]]:
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
        try:
            mol = Chem.MolFromSmiles(
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
        mol_reprotanated: Chem.Mol,
        mol_deprotanated: Chem.Mol,
        list_subs_within_mol: List[str],
    ) -> Optional[Tuple[str, int, Optional[str]]]:
        """
        Try reactions on the molecule.

        Args:
            mol_reprotanated (Chem.Mol): Reprotanated molecule.
            mol_deprotanated (Chem.Mol): Deprotanated molecule.
            list_subs_within_mol (List[str]): List of functional groups within
                the molecule.

        Returns:
            Optional[Tuple[str, int, Optional[str]]]: A tuple containing the
                reaction product SMILES, reaction ID number, and complementary
                molecule name (if applicable). Returns None if all reactions
                fail.
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
            )
            if result is not None:
                return result

        return None

    def _identify_functional_grps(
        self, mol_deprotanated: Chem.Mol, mol_reprotanated: Chem.Mol
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

        for key in list(functional_group_dict.keys()):
            substructure = Chem.MolFromSmarts(functional_group_dict[key])
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
        mol_deprotanated: Chem.Mol,
        mol_reprotanated: Chem.Mol,
        list_subs_within_mol: List[str],
    ) -> Optional[Tuple[str, int, Optional[str]]]:
        """
        Try to perform the reaction specified in a_reaction_dict on the
        molecule.

        Args:
            a_reaction_dict (Dict[str, Any]): The reaction to try.
            mol_deprotanated (Chem.Mol): Deprotanated molecule.
            mol_reprotanated (Chem.Mol): Reprotanated molecule.
            list_subs_within_mol (List[str]): List of functional groups within
                the molecule.

        Returns:
            Optional[Tuple[str, int, Optional[str]]]: A tuple containing the
                reaction product SMILES, reaction ID number, and complementary
                molecule name (if applicable). Returns None if reaction fails.
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

        # Determine whether to react using the protanated or
        # deprotanated form of the ligand
        substructure = Chem.MolFromSmarts(
            self.functional_group_dict[fun_groups_in_rxn[contains_group]]
        )

        if mol_deprotanated.HasSubstructMatch(substructure) is True:
            mol_to_use = copy.deepcopy(mol_deprotanated)
        else:
            mol_to_use = copy.deepcopy(mol_reprotanated)
        substructure = None

        rxn = AllChem.ReactionFromSmarts(str(a_reaction_dict["reaction_string"]))
        rxn.Initialize()

        # if the reaction requires only a single reactant we will attempt
        # to run the reaction
        if a_reaction_dict["num_reactants"] == 1:
            return self._try_single_reactant_reaction(rxn, mol_to_use, a_reaction_dict)
        else:
            return self._try_multi_reactant_reaction(
                rxn, mol_to_use, a_reaction_dict, contains_group,
            )

    def _try_single_reactant_reaction(
        self,
        rxn: AllChem.ChemicalReaction,
        mol_to_use: Chem.Mol,
        a_reaction_dict: Dict[str, Any],
    ) -> Optional[Tuple[str, int, Optional[str]]]:
        """
        Try a single reactant reaction.

        Args:
            rxn (AllChem.ChemicalReaction): The reaction object.
            mol_to_use (Chem.Mol): The molecule to react.
            a_reaction_dict (Dict[str, Any]): The reaction dictionary.

        Returns:
            Optional[Tuple[str, int, Optional[str]]]: A tuple containing the
                reaction product SMILES, reaction ID number, and None (as
                there's no complementary molecule). Returns None if reaction
                fails.
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
                    reaction_product_smiles = self._validate_product(reaction_product)
                    if reaction_product_smiles is not None:
                        reaction_id_number = a_reaction_dict["RXN_NUM"]
                        return reaction_product_smiles, reaction_id_number, None
        return None

    def _try_multi_reactant_reaction(
        self,
        rxn: AllChem.ChemicalReaction,
        mol_to_use: Chem.Mol,
        a_reaction_dict: Dict[str, Any],
        contains_group: int,
    ) -> Optional[Tuple[str, int, Optional[str]]]:
        """
        Try a multi-reactant reaction.

        Args:
            rxn (AllChem.ChemicalReaction): The reaction object.
            mol_to_use (Chem.Mol): The molecule to react.
            a_reaction_dict (Dict[str, Any]): The reaction dictionary.
            contains_group (int): Index of the functional group in the reaction
                that is in the molecule.

        Returns:
            Optional[Tuple[str, int, Optional[str]]]: A tuple containing the
                reaction product SMILES, reaction ID number, and complementary
                molecule name(s). Returns None if reaction fails.
        """
        fun_groups_in_rxn = a_reaction_dict["functional_groups"]
        list_reactant_mols = []
        comp_mol_id = []
        for i in range(len(fun_groups_in_rxn)):
            if i == contains_group:
                # This is where the molecule goes
                list_reactant_mols.append(mol_to_use)
            else:
                # for reactants which need to be taken from the
                # complementary dictionary. Find the reactants
                # functional group
                functional_group_name = str(a_reaction_dict["functional_groups"][i])

                # Determine the substructure
                substructure = Chem.MolFromSmarts(
                    self.functional_group_dict[fun_groups_in_rxn[i]]
                )

                # Let's give up to 100 tries to find a complementary molecule
                for _ in range(100):
                    # Find that group in the complementary dictionary.
                    # comp_molecule = ["cccc", "ZINC123"]
                    comp_molecule = self._get_random_complementary_mol(
                        functional_group_name
                    )

                    # zinc_database name
                    zinc_database_comp_mol_name = comp_molecule[1]

                    # SMILES string of complementary molecule
                    comp_smiles = comp_molecule[0]

                    # Check if this is a sanitizable molecule
                    comp_mol = Chem.MolFromSmiles(comp_smiles, sanitize=False)
                    # Try sanitizing, which is necessary later
                    comp_mol = MOH.check_sanitization(comp_mol)

                    # Try with deprotanated molecule to recognize for the reaction
                    comp_mol = MOH.try_deprotanation(comp_mol)
                    if comp_mol is None:
                        continue

                    if comp_mol.HasSubstructMatch(substructure) is True:
                        list_reactant_mols.append(comp_mol)
                        comp_mol_id.append(zinc_database_comp_mol_name)
                        break

                    # Try with deprotanated molecule
                    comp_mol = MOH.try_deprotanation(comp_mol)
                    if comp_mol is None:
                        continue

                    if comp_mol.HasSubstructMatch(substructure) is True:
                        list_reactant_mols.append(comp_mol)
                        comp_mol_id.append(zinc_database_comp_mol_name)
                        break
                else:
                    # Could not find a complementary molecule
                    return None

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
            return None

        if reaction_products_list:
            for reaction_product in reaction_products_list:
                # Filter and check if the product is valid
                reaction_product_smiles = self._validate_product(reaction_product)
                if reaction_product_smiles is not None:
                    reaction_id_number = a_reaction_dict["RXN_NUM"]
                    if len(comp_mol_id) == 1:
                        zinc_database_comp_mol_names = comp_mol_id[0]
                    else:
                        zinc_database_comp_mol_names = "+".join(comp_mol_id)
                    return (
                        reaction_product_smiles,
                        reaction_id_number,
                        zinc_database_comp_mol_names,
                    )
        return None

    def _validate_product(self, reaction_product):
        """
        Validate the reaction product.

        This function will test whether the product passes all of the
        requirements:

        1) Mol sanitizes
        2) It passes Filters

        Args:
            reaction_product (Chem.Mol): An rdkit molecule to be checked.

        Returns:
            Optional[str]: The SMILES string of the validated product if it
                passes all checks, or None if it fails any check.
        """
        reaction_product = MOH.check_sanitization(reaction_product)
        if reaction_product is None:
            return None

        # Remove any fragments incase 1 made it through
        reaction_product = MOH.handle_frag_check(reaction_product)
        if reaction_product is None:
            return None

        # Make sure there are no unassigned atoms which made it through. These
        # are very unlikely but possible
        reaction_product = MOH.check_for_unassigned_atom(reaction_product)
        if reaction_product is None:
            return None

        reaction_product = MOH.try_reprotanation(reaction_product)
        if reaction_product is None:
            return None

        # Remove H's
        reaction_product = MOH.try_deprotanation(reaction_product)
        if reaction_product is None:
            return None

        reaction_product = MOH.check_sanitization(reaction_product)
        if reaction_product is None:
            return None

        # Check if product SMILE has been made before
        reaction_product_smiles: str = Chem.MolToSmiles(
            reaction_product, isomericSmiles=True
        )

        # NOTE: I think duplicate smiles are eliminated in the function that
        # calls this plugin. No need to do that here.

        # if reaction_product_smiles in self.existing_smiles:
        #     return None

        # Run through filters
        # passed_filter = Filter.run_filter_on_just_smiles(
        #     reaction_product_smiles, self.filter_object_dict
        # )
        assert self.plugin_managers is not None, "Plugin managers not set"
        passed_filter = (
            len(self.plugin_managers.SmilesFilter.run(smiles=reaction_product_smiles))
            > 0
        )

        if not passed_filter:
            return None

        # passes
        return reaction_product_smiles

    def _get_random_complementary_mol(self, functional_group: str) -> List[str]:
        """
        Get a random complementary molecule for a given functional group.

        Args:
            functional_group (str): The functional group of the needed
                complementary molecule for the reaction.

        Returns:
            List[str]: A list containing the SMILES string and name of the
                randomly chosen complementary molecule.
        """
        infile = self.complementary_mol_dict[functional_group]

        with open(infile, "r") as f:
            random_comp_mol_line = random.choice(f.readlines())
            random_comp_mol_line = (
                random_comp_mol_line.replace("\n", "")
                .replace("\t", " ")
                .replace("    ", " ")
            )
            for _ in range(10):
                random_comp_mol_line.replace("  ", " ")
            parts = random_comp_mol_line.split(
                " "
            )  # split line into parts separated by 4-spaces
            # parts = [x for x in random_comp_mol_line.split(" ") if x!= ""]
            # # split line into parts separated by 4-spaces

            smile_list = parts[0]
            zinc_name_list = parts[1]
            random_comp_mol = [smile_list, zinc_name_list]

        return random_comp_mol
