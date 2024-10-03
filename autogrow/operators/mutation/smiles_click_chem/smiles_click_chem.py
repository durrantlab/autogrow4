"""SMILECLICK Class"""
import __future__

import random
import os
import json
import copy
from typing import Any, Dict, List, Optional, Union

from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter
from autogrow.types import PreDockedCompoundInfo
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH
import autogrow.operators.filter.execute_filters as Filter


class SmilesClickChem(object):
    """    This class will take a molecule and Mutate it by reacting it.    """

    def __init__(
        self,
        rxn_library_variables: List[str],
        list_of_already_made_smiles: List[PreDockedCompoundInfo],
        filter_object_dict: Dict[str, ParentFilter],
    ) -> None:
        """
        init for SmilesClickChem. This will set up all the reaction and
        functional dictionaries required to Mutate a molecular

        Inputs:
        :param list rxn_library_variables: a list of user variables which
            define the rxn_library, rxn_library_file,
            complementary_mol_directory, and function_group_library. ie.
            rxn_library_variables = [params['rxn_library'],
            params['rxn_library_file'],
            params['function_group_library'],params['complementary_mol_directory']]
        :param list list_of_already_made_smiles: a list of lists. Each
            sublist contains info about a smiles made in this generation via
            mutation ie.[['O=C([O-])',
            '(Gen_3_Mutant_37_747+ZINC51)Gen_4_Mutant_15_52']]
        :param dict filter_object_dict: a dictionary of all filter objects
            which are to be applied to the newly created ligands.
        """

        # Unpackage the rxn_library_variables

        rxn_library = rxn_library_variables[0]
        rxn_library_file = rxn_library_variables[1]
        function_group_library = rxn_library_variables[2]
        complementary_mol_dir = rxn_library_variables[3]
        self.reaction_dict = self.retrieve_reaction_dict(rxn_library, rxn_library_file)

        # Retrieve the dictionary containing
        # all the possible ClickChem Reactions
        self.list_of_reaction_names = list(self.reaction_dict.keys())

        self.functional_group_dict = self.retrieve_functional_group_dict(
            rxn_library, function_group_library
        )
        self.complementary_mol_dict = self.retrieve_complementary_dictionary(
            rxn_library, complementary_mol_dir
        )

        # List of already predicted smiles
        self.list_of_already_made_smiles = [
            x.smiles for x in list_of_already_made_smiles
        ]
        # Dictionary containing all Filter class
        # objects to be impossed on the ligand
        self.filter_object_dict = filter_object_dict

    def update_list_of_already_made_smiles(
        self, list_of_already_made_smiles_infos: List[PreDockedCompoundInfo]
    ) -> None:
        """
        This updates the list of Smiles which have been made in this
        generation via mutation.

        Inputs:
        :param list list_of_already_made_smiles_infos: a list of lists. Each sublist
            contains info about a smiles made in this generation via mutation.
            ie. [['O=C([O-])',
            '(Gen_3_Mutant_37_747+ZINC51)Gen_4_Mutant_15_52']]
        """
        list_of_already_made_smiles = [
            x.smiles for x in list_of_already_made_smiles_infos
        ]
        self.list_of_already_made_smiles.extend(list_of_already_made_smiles)

    def rxn_lib_format_json_dict_of_dict(
        self, old_dict: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        json dictionaries  import as type unicode. This script converts all
        the keys and items to strings, with a few specific exceptions. It
        takes both the functional group dictionary and the reaction library.

        The reaction library is a dictionary of dictionary and has a few
        exceptions which are not intended to be strings. ie. the num_reactants
        which converts to interger and functional_groups which convert to a
        list of strings.

        The functional_group_dictionary is simply a dictionary with all items
        and keys needing to be strings.

        Inputs:
        :param dic old_dict: a dictionary of the the reaction library or
            functional groups. This is what is importanted from the .json file.

        Returns:
        :returns: dic new_dict: a dictionary of the the reaction library or
            functional groups where the unicode type items have been replaced with
            the proper python data types.
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

    def retrieve_reaction_dict(
        self, rxn_library: str, rxn_library_file: str
    ) -> Dict[str, Dict[str, Any]]:
        """
        This is where all the chemical reactions for SmartClickChem are
        retrieved. If you want to add more just add a Custom set of reactions
        please add a folder to
        PATH/autogrow/operators/mutation/smiles_click_chem/Reaction_libraries/.
        They should be formatted as a dictionary of dictionary using the same
        format as :
        os.path.join(pwd,"reaction_libraries",
                    "click_chem_rxns","ClickChem_rxn_library.json")

        The reactions are written as SMARTS-reaction strings.

        This dictionary uses the reaction name as the key and the Reaction
        Smarts as the value.

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            params['rxn_library_file']
        :param str rxn_library_file: a PATH to a Custom reaction library file
            formatted in a dictionary of dictionaries. in a .json file. This will
            be a blank string if one choses a predefined rxn_library option.

        Returns:
        :returns: dict reaction_dict: A dictionary containing all the
            reactions for ClickChemistry and all the information required to run
            the reaction
        """
        # Get the JSON file to import the proper reaction library
        pwd = os.path.dirname(__file__)
        if not rxn_library_file:
            if rxn_library == "click_chem_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "click_chem_rxns",
                    "ClickChem_rxn_library.json",
                )
            elif rxn_library == "robust_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "robust_rxns",
                    "Robust_Rxns_rxn_library.json",
                )
            elif rxn_library == "all_rxns":
                rxn_library_file = os.path.join(
                    pwd, "reaction_libraries", "all_rxns", "All_Rxns_rxn_library.json"
                )
            elif rxn_library == "Custom":
                if os.path.exists(rxn_library_file) is False:
                    raise Exception(
                        "Custom rxn_library_file cannot be found. "
                        + "Please check the path: ",
                        rxn_library_file,
                    )
            else:
                raise Exception(
                    "rxn_library is not incorporated into smiles_click_chem.py"
                )

            # Import the proper reaction library JSON file
            try:
                with open(rxn_library_file, "r") as rxn_file:
                    reaction_dict_raw = json.load(rxn_file)
            except Exception as e:
                raise Exception(
                    "rxn_library_file json file not able to be imported."
                    + " Check that the rxn_library is formatted correctly"
                ) from e

        elif type(rxn_library_file) == str:
            if os.path.exists(rxn_library_file) is False:
                raise Exception(
                    "Custom specified rxn_library_file directory can not be found"
                )

            if os.path.isfile(rxn_library_file) is False:
                raise Exception("Custom specified rxn_library_file is not a file")

            try:
                extension = os.path.splitext(rxn_library_file)[1]
            except Exception as e:
                raise Exception(
                    "Custom specified rxn_library_file is not .json file."
                    + " It must be a .json dictionary"
                ) from e

            if extension != ".json":
                raise Exception(
                    "Custom specified rxn_library_file is not .json file."
                    + " It must be a .json dictionary"
                )

            # Import the proper reaction library JSON file
            try:
                with open(rxn_library_file, "r") as rxn_file:
                    reaction_dict_raw = json.load(rxn_file)
            except Exception as exc:
                raise Exception(
                    "Custom specified rxn_library_file json file not able to "
                    + "be imported. Check that the rxn_library is "
                    + "formatted correctly"
                ) from exc

        else:
            raise Exception(
                "Custom specified rxn_library_file directory can not be found"
            )

        # Convert the reaction_dict_raw from unicode to the proper
        return self.rxn_lib_format_json_dict_of_dict(reaction_dict_raw)

    def retrieve_functional_group_dict(
        self, rxn_library: str, function_group_library: str
    ) -> Dict[str, str]:
        """
        This retrieves a dictionary of all functional groups required for the
        respective reactions. This dictionary will be used to identify
        possible reactions.

        This is where all the functional groups which will be used in the
        SmartClickChem reactions are retrieved. If you want to add more just
        add a Custom set of reactions please add a folder to
        PATH/autogrow/operators/mutation/smiles_click_chem/Reaction_libraries/.
        They should be formatted as a dictionary of dictionary using the same
        format as :
        os.path.join(pwd,"reaction_libraries","click_chem_rxns",
                     "ClickChem_functional_groups.json")

        IF YOU CHOSE TO DO A Custom REACTION SET YOU MUST PROVIDE A DICTIONARY
        OF ALL FUNCTIONAL GROUPS IT WILL REACT. IF YOU FORGET TO ADD A
        FUNCTIONAL GROUP TO YOUR Custom DICTIONARY, THE REACTION MAY NEVER BE
        UTILIZED.

        Please note if your functional groups involve stereochemistry
            notations such as '\' please replace with '\\' (all functional
            groups should be formatted as SMARTS)

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            params['function_group_library']
        :param str function_group_library: a PATH to a Custom functional group
            dictionary in a .json file. This will be a blank string if one choses
            a predefined functional groups option.

        Returns:
        :returns: dict functional_group_dict: A dictionary containing all
            SMARTS for identifying the functional groups for ClickChemistry
        """

        # Get the JSON file to import the proper reaction library
        pwd = os.path.dirname(__file__)

        if not function_group_library:
            if rxn_library == "click_chem_rxns":
                function_group_library = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "click_chem_rxns",
                    "ClickChem_functional_groups.json",
                )
            elif rxn_library == "robust_rxns":
                function_group_library = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "robust_rxns",
                    "Robust_Rxns_functional_groups.json",
                )
            elif rxn_library == "all_rxns":
                function_group_library = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "all_rxns",
                    "All_Rxns_functional_groups.json",
                )
            elif rxn_library == "Custom":
                if os.path.exists(function_group_library) is False:
                    raise Exception(
                        "Custom function_group_library cannot be found. "
                        + "Please check the path: ",
                        function_group_library,
                    )
            else:
                raise Exception(
                    "rxn_library is not incorporated into smiles_click_chem.py"
                )

            # Import the proper function_group_library JSON file
            try:
                with open(function_group_library, "r") as func_dict_file:
                    functional_group_dict_raw = json.load(func_dict_file)
            except Exception as e:
                raise Exception(
                    "function_group_library json file not able to be imported. "
                    + "Check that the rxn_library is formatted correctly"
                ) from e

        elif type(function_group_library) == str:
            if os.path.exists(function_group_library) is False:
                raise Exception(
                    "Custom specified function_group_library directory can not be found"
                )

            if os.path.isfile(function_group_library) is False:
                raise Exception("Custom specified function_group_library is not a file")

            try:
                extension = os.path.splitext(function_group_library)[1]
            except Exception as e:
                raise Exception(
                    "Custom specified function_group_library is not .json "
                    + "file. It must be a .json dictionary"
                ) from e

            if extension != ".json":
                raise Exception(
                    "Custom specified function_group_library is not .json "
                    + "file. It must be a .json dictionary"
                )

            # Import the proper function_group_library JSON file
            try:
                with open(function_group_library, "r") as func_dict_file:
                    functional_group_dict_raw = json.load(func_dict_file)
            except Exception as exc:
                raise Exception(
                    "function_group_library json file not able to be imported."
                    + " Check that the rxn_library is formatted correctly"
                ) from exc
        else:
            raise Exception(
                "Custom specified function_group_library directory can not be found"
            )

        # Convert the reaction_dict_raw from unicode to the proper
        return self.rxn_lib_format_json_dict_of_dict(functional_group_dict_raw)

    def rand_key_list(self, dictionary: Dict[Any, Any]) -> List[Any]:
        """
        Get a random ordered list of all the keys from  a dictionary.

        Inputs:
        :param dict dictionary: any dictionary

        Returns:
        :returns: list keys: a randomly ordered list containing all the keys
            from the dictionary
        """
        keys = list(dictionary.keys())  # List of keys
        random.shuffle(keys)
        return keys

    def retrieve_complementary_dictionary(
        self, rxn_library: str, complementary_mol_dir: str
    ) -> Dict[str, str]:
        """
        Based on user controlled variables, this definition will retrieve a
        dictionary of molecules separated into classes by their functional
        groups. The sorting of a .smi file into this should be handled in the
        user parameter testing when autogrow is initially started.

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            params['complementary_mol_dir']
        :param dict complementary_mol_dir: the path to the
            complementary_mol_dir directory. It may be an empty string in which
            case the complementary_mol_dir directory will default to those of the
            rxn_library

        Returns:
        :returns: dict complementary_mols_dict: a dictionary of complementary molecules
        """
        script_dir = os.path.dirname(os.path.realpath(__file__))

        if complementary_mol_dir == "":
            if rxn_library == "click_chem_rxns":
                complementary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "click_chem_rxns",
                    "complementary_mol_dir",
                )
            elif rxn_library == "robust_rxns":
                complementary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "robust_rxns",
                    "complementary_mol_dir",
                )
            elif rxn_library == "all_rxns":
                complementary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "all_rxns",
                    "complementary_mol_dir",
                )
            elif rxn_library == "Custom":
                if os.path.isdir(complementary_mol_dir) is False:
                    raise Exception(
                        "Custom complementary_mol_dir cannot be found. "
                        + "Please check the path: ",
                        complementary_mol_dir,
                    )
            else:
                raise Exception(
                    "rxn_library is not incorporated into smiles_click_chem.py"
                )

        elif os.path.isdir(complementary_mol_dir) is False:
            raise Exception(
                "complementary_mol_dir is not a directory. It must be a \
                    directory with .smi files containing SMILES specified by \
                    functional groups.These .smi files must be named the same \
                    as the files in the complementary_mol_dir."
            )

        # Make a list of all the functional groups. These will be the name of
        # the .smi folders already separated by group.
        functional_groups = self.functional_group_dict.keys()

        missing_smi_files = []
        complementary_mols_dict = {}
        for group in functional_groups:
            filepath = f"{complementary_mol_dir}{os.sep}{group}.smi"

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

    def make_reactant_order_list(
        self,
        substructure_search_results: List[int],
        has_substructure_matches_count: int,
    ) -> List[int]:
        """
        make an ordered list of reactants which composed of 0 and 1. This list
        will be used (in later steps) to determine which reactant is the
        ligand and which requires a complementary molecule.

        Inputs:
        :param list substructure_search_results: list composed of 0 and 1. 1
            for if it has the substructure 0 for not
        :param int has_substructure_matches_count: how many substructure
            matches there are
        Returns:
        :returns: list reactant_order_list: an ordered list of reactants which
            composed of 0 and 1.
        """
        # create a list to be used to determine which reactants need
        # complementary mol and which will use the Ligand
        reactant_order_list = []

        # for mols w at least 1 substructure
        if has_substructure_matches_count == 1:
            reactant_order_list = substructure_search_results
        elif has_substructure_matches_count > 1:
            # if more than 1 reactant is found in the ligand than we need to
            # randomly pick 1 to be the molecule in the reaction and the
            # other(s) to be mols chosen from the complementary molecule
            # dictionary

            chosen_as_mol_num = random.randint(0, has_substructure_matches_count - 1)
            counter_of_matches = 0
            for substructure_search_result in substructure_search_results:
                if (
                    substructure_search_result == 1
                    and counter_of_matches == chosen_as_mol_num
                ):
                    reactant_order_list.append(1)
                    counter_of_matches = counter_of_matches + 1
                elif substructure_search_result == 1:
                    reactant_order_list.append(0)
                    counter_of_matches = counter_of_matches + 1
                else:
                    reactant_order_list.append(0)
        return reactant_order_list

    def get_random_complementary_mol(self, functional_group: str) -> List[str]:
        """
        This function will get a dictionary of complementary mols

        Inputs:
        :param str functional_group: the functional group of the needed
            complementary molecule for the reaction

        Returns:
        :returns: list random_comp_mol: list with the SMILES string and name
            of molecule for the randomly chosen comp mol
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

    def determine_functional_groups_in_mol(
        self, mol_deprotanated: Chem.Mol, mol_reprotanated: Chem.Mol
    ) -> List[str]:
        """
        This function will take a molecule and find which functional groups it
        has. This will save time for picking reactions, particularly as
        reaction lists become larger.

        Inputs:
        :param rdkit.Chem.rdchem.Mol mol_deprotanated: an rdkit molecule which
            has been sanitized and deprotanated
        :param rdkit.Chem.rdchem.Mol mol_reprotanated: an rdkit molecule which
            has been sanitized and fully protanated

        Returns:
        :returns: list list_subs_within_mol: a list of the name of every
            functional group found within the molecule. these will be used later
            to filter for reactions.
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

    def run_smiles_click(
        self, ligand_smiles_string: str
    ) -> Optional[List[Union[str, int, None]]]:
        """
        This will take the shuffled list of reaction names
        (self.shuffled_reaction_list) and test the Ligand to see if it is
        capable of being used in the reaction. If the ligand is unable to be
        used in the reaction, then we move on to the next reaction in the
        list. If none work, we return a  None.

        Inputs:
        :param str ligand_smiles_string: SMILES string of a molecule to be
            reacted

        Returns:
        :returns: list product_info: list containing the reaction product, the
            id_number of the reaction as found in the reaction_dict and the id for
            the complementary mol (None if it was a single reactant reaction)
            [reaction_product_smilestring, reaction_id_number,
            zinc_database_comp_mol_name]. returns None if all reactions failed or
            input failed to convert to a sanitizable rdkit mol.
        """
        try:
            mol = Chem.MolFromSmiles(
                ligand_smiles_string, sanitize=False
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
        list_subs_within_mol = self.determine_functional_groups_in_mol(
            mol_deprotanated, mol_reprotanated
        )
        if len(list_subs_within_mol) == 0:
            print(f"{ligand_smiles_string} had no functional groups to react with.")
            return None

        shuffled_reaction_list = self.rand_key_list(
            self.reaction_dict
        )  # Randomize the order of the list of reactions

        tries = 0
        is_rxn_complete = False
        # go through all possible rxns in dictionary of rxns using the random
        # order of rxns loop ends when a rxn is successful or when it runs out
        # of reactions
        reaction_product_smilestring = None
        a_reaction_dict = None
        comp_mol_id = None
        reaction_product = None

        while tries < len(shuffled_reaction_list) and not is_rxn_complete:
            reaction_name = shuffled_reaction_list[tries]
            a_reaction_dict = self.reaction_dict[reaction_name]

            fun_groups_in_rxn = a_reaction_dict["functional_groups"]
            contains_group = None
            for i in range(len(fun_groups_in_rxn)):
                if fun_groups_in_rxn[i] in list_subs_within_mol:
                    contains_group = i
                    # The number i which contains_group is now equal to will
                    # be used to remember the placement of the molecule later
                    # in the reaction.
                    break

                continue

            if contains_group is None:
                # Reaction doesn't contain a functional group found in the
                # reactant molecule. So lets move on to the next molecule
                tries += 1
                continue

            # Determine whether to react using the protanated or
            # deprotanated form of the ligand
            substructure = Chem.MolFromSmarts(
                self.functional_group_dict[fun_groups_in_rxn[i]]
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
                # "Try reaction"
                zinc_database_comp_mol_name = None
                comp_mol_id = None
                try:
                    # if reaction works keep it
                    reaction_products_list = [
                        x[0] for x in rxn.RunReactants((mol_to_use,))
                    ]

                    # randomly shuffle the lists of products so that we don't
                    # bias a single product type. ie ClickChem Reactions
                    # 5_Alkyne_and_Azide produces two products: a 1,5 isomer
                    # and a 1,4 isomer; This will shuffle the list and try
                    # each option
                    random.shuffle(reaction_products_list)

                    if (
                        reaction_products_list not in [(), []]
                        and reaction_products_list
                    ):
                        is_rxn_complete = False
                        for reaction_product in reaction_products_list:
                            # Filter and check the product is valid
                            reaction_product_smilestring = self.check_if_product_is_good(
                                reaction_product
                            )
                            if reaction_product_smilestring is None:
                                is_rxn_complete = False
                            else:
                                # REACTION WORKED!
                                is_rxn_complete = True
                                break
                        if (
                            reaction_product_smilestring is not None
                            and is_rxn_complete is True
                        ):
                            # REACTION WORKED!
                            break
                    # if reaction fails then lets move on to the next
                    # reaction
                    tries += 1
                except Exception:
                    # if reaction fails then lets move on to the next reaction
                    mol_to_use = None
                    tries += 1
                    break
            else:
                # for each functional group in the reaction, test
                # if the ligand has that as a substructure

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
                        functional_group_name = str(
                            a_reaction_dict["functional_groups"][i]
                        )

                        # Determine whether to react using the protanated or
                        # deprotanated form of the ligand
                        substructure = Chem.MolFromSmarts(
                            self.functional_group_dict[fun_groups_in_rxn[i]]
                        )

                        # lets give up to 100 tries to find a comp molecule
                        # which is viable
                        for _ in range(100):
                            # find that group in the complementary dictionary.
                            # comp_molecule = ["cccc", "ZINC123"]
                            comp_molecule = self.get_random_complementary_mol(
                                functional_group_name
                            )

                            # zinc_database name
                            zinc_database_comp_mol_name = comp_molecule[1]

                            # Smiles String of complementary molecule
                            comp_smiles_string = comp_molecule[0]

                            # check this is a santizable molecule
                            comp_mol = Chem.MolFromSmiles(
                                comp_smiles_string, sanitize=False
                            )
                            # try sanitizing, which is necessary later
                            comp_mol = MOH.check_sanitization(comp_mol)

                            # Try with deprotanated molecule rdkit to
                            # recognize for the reaction
                            comp_mol = MOH.try_deprotanation(comp_mol)
                            if comp_mol is None:
                                continue

                            if comp_mol.HasSubstructMatch(substructure) is True:
                                comp_mol = comp_mol
                                # append to ordered list
                                list_reactant_mols.append(comp_mol)
                                comp_mol_id.append(zinc_database_comp_mol_name)
                                break

                            # Try with deprotanated molecule rdkit to
                            # recognize for the reaction
                            comp_mol = MOH.try_deprotanation(comp_mol)
                            if comp_mol is None:
                                continue

                            if comp_mol.HasSubstructMatch(substructure) is True:
                                comp_mol = comp_mol
                                # append to ordered list
                                list_reactant_mols.append(comp_mol)
                                comp_mol_id.append(zinc_database_comp_mol_name)
                                break

                            comp_mol = None
                            continue

                # we will make a tuple of the molecules as rdkit mol objects
                # 1st we generate a list of reactant mol objects then we
                # convert to tuple

                # convert list to tuple
                tuple_reactant_mols = tuple(list_reactant_mols)

                # Run the reaction: We use a try/except statement incase an
                # error occurs and rdkit is unable to complete the reaction.
                # without this a failure to complete the reaction would result
                # in the terminating.

                # Try to run reaction
                try:
                    # if reaction works keep it
                    reaction_products_list = [
                        x[0] for x in rxn.RunReactants(tuple_reactant_mols)
                    ]

                    # randomly shuffle the lists of products so that we don't
                    # bias a single product type. ie ClickChem Reactions
                    # 5_Alkyne_and_Azide produces two products: a 1,5 isomer
                    # and a 1,4 isomer; This will shuffle the list and try
                    # each option
                    random.shuffle(reaction_products_list)

                except Exception:
                    reaction_product = None
                    tries += 1
                    continue

                if reaction_products_list in [(), []] or not reaction_products_list:
                    reaction_id_number = a_reaction_dict["RXN_NUM"]
                    tries += 1
                    continue
                else:
                    is_rxn_complete = False
                    for reaction_product in reaction_products_list:
                        # Filter and check the product is valid
                        reaction_product_smilestring = self.check_if_product_is_good(
                            reaction_product
                        )
                        if reaction_product_smilestring is None:
                            is_rxn_complete = False
                        else:
                            # REACTION WORKED!
                            is_rxn_complete = True
                            break
                    if (
                        reaction_product_smilestring is not None
                        and is_rxn_complete is True
                    ):
                        # REACTION WORKED!
                        break
                    # try again
                    tries += 1

        # end of the big while loop (while tries < len(shuffled_reaction_list)
        # and is_rxn_complete is False)

        # check that a reaction was successful
        if is_rxn_complete is True:
            # NOTE: Don't extract below to new function.
            reaction_product = MOH.check_sanitization(reaction_product)
            if reaction_product is None:
                return None

            reaction_product_smilestring = Chem.MolToSmiles(
                reaction_product, isomericSmiles=True
            )
            assert a_reaction_dict is not None, "a_reaction_dict is None"
            reaction_id_number = a_reaction_dict["RXN_NUM"]

            # RETURNS THE NEW PRODUCTS SMILESTRING, THE REACTION ID NUMBER (SO
            # ONE CAN TRACK THE MOLS LINEAGE). THE COMP_MOL ZINC DATABASE ID
            # NUMBER (IF IT WAS A RXN WITH ONLY 1 REACTANT THIS IS None)
            if comp_mol_id is None:
                zinc_database_comp_mol_names = None
            elif len(comp_mol_id) == 1:
                zinc_database_comp_mol_names = comp_mol_id[0]
            else:
                zinc_database_comp_mol_names = "+".join(comp_mol_id)
            return [
                reaction_product_smilestring,
                reaction_id_number,
                zinc_database_comp_mol_names,
            ]
        # reaction failed
        return None

    def check_if_product_is_good(self, reaction_product):
        """
        This function will test whether the product passes all of the
            requirements:
            1) Mol sanitizes
            2) It isn't in the self.list_of_already_made_smiles
            3) It passes Filters
        Returns the smile if it passes; returns None if it fails.

        Inputs:
        :param rdkit.Chem.rdchem.Mol reaction_product: an rdkit
            molecule to be checked.
        Returns:
        :returns: str reaction_product_smilestring:
            this will return either a SMILES string if it is a good molecule
            or None if it can not sanitize and be cleaned
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
        reaction_product_smilestring = Chem.MolToSmiles(
            reaction_product, isomericSmiles=True
        )
        if reaction_product_smilestring in self.list_of_already_made_smiles:
            return None

        # Run through filters
        passed_filter = Filter.run_filter_on_just_smiles(
            reaction_product_smilestring, self.filter_object_dict
        )
        if passed_filter == False:  # NOTE: Keep as "== False"
            return None
        # passes
        return reaction_product_smilestring
