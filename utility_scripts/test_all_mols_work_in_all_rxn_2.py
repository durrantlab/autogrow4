"""SMILECLICK Class"""
import __future__

import random
import os
import json
import copy
import argparse
import sys

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


import support_scripts.Multiprocess as mp
import support_scripts.mol_object_handling as MOH

class SmilesClickChem(object):
    """
    This class will take a molecule and Mutate it by reacting it.
    
    This is modified from the AutoGrow source code file: 
        /autogrow4/autogrow/operators/mutation/smiles_click_chem/smiles_click_chem.py
    Filters were removed for simplicity
    """

    def __init__(self, rxn_library_variables, list_of_already_made_smiles):
        """
        init for SmilesClickChem. This will set up all the reaction and
        functional dictionaries required to Mutate a molecular

        Inputs:
        :param list rxn_library_variables: a list of user variables which
            define the rxn_library, rxn_library_file,
            complimentary_mol_directory, and function_group_library. ie.
            rxn_library_variables = [vars['rxn_library'],
            vars['rxn_library_file'],
            vars['function_group_library'],vars['complimentary_mol_directory']]
        :param list list_of_already_made_smiles: a list of lists. Each
            sublist contains info about a smiles made in this generation via
            mutation ie.[['O=C([O-])',
            '(Gen_3_Mutant_37_747+ZINC51)Gen_4_Mutant_15_52']]
        """

        # Unpackage the rxn_library_variables

        rxn_library = rxn_library_variables[0]
        rxn_library_file = rxn_library_variables[1]
        function_group_library = rxn_library_variables[2]
        complimentary_mol_dir = rxn_library_variables[3]
        self.reaction_dict = self.retrieve_reaction_dict(
            rxn_library, rxn_library_file
        )
        # Retrieve the dictionary containing
        # all the possible ClickChem Reactions
        self.list_of_reaction_names = list(self.reaction_dict.keys())

        self.functional_group_dict = self.retrieve_functional_group_dict(
            rxn_library, function_group_library
        )
        self.complimentary_mol_dict = self.retrieve_complimentary_dictionary(
            rxn_library, complimentary_mol_dir
        )

        # List of already predicted smiles
        self.list_of_already_made_smiles = [x[0] for x in list_of_already_made_smiles]

    def update_list_of_already_made_smiles(self, list_of_already_made_smiles):
        """
        This updates the list of Smiles which have been made in this
        generation via mutation.

        Inputs:
        :param list list_of_already_made_smiles: a list of lists. Each sublist
            contains info about a smiles made in this generation via mutation.
            ie. [['O=C([O-])',
            '(Gen_3_Mutant_37_747+ZINC51)Gen_4_Mutant_15_52']]
        """
        list_of_already_made_smiles = [x[0] for x in list_of_already_made_smiles]
        self.list_of_already_made_smiles.extend(list_of_already_made_smiles)

    def rxn_lib_format_json_dict_of_dict(self, old_dict):
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
            functional groups. This is what is importated from the .json file.

        Returns:
        :returns: dic new_dict: a dictionary of the the reaction library or
            functional groups where the unicode type items have been replaced with
            the proper python data types.
        """
        new_dict = {}
        for rxn_key in old_dict.keys():
            rxn_dic_old = old_dict[rxn_key]
            key_str = str(rxn_key)

            # For reaction libraries
            if type(rxn_dic_old) == dict:
                new_sub_dict = {}
                for key in rxn_dic_old.keys():
                    sub_key_str = str(key)
                    item = rxn_dic_old[key]

                    if sub_key_str == "num_reactants":
                        item = int(item)
                    elif sub_key_str == "functional_groups":
                        new_list = []
                        for i in item:
                            i_str = str(i)
                            new_list.append(i_str)

                        item = new_list
                    else:
                        item = str(item)

                    new_sub_dict[sub_key_str] = item
                new_dict[key_str] = new_sub_dict

            # For functional groups
            else:
                item = old_dict[rxn_key]
                new_dict[key_str] = str(item)

        return new_dict

    def retrieve_reaction_dict(self, rxn_library, rxn_library_file):
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
            vars['rxn_library_file']
        :param str rxn_library_file: a PATH to a Custom reaction library file
            formated in a dictionary of dictionaries. in a .json file. This will
            be a blank string if one choses a predefined rxn_library option.

        Returns:
        :returns: dict reaction_dict: A dictionary containing all the
            reactions for ClickChemistry and all the information required to run
            the reaction
        """
        # Get the JSON file to import the proper reaction library
        pwd = os.path.dirname(__file__)
        if rxn_library_file == "":

            if rxn_library == "click_chem_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "click_chem_rxns",
                    "ClickChem_rxn_library.json"
                )
            elif rxn_library == "robust_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "robust_rxns",
                    "Robust_Rxns_rxn_library.json"
                )
            elif rxn_library == "all_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "all_rxns",
                    "All_Rxns_rxn_library.json"
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
            except:
                raise Exception(
                    "rxn_library_file json file not able to be imported."
                    + " Check that the rxn_library is formated correctly"
                )

        elif type(rxn_library_file) == str:
            if os.path.exists(rxn_library_file) is False:
                raise Exception(
                    "Custom specified rxn_library_file directory can not be found"
                )

            if os.path.isfile(rxn_library_file) is False:
                raise Exception(
                    "Custom specified rxn_library_file is not a file"
                )

            try:
                extension = os.path.splitext(rxn_library_file)[1]
            except:
                raise Exception(
                    "Custom specified rxn_library_file is not .json file."
                    + " It must be a .json dictionary"
                )

            if extension != ".json":
                raise Exception(
                    "Custom specified rxn_library_file is not .json file."
                    + " It must be a .json dictionary"
                )

            # Import the proper reaction library JSON file
            try:
                with open(rxn_library_file, "r") as rxn_file:
                    reaction_dict_raw = json.load(rxn_file)
            except:
                raise Exception(
                    "Custom specified rxn_library_file json file not able to "
                    + "be imported. Check that the rxn_library is "
                    + "formated correctly"
                )

        else:
            raise Exception(
                "Custom specified rxn_library_file directory can not be found"
            )

        # Convert the reaction_dict_raw from unicode to the proper
        reaction_dict = self.rxn_lib_format_json_dict_of_dict(reaction_dict_raw)

        return reaction_dict

    def retrieve_functional_group_dict(self, rxn_library, function_group_library):
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
            groups should be formated as SMARTS)

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            vars['function_group_library']
        :param str function_group_library: a PATH to a Custom functional group
            dictionary in a .json file. This will be a blank string if one choses
            a predefined functional groups option.

        Returns:
        :returns: dict functional_group_dict: A dictionary containing all
            SMARTS for identifying the functional groups for ClickChemistry
        """

        # Get the JSON file to import the proper reaction library
        pwd = os.path.dirname(__file__)

        if function_group_library == "":

            if rxn_library == "click_chem_rxns":
                function_group_library = os.path.join(
                    pwd, "reaction_libraries",
                    "click_chem_rxns",
                    "ClickChem_functional_groups.json",
                )
            elif rxn_library == "robust_rxns":
                function_group_library = os.path.join(
                    pwd, "reaction_libraries",
                    "robust_rxns",
                    "Robust_Rxns_functional_groups.json",
                )
            elif rxn_library == "all_rxns":
                function_group_library = os.path.join(
                    pwd, "reaction_libraries",
                    "all_rxns", "All_Rxns_functional_groups.json",
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
            except:
                raise Exception(
                    "function_group_library json file not able to be imported. "
                    + "Check that the rxn_library is formated correctly"
                )

        elif type(function_group_library) == str:
            if os.path.exists(function_group_library) is False:
                raise Exception(
                    "Custom specified function_group_library directory can not be found"
                )

            if os.path.isfile(function_group_library) is False:
                raise Exception("Custom specified function_group_library is not a file")

            try:
                extension = os.path.splitext(function_group_library)[1]
            except:
                raise Exception(
                    "Custom specified function_group_library is not .json "
                    + "file. It must be a .json dictionary"
                )

            if extension != ".json":
                raise Exception(
                    "Custom specified function_group_library is not .json "
                    + "file. It must be a .json dictionary"
                )

            # Import the proper function_group_library JSON file
            try:
                with open(function_group_library, "r") as func_dict_file:
                    functional_group_dict_raw = json.load(func_dict_file)
            except:
                raise Exception(
                    "function_group_library json file not able to be imported."
                    + " Check that the rxn_library is formated correctly"
                )
        else:
            raise Exception(
                "Custom specified function_group_library directory can not be found"
            )

        # Convert the reaction_dict_raw from unicode to the proper
        functional_group_dict = self.rxn_lib_format_json_dict_of_dict(
            functional_group_dict_raw
        )

        return functional_group_dict

    def rand_key_list(self, dictionary):
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

    def retrieve_complimentary_dictionary(self, rxn_library, complimentary_mol_dir):
        """
        Based on user controled variables, this definition will retrieve a
        dictionary of molecules seperated into classes by their functional
        groups. The sorting of a .smi file into this should be handled in the
        user parameter testing when autogrow is initailly started.

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            vars['complimentary_mol_dir']
        :param dict complimentary_mol_dir: the path to the
            complimentary_mol_dir directory. It may be an empty string in which
            case the complimentary_mol_dir directory will default to those of the
            rxn_library

        Returns:
        :returns: dict complimentary_mols_dict: a dictionary of complimentary molecules
        """
        script_dir = os.path.dirname(os.path.realpath(__file__))

        if complimentary_mol_dir == "":
            if rxn_library == "click_chem_rxns":
                complimentary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "click_chem_rxns",
                    "complimentary_mol_dir",
                )
            elif rxn_library == "robust_rxns":
                complimentary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "robust_rxns",
                    "complimentary_mol_dir",
                )
            elif rxn_library == "all_rxns":
                complimentary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "all_rxns",
                    "complimentary_mol_dir",
                )
            elif rxn_library == "Custom":
                if os.path.isdir(complimentary_mol_dir) is False:
                    raise Exception(
                        "Custom complimentary_mol_dir cannot be found. "
                        + "Please check the path: ",
                        complimentary_mol_dir,
                    )
            else:
                raise Exception(
                    "rxn_library is not incorporated into smiles_click_chem.py"
                )

        else:
            if os.path.isdir(complimentary_mol_dir) is False:
                raise Exception(
                    "complimentary_mol_dir is not a directory. It must be a \
                    directory with .smi files containing SMILES specified by \
                    functional groups.These .smi files must be named the same \
                    as the files in the complimentary_mol_dir."
                )

        # Make a list of all the functional groups. These will be the name of
        # the .smi folders already seperated by group.
        functional_groups = self.functional_group_dict.keys()

        missing_smi_files = []
        complimentary_mols_dict = {}
        for group in functional_groups:
            filepath = "{}{}{}.smi".format(complimentary_mol_dir, os.sep, group)

            if os.path.isfile(filepath) is True:
                complimentary_mols_dict[group] = filepath

            else:
                missing_smi_files.append(filepath)
                print(
                    "Could not find the following .smi file for complimentary "
                    + " molecules for Mutation: {}".format(filepath)
                )

        if len(missing_smi_files) != 0:
            raise Exception(
                "The following .smi file for complimentary molecules "
                + "for Mutation is missing: ",
                missing_smi_files,
            )

        return complimentary_mols_dict

    def make_reactant_order_list(self, substructure_search_result,
                                 has_substructure_matches_count):
        """
        make an ordered list of reactants which composed of 0 and 1. This list
        will be used (in later steps) to determine which reactant is the
        ligand and which requires a complimentary molecule.

        Inputs:
        :param list substructure_search_result: list composed of 0 and 1. 1
            for if it has the substructure 0 for not
        :param int has_substructure_matches_count: how many substructure
            matches there are
        Returns:
        :returns: list reactant_order_list: an ordered list of reactants which
            composed of 0 and 1.
        """
        # for mols w atleast 1 substructure
        if has_substructure_matches_count == 1:
            reactant_order_list = substructure_search_result
        elif has_substructure_matches_count > 1:
            # if more than 1 reactant is found in the ligand than we need to
            # randomly pick 1 to be the molecule in the reaction and the
            # other(s) to be mols chosen from the complimentary molecule
            # dictionary

            # create a list to be used to determine which reactants need
            # complimentary mol and which will use the Ligand
            reactant_order_list = []

            chosen_as_mol_num = random.randint(0, has_substructure_matches_count - 1)
            counter_of_matches = 0
            for i in range(0, len(substructure_search_result)):
                if (
                        substructure_search_result[i] == 1
                        and counter_of_matches == chosen_as_mol_num
                ):
                    reactant_order_list.append(1)
                    counter_of_matches = counter_of_matches + 1
                elif (
                        substructure_search_result[i] == 1
                        and counter_of_matches != chosen_as_mol_num
                ):
                    reactant_order_list.append(0)
                    counter_of_matches = counter_of_matches + 1
                else:
                    reactant_order_list.append(0)
        return reactant_order_list

    def get_random_complimentary_mol(self, functional_group):
        """
        This function will get a dictionary of complimentary mols

        Inputs:
        :param str functional_group: the functional group of the needed
            complimentary molecule for the reaction

        Returns:
        :returns: list random_comp_mol: list with the SMILES string and name
            of molecule for the randomly chosen comp mol
        """
        infile = self.complimentary_mol_dict[functional_group]

        with open(infile, "r") as f:
            random_comp_mol_line = random.choice(f.readlines())
            random_comp_mol_line = (
                random_comp_mol_line.replace("\n", "")
                .replace("\t", " ")
                .replace("    ", " ")
            )
            for i in range(10):
                random_comp_mol_line.replace("  ", " ")
            parts = random_comp_mol_line.split(
                " "
            )  # split line into parts seperated by 4-spaces
            # parts = [x for x in random_comp_mol_line.split(" ") if x!= ""]
            # # split line into parts seperated by 4-spaces

            smile_list = parts[0]
            zinc_name_list = parts[1]
            random_comp_mol = [smile_list, zinc_name_list]

        return random_comp_mol

    def determine_functional_groups_in_mol(self, mol_deprotanated, mol_reprotanated):
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
            else:
                if mol_deprotanated.HasSubstructMatch(substructure):
                    list_subs_within_mol.append(key)
                else:
                    continue
        return list_subs_within_mol

    def run_smiles_click(self, ligand_smiles_string):
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
            the complimentary mol (None if it was a single reactant reaction)
            [reaction_product_smilestring, reaction_id_number,
            zinc_database_comp_mol_name]. returns None if all reactions failed or
            input failed to convert to a sanitizable rdkit mol.
        """
        try:
            mol = Chem.MolFromSmiles(
                ligand_smiles_string, sanitize=False
            )  # This is the input molecule which serves as the parent molecule
        except:
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
            print(
                "{} had no functional groups to react with.".format(
                    ligand_smiles_string
                )
            )
            return None

        shuffled_reaction_list = self.rand_key_list(
            self.reaction_dict
        )  # Randomize the order of the list of reactions

        tries = 0
        is_rxn_complete = False
        # go through all possible rxns in dicitonary of rxns using the random
        # order of rxns loop ends when a rxn is successful or when it runs out
        # of reactions
        while tries < len(shuffled_reaction_list) and is_rxn_complete is False:
            reaction_name = shuffled_reaction_list[tries]
            a_reaction_dict = self.reaction_dict[reaction_name]

            fun_groups_in_rxn = a_reaction_dict["functional_groups"]
            contains_group = None
            for i in range(0, len(fun_groups_in_rxn)):
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
                tries = tries + 1
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
                            reaction_products_list in [(), []]
                            or len(reaction_products_list) == 0
                    ):
                        # if reaction fails then lets move on to the next
                        # reaction
                        tries = tries + 1
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
                        # else:
                        tries = tries + 1

                except:
                    # if reaction fails then lets move on to the next reaction
                    mol_to_use = None
                    tries = tries + 1
                    break
            else:
                # for each functional group in the reaction, test
                # if the ligand has that as a substructure

                list_reactant_mols = []
                comp_mol_id = []
                for i in range(0, len(fun_groups_in_rxn)):
                    if i == contains_group:
                        # This is where the molecule goes
                        list_reactant_mols.append(mol_to_use)

                    else:
                        # for reactants which need to be taken from the
                        # complimentary dictionary. Find the reactants
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
                        for find_mol_tries in range(0, 100):

                            # find that group in the complimentary dictionary.
                            # comp_molecule = ["cccc", "ZINC123"]
                            comp_molecule = self.get_random_complimentary_mol(
                                functional_group_name
                            )

                            # zinc_database name
                            zinc_database_comp_mol_name = comp_molecule[1]

                            # Smiles String of complimentary molecule
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

                except:
                    reaction_product = None
                    tries = tries + 1
                    continue

                if (
                        reaction_products_list in [(), []]
                        or len(reaction_products_list) == 0
                ):
                    reaction_id_number = a_reaction_dict["RXN_NUM"]
                    tries = tries + 1
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
                    if reaction_product_smilestring is not None and is_rxn_complete is True:
                        # REACTION WORKED!
                        break
                    # try again
                    tries = tries + 1

        # end of the big while loop (while tries < len(shuffled_reaction_list)
        # and is_rxn_complete is False)

        # check that a reaction was sucessful
        if is_rxn_complete is True:
            reaction_product = MOH.check_sanitization(reaction_product)
            if reaction_product is None:
                return None

            reaction_product_smilestring = Chem.MolToSmiles(
                reaction_product, isomericSmiles=True
            )
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
            product_info = [
                reaction_product_smilestring,
                reaction_id_number,
                zinc_database_comp_mol_names,
            ]
            return product_info
        # reaction failed
        return None

    def check_if_product_is_good(self, reaction_product):
        """
        This function will test whether the product passes all of the
            requirements:
            1) Mol sanitizes
            2) It isn't in the self.list_of_already_made_smiles
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

        # passes
        return reaction_product_smilestring
#
def get_usable_fomat(infile):
    """
    This code takes a string for an file which is formated as an .smi file. It
    opens the file and reads in the components into a usable list.

    The .smi must follow the following format for each line:
        MANDATORY INFO
            part 1 is the SMILES string
            part 2 is the SMILES name/ID

        Optional info
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating
                - For default setting this is the Docking score
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)

    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_list_of_smiles: list of SMILES and their associated
        information formated into a list which is usable by the rest of Autogrow
    """

    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_list_of_smiles = []

    if os.path.exists(infile) is False:
        print("\nFile of Source compounds does not exist: {}\n".format(infile))
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts seperated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts seperated by 4-spaces

            choice_list = []
            for i in range(0, len(parts)):
                choice_list.append(parts[i])
            usable_list_of_smiles.append(choice_list)

    return usable_list_of_smiles
#
def react_with_multiple_reactants(mol_set, mol_name, rxn_obj):
    """
    This will run a single molecule through a 1-reactant reaction.

    If it fails it will return the name of mol (mol_info[1])
    If it passes it will return None
    Inputs:
    :param list mol_info: list of mol info
        mol_info[0] is the SMILES, 
        mol_info[1] is the name, 
        mol_info[-1] is the rdkit mol obj,
    :param rdkit.Chem.rdChemReactions.ChemicalReaction rxn_obj: the reaction object to use

    Returns:
    :returns: str mol_name: returns the mol_name if it fails to react;
        returns None if it passes reaction
    """
    # mol_name = mol_info[1]
    # mol_1 = mol_info[-1]
    try:
        # if reaction works keep it
        reaction_products_list = [
                        x[0] for x in rxn_obj.RunReactants(mol_set)
                    ]
    except:
        # return mol_name
        print("FAIL")
    sys.exit(0)

    if len(reaction_products_list) == 0:
        return mol_name
    # created a new compound so it passes
    return None

# 
def run_a_single_reactant_reaction(mol_info, rxn_obj):
    """
    This will run a single molecule through a 1-reactant reaction.

    If it fails it will return the name of mol (mol_info[1])
    If it passes it will return None
    Inputs:
    :param list mol_info: list of mol info
        mol_info[0] is the SMILES, 
        mol_info[1] is the name, 
        mol_info[-1] is the rdkit mol obj,
    :param rdkit.Chem.rdChemReactions.ChemicalReaction rxn_obj: the reaction object to use

    Returns:
    :returns: str mol_name: returns the mol_name if it fails to react;
        returns None if it passes reaction
    """
    mol_name = mol_info[1]
    mol_1 = mol_info[-1]
    
    try:
        # if reaction works keep it
        reaction_products_list = rxn_obj.RunReactants((mol_1,))
    except:
        return mol_name
    if len(reaction_products_list) == 0:
        return mol_name
    # created a new compound so it passes
    return None

def run_main(vars):
    """
    This runs the main testing. 
    
    Inputs:
    :param dict vars: Dictionary of User variables
    """

    # Force rxn_library to be custom because why else run this
    rxn_library = "Custom"

    rxn_library_file = vars["rxn_library_file"]
    function_group_library = vars["function_group_library"]
    complimentary_mol_dir = vars["complimentary_mol_directory"]
    number_of_processors = vars["number_of_processors"]

    rxn_library_variables = [
        rxn_library, 
        rxn_library_file,
        function_group_library,
        complimentary_mol_dir
    ]
    new_mutation_smiles_list = []

    a_smiles_click_chem_object = SmilesClickChem(
        rxn_library_variables, new_mutation_smiles_list
    )

    complimentary_mol_dict = a_smiles_click_chem_object.complimentary_mol_dict
    list_of_reaction_names = a_smiles_click_chem_object.list_of_reaction_names
    functional_group_dict = a_smiles_click_chem_object.functional_group_dict
    reaction_dict = a_smiles_click_chem_object.reaction_dict

    rxns_by_fun_group = {}
    for fun_group in functional_group_dict.keys():
        rxns_by_fun_group[fun_group] = []
    for rxn_name in list_of_reaction_names:
        current_rxn_dict = reaction_dict[rxn_name]
        for fun_group in current_rxn_dict["functional_groups"]:
            temp_list = rxns_by_fun_group[fun_group]
            temp_list.append(rxn_name)
            rxns_by_fun_group[fun_group] = temp_list

    failed_to_sanitize_by_fun_group = {}
    failed_to_react_by_fun_group = {}
    for fun_group in rxns_by_fun_group.keys():
        smi_comp_file = complimentary_mol_dict[fun_group]
        fun_group_list = get_usable_fomat(smi_comp_file)
        fun_group_mol_list = []
        failed_to_sanitize = []
        for info in fun_group_list:
            mol = Chem.MolFromSmiles(info[0])
            mol = MOH.check_sanitization(mol)
            if mol is None:
                failed_to_sanitize.append(info)
                continue
            info.append(mol)
            fun_group_mol_list.append(info)
        failed_to_sanitize_by_fun_group[fun_group] = failed_to_sanitize

        failed_to_react = []
        for rxn_name in rxns_by_fun_group[fun_group]:
            
            current_rxn_dict = reaction_dict[rxn_name]
            # Test example reactants
            example_smiles_rxn_reactants = current_rxn_dict["example_rxn_reactants"]
            example_smiles_rxn_reactants = example_smiles_rxn_reactants.replace("[","").replace("]","")
            example_smiles_rxn_reactants = example_smiles_rxn_reactants.replace(" ","").replace('"',"")
            example_smiles_rxn_reactants = example_smiles_rxn_reactants.split(",")
            
            example_rxn_reactants = []
            for smile_str in example_smiles_rxn_reactants:
                example_smile = Chem.MolFromSmiles(smile_str)
                
                example_mol = MOH.check_sanitization(mol)
                if example_mol is None:
                    printout = "example mol from rxn: {}".format(rxn_name)
                    printout = printout + " failed to sanitize in RDKit"
                    print(printout)
                    raise Exception(printout)
                example_rxn_reactants.append(example_mol)
                
            # convert example_rxn_reactants to tuple
            tuple_example_rxn_reactants = tuple(example_rxn_reactants)
            reaction_string = current_rxn_dict["reaction_string"]
            try:
                rxn_obj = AllChem.ReactionFromSmarts(reaction_string)
                rxn_obj.Initialize()
            except:
                printout = "rxn {} failed to be created.".format(rxn_name)
                printout = printout + "Rxn SMART is flawed"
                print(printout)
                raise Exception(printout)

            # Demo on example reactants
            print(rxn_name)
            print(reaction_string)
            print(tuple_example_rxn_reactants)
            reaction_products_list =  rxn_obj.RunReactants(tuple_example_rxn_reactants)
            print(reaction_products_list)
            sys.exit(0)
            example_results = react_with_multiple_reactants(example_rxn_reactants, "test_reactions", rxn_obj)
            print(example_results)
            if example_results is not None:
                printout = "rxn {} failed to run on example compounds.".format(rxn_name)
                printout = printout + "\nPlease check example compounds"
                print(printout)
                raise Exception(printout)

            num_reactants = current_rxn_dict["num_reactants"]
            list_of_reactants = []
            if num_reactants == 1:

                for mol_info in fun_group_mol_list:
                    list_of_reactants.append([mol_info, rxn_obj])

                list_of_reactants = [tuple(x) for x in list_of_reactants]
                list_of_reactants = tuple(list_of_reactants)

                output = mp.multi_threading(list_of_reactants, number_of_processors,
                                run_a_single_reactant_reaction)
                output = [x for x in output if x is not None]
                failed_to_react.append([rxn_name, output])
            else:
                functional_groups_rxn = current_rxn_dict["functional_groups"] 
                
                for f_group in functional_groups_rxn:
                    num_reactants = current_rxn_dict["num_reactants"]

def get_arguments_from_argparse(args_dict):
    """
    This function handles the arg parser arguments for the script.

    Inputs:
    :param dict args_dict: dictionary of parameters
    Returns:
    :returns: dict args_dict: dictionary of parameters
    """
    # Argument handling
    if args_dict["rxn_library_file"] == "" or args_dict["function_group_library"] == "":
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
                 THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
            )
    if os.path.exists(args_dict["rxn_library_file"]) is False:
        raise ValueError(
            "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
            THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
        )

    if args_dict["complimentary_mol_directory"] == "":
        raise ValueError(
            "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
            TO THE REACTION LIBRARY USING INPUT PARAMETER function_group_library"
        )
    else:
        if os.path.isdir(args_dict["complimentary_mol_directory"]) is False:
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
                TO THE REACTION LIBRARY USING INPUT PARAMETER complimentary_mol_directory"
            )

    if "number_of_processors" not in args_dict.keys():
        args_dict["number_of_processors"] = -1
    try:
        args_dict["number_of_processors"] = int(args_dict["number_of_processors"])
    except:
        raise ValueError(
                "number_of_processors must be an int. To use all processors set to -1."
            )

    return args_dict
#


# Argment parsing
PARSER = argparse.ArgumentParser()
# Mutation Settings
PARSER.add_argument(
    "--rxn_library_file",
    type=str,
    default="",
    required=True,
    help="This PATH to a Custom json file of SMARTS reactions to use for Mutation."
)
PARSER.add_argument(
    "--function_group_library",
    type=str,
    default="",
    required=True,
    help="This PATH for a dictionary of functional groups to be used for Mutation.",
)
PARSER.add_argument(
    "--complimentary_mol_directory",
    type=str,
    default="",
    required=True,
    help="This PATH to the directory containing all the molecules being used \
    to react with. The directory should contain .smi files contain SMILES of \
    molecules containing the functional group represented by that file. Each file \
    should be named with the same title as the functional groups described in \
    rxn_library_file & function_group_library +.smi \
    All Functional groups specified function_group_library must have its \
    own .smi file. We recommend you filter these dictionaries prior to Autogrow \
    for the Drug-likeliness and size filters you will Run Autogrow with.",
)
# processors and multithread mode
PARSER.add_argument(
    "--number_of_processors",
    "-p",
    type=int,
    default=-1,
    help="Number of processors to use for parallel calculations. \
    Set to -1 for all availble CPUs.",
)



ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)
run_main(ARGS_DICT)
print("done")