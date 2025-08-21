"""
This script will test a complementary molecule library to ensure all compounds
react in all reactions they may be used in.

Example submit:

python autogrow/accessory_scripts/test_complementary_mol_library.py \
--rxn_library_path autogrow4/plugins/mutation/reaction_libraries/all_rxns \
--output_folder autogrow/accessory_scripts/output/
"""

import __future__

import os
import json
import copy
import argparse
import contextlib
from typing import Any, Dict

import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


import support_scripts.Multiprocess as mp
import support_scripts.mol_object_handling as MOH
import gzip
from autogrow.plugins.registry_base import plugin_managers
from autogrow.plugins.mutation import MutationPluginManager, MutationBase
from autogrow.plugins.mutation.utils import validate_rxn_library_path

def react_with_multiple_reactants(mol_tuple, mol_name, rxn_obj):
    """
    This will run a single molecule through a 1-reactant reaction.

    If it fails it will return the name of mol (mol_info[1])
    If it passes it will return None
    Inputs:
    :param tuple mol_tuple: a tuple of all mols to react
    :param str mol_name: name of the molecule being tested
    :param rdkit.Chem.rdChemReactions.ChemicalReaction rxn_obj: the reaction object to use

    Returns:
    :returns: str mol_name: returns the mol_name if it fails to react;
        returns None if it passes reaction
    """
    try:
        # if reaction works keep it
        reaction_products_list = [x[0] for x in rxn_obj.RunReactants(mol_tuple)]
    except Exception:
        return mol_name

    if not reaction_products_list:
        return mol_name
    # created a new compound so it passes
    return None


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
    except Exception:
        return mol_name
    if not reaction_products_list:
        return mol_name
    # created a new compound so it passes
    return None


def get_rxn_and_examples(current_rxn_dict):
    """
    get the example reaction molecules from current_rxn_dict, create the rxn_obj,
    and test examples in the rxn.

    Inputs:
    :param dict current_rxn_dict: a dictionary of information about a reaction

    Returns:
    :returns: tuple example_rxn_reactants: a tuple of rdkit
            mol objects that are example compounds
    :returns: rdkit.Chem.rdChemReactions.ChemicalReaction rxn_obj: the
        reaction object to use
    """
    rxn_name = current_rxn_dict["reaction_name"]
    # Test example reactants
    example_smiles_rxn_reactants = current_rxn_dict["example_rxn_reactants"]
    example_rxn_reactants = []
    for smiles_str in example_smiles_rxn_reactants:
        example_mol = Chem.MolFromSmiles(smiles_str)

        example_mol = MOH.check_sanitization(example_mol)
        if example_mol is None:
            print(smiles_str)
            printout = f"example mol from rxn: {rxn_name} failed to sanitize in RDKit"
            print(printout)
            raise Exception(printout)
        example_rxn_reactants.append(example_mol)

    # convert example_rxn_reactants to tuple
    example_rxn_reactants = tuple(example_rxn_reactants)
    reaction_string = current_rxn_dict["reaction_string"]
    try:
        rxn_obj = AllChem.ReactionFromSmarts(reaction_string)
        rxn_obj.Initialize()
    except Exception as e:
        printout = f"rxn {rxn_name} failed to be created. Rxn SMART is flawed"
        print(printout)
        raise Exception(printout) from e

    # Demo on example reactants
    example_results = react_with_multiple_reactants(
        example_rxn_reactants, "test_reactions", rxn_obj
    )
    if example_results is not None:
        printout = (
            f"rxn {rxn_name} failed to run on example compounds."
            + "\nPlease check example compounds"
        )
        print(printout)
        raise Exception(printout)

    return example_rxn_reactants, rxn_obj


def run_all_for_fun_group(params, fun_group, rxns_by_fun_group, fragment_addition_object):
    """
    This runs the all testing for a single functional group.

    This will also write the compounds which pass to a .smi file.
    Args:
        params (Dict): Dictionary of User variables
        fun_group (str): functional group name
        rxns_by_fun_group (Dict): Dictionary of rxns names organized by
            functional groups
        fragment_addition_object (object): a FragmentAddition class object.
            This provides useful pathing information.
    Returns:
        list: a list of mol names which failed to react
        list: a list of mol names which failed to sanitize
    """
    # unpack variables
    complementary_mol_dict = fragment_addition_object.complementary_mol_dict
    reaction_dict = fragment_addition_object.reaction_dict
    number_of_processors = params["number_of_processors"]
    output_folder = params["output_folder"]

    fun_group_list = complementary_mol_dict.get(fun_group, [])

    fun_group_mol_list = []
    failed_to_sanitize = []
    for info in fun_group_list:
        mol = Chem.MolFromSmiles(info[0])
        mol = MOH.check_sanitization(mol)
        if mol is None:
            failed_to_sanitize.append(info)
            continue
        temp = copy.deepcopy(info)
        temp.append(mol)
        fun_group_mol_list.append(temp)

    # print info about failures
    if failed_to_sanitize:
        printout = (
            f"{len(failed_to_sanitize)} compounds "
            + f"failed to sanitize from: {fun_group}"
        )
        print(printout)

    failed_to_react = []
    for rxn_name in rxns_by_fun_group[fun_group]:

        current_rxn_dict = reaction_dict[rxn_name]
        example_reactants, rxn_obj = get_rxn_and_examples(current_rxn_dict)

        list_of_reactants = []
        functional_groups_rxn = current_rxn_dict["functional_groups"]
        i_count_to_use = None
        for i_count in range(len(functional_groups_rxn)):
            f_group = functional_groups_rxn[i_count]

            if fun_group == f_group:
                i_count_to_use = i_count
            else:
                continue
        if i_count_to_use is None:
            raise Exception("This is a code error.")

        list_of_reactants = []
        for mol_info in fun_group_mol_list:
            mol_tuple_temp = []
            for i_count in range(len(functional_groups_rxn)):
                if i_count == i_count_to_use:
                    mol_tuple_temp.append(mol_info[-1])
                else:
                    mol_tuple_temp.append(example_reactants[i_count])

            list_of_reactants.append((tuple(mol_tuple_temp), mol_info[1], rxn_obj))

        output = mp.multi_threading(
            list_of_reactants, number_of_processors, react_with_multiple_reactants
        )
        output = [x for x in output if x is not None]
        failed_to_react.append([rxn_name, output])

        # print info about failures
        if output:
            printout = "{} compounds failed to react from ".format(len(output))
            printout += "react from {} ".format(fun_group)
            printout += "in rxn: {}".format(rxn_name)
            print(printout)

    master_failed_to_react = []
    master_passes_reactions = []
    for fail_mol_list in failed_to_react:
        master_failed_to_react.extend(fail_mol_list[1])
    for mol_info in fun_group_list:
        if mol_info[1] in master_failed_to_react:
            continue
        master_passes_reactions.append(" ".join(mol_info))
    # write to output .smi file
    with open(output_folder + fun_group + ".smi", "w") as f:
        f.write("\n".join(master_passes_reactions))
    return failed_to_react, failed_to_sanitize


def run_main(params: Dict[str, Any]):
    """
    This runs the main testing.
    Args:
        params (Dict[str, Any]): Dictionary of User variables
    """
    # Setup the plugin managers
    plugin_managers.setup_plugin_managers(params)
    # Instantiate and setup mutation manager to get FragmentAddition plugin
    mutation_manager = MutationPluginManager(MutationBase)
    mutation_manager.setup_plugin_manager(params, plugin_managers)
    fragment_addition_object = mutation_manager.plugins.get('FragmentAddition')
    if fragment_addition_object is None:
        raise ValueError("FragmentAddition plugin could not be initialized. "
                         "Ensure it is enabled and configured correctly.")

    list_of_reaction_names = list(fragment_addition_object.reaction_dict.keys())
    functional_group_dict = fragment_addition_object.functional_group_dict
    reaction_dict = fragment_addition_object.reaction_dict
    rxns_by_fun_group = {fun_group: [] for fun_group in functional_group_dict.keys()}
    for rxn_name in list_of_reaction_names:
        current_rxn_dict = reaction_dict[rxn_name]
        for fun_group in current_rxn_dict["functional_groups"]:
            temp_list = rxns_by_fun_group[fun_group]
            temp_list.append(rxn_name)
            rxns_by_fun_group[fun_group] = temp_list

    failed_to_sanitize_by_fun_group = {}
    failed_to_react_by_fun_group = {}

    for fun_group in rxns_by_fun_group:
        failed_to_react, failed_to_sanitize = run_all_for_fun_group(
            params, fun_group, rxns_by_fun_group, fragment_addition_object
        )
        failed_to_react_by_fun_group[fun_group] = failed_to_react
        failed_to_sanitize_by_fun_group[fun_group] = failed_to_sanitize

    # Handle saving log
    output_folder = params["output_folder"]
    with open(f"{output_folder}failed_to_sanitize_mol_by_fun_group.json", "w") as fp:
        json.dump(failed_to_sanitize_by_fun_group, fp, indent=4)

    with open(f"{output_folder}failed_to_react_by_fun_group.json", "w") as fp:
        json.dump(failed_to_react_by_fun_group, fp, indent=4)

    master_failed_list = []
    for fun_group in failed_to_react_by_fun_group:
        temp = [x[1] for x in failed_to_react_by_fun_group[fun_group]]
        for x in temp:
            master_failed_list.extend(x)
    master_failed_list = list(set(master_failed_list))
    if not master_failed_list:
        print("All compounds passed!")
    else:
        print(f"{len(master_failed_list)} compounds failed. Please check logs")


def get_arguments_from_argparse(args_dict):
    """
    This function handles the arg parser arguments for the script.
    Args:
        args_dict (dict): dictionary of parameters
    Returns:
        dict: dictionary of parameters
    """
    validate_rxn_library_path(args_dict)

    if "number_of_processors" not in args_dict.keys():
        args_dict["number_of_processors"] = -1
    try:
        args_dict["number_of_processors"] = int(args_dict["number_of_processors"])
    except Exception as e:
        raise ValueError(
            "number_of_processors must be an int. \
   To use all processors set to -1."
        ) from e

    if "output_folder" not in args_dict.keys():
        printout = (
            "output_folder is a required variable. it is the PATH to where "
            + "filtered .smi file and log files will be placed. Will save a file "
            + "in this directory for mols which failed sanitization, mols which "
            + "failed to react in specific reactions, and .smi files that contain "
            + "all mols that reacted properly."
        )
        raise ValueError(printout)

    if type(args_dict["output_folder"]) != str or args_dict["output_folder"] == "":
        printout = (
            "output_folder is a required variable. it is the PATH to where "
            + "filtered .smi file and log files will be placed. Will save a file "
            + "in this directory for mols which failed sanitization, mols which "
            + "failed to react in specific reactions, and .smi files that contain "
            + "all mols that reacted properly."
        )
        raise ValueError(printout)

    args_dict["output_folder"] = os.path.abspath(args_dict["output_folder"]) + os.sep
    if os.path.exists(args_dict["output_folder"]) is True:
        if os.path.isdir(args_dict["output_folder"]) is False:
            print(args_dict["output_folder"])
            printout = "output_folder must be a directory. Please check input arguments"
            raise ValueError(printout)
    else:
        with contextlib.suppress(Exception):
            os.mkdir(args_dict["output_folder"])
        if os.path.exists(args_dict["output_folder"]) is False:
            raise Exception("output_folder could not be made or found.")

    # Setup params for plugin system
    args_dict['RDKitToolkit'] = True
    args_dict['FragmentAddition'] = True

    return args_dict


# Argument parsing
PARSER = argparse.ArgumentParser()
# Mutation Settings
PARSER.add_argument(
    "--rxn_library_path",
    type=str,
    default="all_rxns",
    required=True,
    help="Path to a reaction library directory, or a built-in library name (e.g., all_rxns).",
)
PARSER.add_argument(
    "--output_folder",
    type=str,
    default="",
    required=True,
    help="This PATH to where filtered .smi file and log files will be placed. \
  Will save a file in this directory for mols which failed sanitization, \
  mols which failed to react in specific reactions, and .smi files \
  that contain all mols that reacted properly.",
)
# processors and multithread mode
PARSER.add_argument(
    "--number_of_processors",
    "-p",
    type=int,
    default=-1,
    help="Number of processors to use for parallel calculations. \
 Set to -1 for all available CPUs.",
)


ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)
run_main(ARGS_DICT)
print("done")
