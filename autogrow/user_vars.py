"""user_vars
This should contain the functions for defining input variables.
Both the default variables and the user input variables.
This should also validate them.
"""


import __future__

import contextlib
import os
import copy
import datetime
import json
import sys
import platform
from shutil import copyfile
from typing import Any, Dict, List, Optional, Tuple, Union

from autogrow.config.config_run_directory import set_run_directory


def program_info() -> str:
    """
    Get the program version number, etc.

    Returns:
    :returns: str program_output: a string for the print of the program information
    """
    program_output = "\nAutoGrow Version 4.0.3\n" + " ================== \n"
    program_output += "If you use AutoGrow 4.0.3 in your research, please cite the following reference:\n"
    program_output += "Spiegel, J.O., Durrant, J.D. \n"
    program_output += "AutoGrow4: an open-source genetic algorithm "
    program_output += "for de novo drug design and lead optimization. \n"
    program_output += "J Cheminform 12, 25 (2020). \n"
    program_output += "[doi: 10.1186/s13321-020-00429-4]\n"
    program_output += " ================== \n\n"

    return program_output


#
def save_vars_as_json(params: Dict[str, Any]) -> None:
    """
    This function saves the params dictionary as a json file. This can be used
    later to track experiments and is necessary for several of the utility
    scripts.
    It saves all variables except the parallelizer class object.

    It saves the file to the output_directory + "vars.json"
        -If AutoGrow has been run multiple times for the same directory it
        will save the new vars file as append a number to the file name
        starting with 2. The util scripts will only look at the original "vars.json"
            ie) output_directory + "vars_2.json"

    Inputs:
    :param dict params: dict of user variables which will govern how the programs runs
    """
    output_directory = params["output_directory"]

    vars_file = output_directory + os.sep + "vars.json"
    if os.path.exists(vars_file):
        # vars.json already exist
        # lets make the next file
        path_exists = True
        i = 2
        while path_exists:
            vars_file = f"{output_directory}{os.sep}vars_{i}.json"
            if os.path.exists(vars_file):
                i = i + 1
            else:
                path_exists = False

    temp_vars = {k: copy.deepcopy(params[k]) for k in params if "parallelizer" not in k}
    with open(vars_file, "w") as fp:
        json.dump(temp_vars, fp, indent=4)


def multiprocess_handling(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    This function handles the multiprocessing functions. It establishes a Paralellizer object
    and adds it to the params dictionary.
    Inputs:
    :param dict params: dict of user variables which will govern how the programs runs
    Returns:
    :returns: dict params: dict of user variables which will govern how the programs runs
    """

    # Handle Serial overriding number_of_processors
    # serial fixes it to 1 processor
    if params["multithread_mode"].lower() == "serial":
        params["multithread_mode"] = "serial"
        if params["number_of_processors"] != 1:
            print(
                "Because --multithread_mode was set to serial, "
                + "this will be run on a single processor."
            )
        params["number_of_processors"] = 1

    # Avoid EOF error
    from autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import (
        Parallelizer,
    )

    params["parallelizer"] = Parallelizer(
        params["multithread_mode"], params["number_of_processors"], True
    )

    return params


def test_docking_executables(
    params: Dict[str, Any], vina_exe: str, qvina2_exe: str
) -> bool:
    """
    This will test if docking executables are compatible with OS.
    This is only required for MacOS.

    Test will output the version of Vina and QVina2.1 executables to txt file
    in the root_output_folder (docking_exe_MACOS_test.txt)
    If both executables are compatible with this MacOS there should be the following
    2 lines in the txt file:
        AutoDock Vina 1.1.2 (May 11, 2011)
        QuickVina 2.1 (24 Dec, 2017)

    Returns True if both work and returns False.

    Inputs:
    :param dict params: dict of user variables which will govern how the programs runs
    :param str vina_exe: path to vina executable
    :param str qvina2_exe: path to quick vina 2 executable
    Returns:
    :returns: bool bool: returns True if both docking executables work; False if either fails
    """
    test_vina_outfile = (
        params["root_output_folder"] + os.sep + "docking_exe_MACOS_test.txt"
    )
    try:
        command = "{} --version > {arg_2} 2>> {arg_2}".format(
            vina_exe, arg_2=test_vina_outfile
        )
        os.system(command)
        command = "{} --version >> {arg_2} 2>> {arg_2}".format(
            qvina2_exe, arg_2=test_vina_outfile
        )
        os.system(command)
    except Exception:
        printout = "Docking executables could not be found."
        # is not compatible on this OS. \nPlease use docker "
        return False

    with open(test_vina_outfile, "r") as test_file:
        lines = test_file.readlines()

    return True


def _process_mac_execs(arg0, vina_exe, qvina2_exe):
    command = f"{arg0}{vina_exe}"
    os.system(command)
    command = f"{arg0}{qvina2_exe}"
    os.system(command)


def _notify_docker_requirement(arg0):
    printout = f"{arg0}Please run using docker version of AutoGrow."
    print(printout)
    raise Exception(printout)


############################################
###### Variables Handlining Settings #######
############################################
def check_for_required_inputs(input_params):
    """
    Confirm all the required inputs were provided.

    Required Variables go here.

    Inputs:
    :param dict input_params: The parameters. A dictionary of {parameter name: value}.
    """
    keys_from_input = list(input_params.keys())

    list_of_required_inputs = [
        "filename_of_receptor",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "root_output_folder",
        "source_compound_file",
    ]

    missing_variables = []
    for variable in list_of_required_inputs:
        if variable in keys_from_input:
            continue
        missing_variables.append(variable)

    if missing_variables:
        _extracted_from_check_for_required_inputs_31(missing_variables)
    # Make sure the dimmensions are in floats. If in int convert to float.
    for x in ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]:
        if type(input_params[x]) == float:
            continue
        if type(input_params[x]) == int:
            input_params[x] = float(input_params[x])
        else:
            printout = f"\n{x} must be a float value.\n"
            print(printout)
            raise Exception(printout)

    # Check Docking Exhaustiveness and modes...
    if "docking_exhaustiveness" in list(input_params.keys()):
        if input_params["docking_exhaustiveness"] == "None":
            input_params["docking_exhaustiveness"] = None
        if input_params["docking_exhaustiveness"] is not None:

            with contextlib.suppress(Exception):
                input_params["docking_exhaustiveness"] = int(
                    input_params["docking_exhaustiveness"]
                )
            if type(input_params["docking_exhaustiveness"]) not in [
                int,
                float,
            ]:
                raise Exception(
                    "docking_exhaustiveness needs to be an interger. \
                    If you do not know what to use, leave this blank and the \
                    default for the docking software will be used."
                )
    if "docking_num_modes" in list(input_params.keys()):
        if input_params["docking_num_modes"] == "None":
            input_params["docking_num_modes"] = None
        if input_params["docking_num_modes"] is not None:
            with contextlib.suppress(Exception):
                input_params["docking_num_modes"] = int(
                    input_params["docking_num_modes"]
                )
            if type(input_params["docking_num_modes"]) not in [int, float]:
                raise Exception(
                    "docking_num_modes needs to be an interger. \
                    If you do not know what to use, leave this blank and the \
                    default for the docking software will be used."
                )

    # Check numbers which may be defined by first generation
    if "top_mols_to_seed_next_generation_first_generation" not in list(
        input_params.keys()
    ):
        if "top_mols_to_seed_next_generation" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["top_mols_to_seed_next_generation"] = 10
            input_params["top_mols_to_seed_next_generation_first_generation"] = 10
        else:
            input_params[
                "top_mols_to_seed_next_generation_first_generation"
            ] = input_params["top_mols_to_seed_next_generation"]

    if "number_of_crossovers_first_generation" not in list(input_params.keys()):
        if "number_of_crossovers" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_of_crossovers"] = 10
            input_params["number_of_crossovers_first_generation"] = 10
        else:
            input_params["number_of_crossovers_first_generation"] = input_params[
                "number_of_crossovers"
            ]

    if "number_of_mutants_first_generation" not in list(input_params.keys()):
        if "number_of_mutants" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_of_mutants"] = 10
            input_params["number_of_mutants_first_generation"] = 10
        else:
            input_params["number_of_mutants_first_generation"] = input_params[
                "number_of_mutants"
            ]

    if "number_elitism_advance_from_previous_gen_first_generation" not in list(
        input_params.keys()
    ):
        if "number_elitism_advance_from_previous_gen" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_elitism_advance_from_previous_gen"] = 10
            input_params[
                "number_elitism_advance_from_previous_gen_first_generation"
            ] = 10
        else:
            input_params[
                "number_elitism_advance_from_previous_gen_first_generation"
            ] = input_params["number_elitism_advance_from_previous_gen"]

    #######################################
    # Check that all required files exist #
    #######################################

    # convert paths to abspath, in case necessary
    input_params["filename_of_receptor"] = os.path.abspath(
        input_params["filename_of_receptor"]
    )
    input_params["root_output_folder"] = os.path.abspath(
        input_params["root_output_folder"]
    )
    input_params["source_compound_file"] = os.path.abspath(
        input_params["source_compound_file"]
    )

    # Check filename_of_receptor exists
    if os.path.isfile(input_params["filename_of_receptor"]) is False:
        raise NotImplementedError(
            "Receptor file can not be found. File must be a .PDB file."
        )
    if ".pdb" not in input_params["filename_of_receptor"]:
        raise NotImplementedError("filename_of_receptor must be a .PDB file.")

    # Check root_output_folder exists
    if os.path.exists(input_params["root_output_folder"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(input_params["root_output_folder"])
        except Exception as e:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            ) from e

        if os.path.exists(input_params["root_output_folder"]) is False:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

    if os.path.isdir(input_params["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )

    # Check source_compound_file exists
    if os.path.isfile(input_params["source_compound_file"]) is False:
        raise NotImplementedError(
            "source_compound_file can not be found. \
            File must be a tab delineated .smi file."
        )
    if ".smi" not in input_params["source_compound_file"]:
        raise NotImplementedError(
            "source_compound_file must be a \
            tab delineated .smi file."
        )


# TODO Rename this here and in `check_for_required_inputs`
def _extracted_from_check_for_required_inputs_31(missing_variables):
    printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunAutogrow -h"
    printout += "\nThe following required variables are missing: "
    for variable in missing_variables:
        printout = printout + "\n\t" + variable
    print("")
    print(printout)
    print("")
    raise NotImplementedError("\n" + printout + "\n")


def determine_bash_timeout_vs_gtimeout() -> str:
    """
    This function tests whether we should use the BASH command "timeout" (for linux)
     or the coreutils function "gtimeout" for MacOS which can be obtained
     through homebrew

    Returns:
    :returns: str timeout_option: A string either "timeout" or "gtimeout" describing
     whether the bash terminal is able to use the bash function timeout or gtimeout
    """

    if sys.platform.lower() in ["linux", "linux2"]:
        # Should be true and default installed in all Linux machines
        return "timeout"

    command = 'timeout 1 echo " "'
    # Running the os.system command for command will return 0,1, or 32512
    # 0 means that the timeout function works (most likely this is a linux os)
    # 32512 means that the timeout function DOES NOT Work (most likely this is MacOS)

    try:  # timeout or gtimeout
        timeout_result = os.system(f"g{command}")
    except Exception as e:
        raise Exception(
            "Something is very wrong. This OS may not be supported \
            by Autogrow or you may need to execute through Bash."
        ) from e
    if timeout_result == 0:
        return "gtimeout"
    print("gtimeout failed to run, we will check timeout")

    try:  # timeout or gtimeout
        timeout_result = os.system(command)
    except Exception as e:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
            Autogrow or you may need to execute through Bash."
        ) from e

    if timeout_result == 0:
        return "timeout"
    printout = (
        "Need to install GNU tools for Bash to work. \n"
        + "This is essential to use Bash Timeout function in Autogrow. \n"
    )
    printout += "\t This will require 1st installing homebrew. \n"
    printout += "\t\t Instructions found at: https://brew.sh/ \n"
    printout += "\t Once brew is installed, please run:"
    printout += " sudo brew install coreutils \n\n"
    print(printout)
    raise Exception(printout)


def check_dependencies() -> None:
    """
    This function will try to import all the installed dependencies that will be
    used in Autogrow. If it fails to import it will raise an ImportError
    """

    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option not in ["timeout", "gtimeout"]:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
        Autogrow or you may need to execute through Bash."
        )

    try:
        import rdkit  # type: ignore
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem  # type: ignore
        from rdkit.Chem import rdDepictor  # type: ignore
        from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore
        from rdkit.Chem.Draw import PrepareMolForDrawing  # type: ignore
        from rdkit.Chem import rdFMCS  # type: ignore
        from rdkit.Chem import FilterCatalog  # type: ignore
        from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
        import rdkit.Chem.Lipinski as Lipinski  # type: ignore
        import rdkit.Chem.Crippen as Crippen  # type: ignore
        import rdkit.Chem.Descriptors as Descriptors  # type: ignore
        import rdkit.Chem.MolSurf as MolSurf  # type: ignore

    except Exception as e:
        print("You need to install rdkit and its dependencies.")
        raise ImportError("You need to install rdkit and its dependencies.") from e

    # molvs is prepackaged within gypsum_dl
    # try:
    #     from molvs import standardize_smiles as ssmiles
    # except Exception:
    #     print("You need to install molvs and its dependencies.")
    #     raise ImportError("You need to install molvs and its dependencies.")

    try:
        import numpy
    except Exception as exc:
        print("You need to install numpy and its dependencies.")
        raise ImportError("You need to install numpy and its dependencies.") from exc

    try:
        from scipy.cluster.vq import kmeans2  # type: ignore
    except Exception as err:
        print("You need to install scipy and its dependencies.")
        raise ImportError("You need to install scipy and its dependencies.") from err

    try:
        import os
        import sys
        import glob
        import subprocess
        import multiprocessing
        import time
    except Exception as exception:
        print(
            "Missing a Python Dependency. Could be import: os,sys,glob,\
                subprocess,multiprocess, time."
        )
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
                os,sys,glob,subprocess,multiprocess, time."
        ) from exception

    try:
        import copy
        import random
        import string
        import math
    except Exception as e:
        print("Missing a Python Dependency. Could be import: copy,random, string,math")
        raise ImportError(
            "Missing a Python Dependency. Could be import: copy,random, string,math"
        ) from e

    try:
        from collections import OrderedDict
        import webbrowser
        import argparse
        import itertools
        import unittest
    except Exception as exc:
        print(
            "Missing a Python Dependency. Could be import: collections,\
            webbrowser,argparse,itertools,unittest"
        )
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
            collections,webbrowser,argparse,itertools,unittest"
        ) from exc

    try:
        import textwrap
        import pickle
        import json
    except Exception as error:
        print("Missing a Python Dependency. Could be import: textwrap, pickle,json")
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
            textwrap, pickle,json"
        ) from error


def define_defaults() -> Dict[str, Any]:
    """
    Sets the command-line parameters to their default values.

    Returns:
    :returns: dict params: a dictionary of all default variables
    """

    # where we are currently (absolute filepath from route)
    # used for relative pathings
    script_dir = os.path.dirname(os.path.realpath(__file__))

    params: Dict[str, Any] = {}

    #### OPTIONAL FILE-LOCATION VARIABLES ####
    # (RECOMMEND SETTING TO "" SO AUTOGROW CAN AUTOLOCATE THESE FILES)#

    # PARSER.add_argument('--conversion_choice', choices
    #    = ["MGLTools","obabel"], default="MGLTools",
    params["conversion_choice"] = "obabel"
    params["obabel_path"] = "obabel"

    # Crossover function
    params["start_a_new_run"] = False
    params["max_time_mcs_prescreen"] = 1
    params["max_time_mcs_thorough"] = 1
    params["min_atom_match_mcs"] = 4
    params["protanate_step"] = False

    # Mutation Settings
    params["rxn_library"] = "click_chem_rxns"
    params["rxn_library_file"] = ""
    params["function_group_library"] = ""
    params["complementary_mol_directory"] = ""

    # processors
    params["number_of_processors"] = 1
    params["multithread_mode"] = "multithreading"

    # Genetic Algorithm Components
    # params["selector_choice"] = "Roulette_Selector"
    # params["tourn_size"] = 0.1

    # Seeding next gen and diversity
    params["top_mols_to_seed_next_generation_first_generation"] = 10
    params["top_mols_to_seed_next_generation"] = 10
    params["diversity_mols_to_seed_first_generation"] = 10
    params["diversity_seed_depreciation_per_gen"] = 2

    # Populations settings
    params["filter_source_compounds"] = True
    params["dock_source_compounds_first"] = True
    params["num_generations"] = 10
    params["number_of_crossovers_first_generation"] = 10
    params["number_of_mutants_first_generation"] = 10
    params["number_of_crossovers"] = 10
    params["number_of_mutants"] = 10
    params["number_elitism_advance_from_previous_gen"] = 10
    params["number_elitism_advance_from_previous_gen_first_generation"] = 10
    params["redock_elite_from_previous_gen"] = False

    # Filters
    params["LipinskiStrictFilter"] = False
    params["LipinskiLenientFilter"] = False
    params["GhoseFilter"] = False
    params["GhoseModifiedFilter"] = False
    params["MozziconacciFilter"] = False
    params["VandeWaterbeemdFilter"] = False
    params["PAINSFilter"] = False
    params["NIHFilter"] = False
    params["BRENKFilter"] = False
    params["No_Filters"] = False
    params["alternative_filter"] = None

    # docking
    # params["dock_choice"] = "QuickVina2Docking"
    # params["vina_like_executable"] = None
    # params["docking_exhaustiveness"] = None
    # params["docking_num_modes"] = None
    # params["docking_timeout_limit"] = 120

    # scoring
    params["scoring_choice"] = "VINA"
    params["rescore_lig_efficiency"] = False

    # gypsum # max variance is the number of conformers made per ligand
    params["max_variants_per_compound"] = 3
    params["gypsum_thoroughness"] = 3
    params["min_ph"] = 6.4
    params["max_ph"] = 8.4
    params["pka_precision"] = 1.0
    params["gypsum_timeout_limit"] = 10

    # Other params
    params["debug_mode"] = False
    params["generate_plot"] = True
    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option in ["timeout", "gtimeout"]:
        params["timeout_vs_gtimeout"] = timeout_option
    else:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
             Autogrow or you may need to execute through Bash."
        )

    return params


############################################
######## Input Handlining Settings #########
############################################
def convert_json_params_from_unicode(params_unicode: Dict[str, Any]) -> Dict[str, Any]:
    """
    Set the parameters that will control this ConfGenerator object.

    :param dict params_unicode: The parameters. A dictionary of {parameter name:
                value}.
    Returns:
    :returns: dict params: Dictionary of User variables
    """
    # Also, rdkit doesn't play nice with unicode, so convert to ascii

    # Because Python2 & Python3 use different string objects, we separate their
    # usecases here.
    params = {}
    if sys.version_info < (3,):
        for param in params_unicode:
            val = params_unicode[param]
            if isinstance(val, unicode):
                val = str(val).encode("utf8")
            key = param.encode("utf8")
            params[key] = val
    else:
        for param in params_unicode:
            val = params_unicode[param]
            key = param
            params[key] = val
    return params


def check_value_types(
    params: Dict[str, Any], argv: Dict[str, Any]
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    This checks that all the user variables loaded in use that same or comparable
    datatypes as the defaults in params. This prevents type issues later in the
    simulation.

    Given the many uservars and the possibility for intentional differences,
    especially as the program is developed, this function tries to be
    NOT OPINIONATED, only correcting for several obvious and easy to correct issues
    of type discrepancies occur between argv[key] and params[key]
        ie
            1) argv[key] = "true" and params[key] = False
                this script will not change argv[key] to False... it will
                convert "true" to True
                ---> argv[key]=True
            2) argv[key] = "1.01" and params[key] = 2.1
                this script will change argv[key] from "1.01" to float(1.01)

    Inputs:
    :param dict params: Dictionary of program defaults, which will later be
        overwritten by argv values
    :param dict argv: Dictionary of User specified variables
    Returns:
    :returns: dict params: Dictionary of program defaults, which will later
        be overwritten by argv values
    :returns: dict argv: Dictionary of User specified variables
    """
    for key in list(argv.keys()):
        if key not in list(params.keys()):
            # Examples may be things like filename_of_receptor or
            # dimensions of the docking box
            #   Just skip these
            continue

        if type(argv[key]) != type(params[key]):
            # Several variable default is None which means checks are
            # processed elsewhere...
            if params[key] is None:
                if type(argv[key]) != str:
                    continue
                # check argv[key] is "none" or "None"
                if argv[key].lower() == "none":
                    argv[key] = None
            elif type(params[key]) in [int, float]:
                if type(argv[key]) in [int, float]:
                    # this is fine
                    continue
                elif type(argv[key]) == str:
                    try:
                        temp_item = float(argv[key])
                        if type(temp_item) == float:
                            argv[key] = temp_item
                        else:
                            printout = (
                                "This parameter is the wrong type.\n \t Check : "
                                + f"{key} type={type(argv[key])}\n"
                            )
                            printout += f"\t Should be type={type(params[key])}\n\t"
                            printout += "Please check Autogrow documentation using -h"
                            raise IOError(printout)
                    except Exception as e:
                        printout = (
                            "This parameter is the wrong type. \n \t Check :"
                            + f" {key} type={type(argv[key])}\n"
                        )
                        printout += f"\t Should be type={type(params[key])}\n\t"
                        printout += "Please check Autogrow documentation using -h"
                        raise IOError(printout) from e
                else:
                    printout = (
                        "This parameter is the wrong type. \n \t Check :"
                        + f" {key} type={type(argv[key])}\n"
                    )
                    printout += f"\t Should be type={type(params[key])}\n\t"
                    printout += "Please check Autogrow documentation using -h"
                    raise IOError(printout)
            elif type(params[key]) == bool:
                if argv[key] is None:
                    # Do not try to handle this. May make sense.
                    continue
                if type(argv[key]) == str:
                    if argv[key].lower() in ["true", "1"]:
                        argv[key] = True
                    elif argv[key].lower() in ["false", "0"]:
                        argv[key] = False
                    elif argv[key].lower() in ["none"]:
                        argv[key] = None
                    else:
                        printout = (
                            "This parameter is the wrong type. \n \t Check :"
                            + f" {key} type={type(argv[key])}\n"
                        )
                        printout += f"\t Should be type={type(params[key])}\n\t"
                        printout += "Please check Autogrow documentation using -h"
                        raise IOError(printout)
                else:
                    printout = (
                        "This parameter is the wrong type. \n \t Check :"
                        + f" {key} type={type(argv[key])}\n"
                    )
                    printout += f"\t Should be type={type(params[key])}\n\t"
                    printout += "Please check Autogrow documentation using -h"
                    raise IOError(printout)
    return params, argv


def load_in_commandline_parameters(argv: Dict[str, Any]) -> Tuple[Dict[str, Any], str]:
    """
    Load in the command-line parameters

    Inputs:
    :param dict argv: Dictionary of User specified variables

    Returns:
    :returns: dict params: Dictionary of User variables
    :returns: str printout: a string to be printed to screen and saved to output file
    """

    params = define_defaults()

    # Load the parameters from the json
    if "json" in argv:
        json_vars = json.load(open(argv["json"]))
        json_vars = convert_json_params_from_unicode(json_vars)
        check_for_required_inputs(json_vars)
        params, json_vars = check_value_types(params, json_vars)
        for key in list(json_vars.keys()):
            params[key] = json_vars[key]

    else:
        check_for_required_inputs(argv)
        params, argv = check_value_types(params, argv)
        for key in list(argv.keys()):
            params[key] = argv[key]

    params = multiprocess_handling(params)

    printout = f"(RE)STARTING AUTOGROW 4.0: {str(datetime.datetime.now())}"
    printout = printout + program_info()
    printout = (
        printout + "\nUse the -h tag to get detailed help regarding program usage.\n"
    )
    print(printout)
    sys.stdout.flush()
    # Check all Dependencies are installed
    check_dependencies()

    params = filter_choice_handling(params)

    ###########################################
    ########## Check variables Exist ##########
    ###########################################

    # Check if the Operating System is Windows, if so turn off Multiprocessing.
    if os.name in {"nt", "ce"}:
        # so it's running under windows. multiprocessing disabled
        params["number_of_processors"] = 1
        printout = (
            printout
            + "\nWARNING: Multiprocessing is disabled on\
            windows machines.\n"
        )

    # make sure directories end in os.sep
    if params["root_output_folder"][-1] != os.sep:
        params["root_output_folder"] = params["root_output_folder"] + os.sep

    # More Handling for Windows OS
    # convert path names with spaces if this is windows
    if os.name in {"nt", "ce"}:
        # so it's running under windows. multiprocessing disabled

        if " " in params["filename_of_receptor"]:
            params["filename_of_receptor"] = '"' + params["filename_of_receptor"] + '"'
        if " " in params["root_output_folder"]:
            params["root_output_folder"] = '"' + params["root_output_folder"] + '"'

    # output the paramters used
    printout = printout + "\nPARAMETERS" + "\n"
    printout = f"{printout} ========== " + "\n"

    # Make sure scripts and executables exist
    # If MGLTools is being used handle its paths
    if params["conversion_choice"] == "MGLToolsConversion":
        if not os.path.exists(params["prepare_ligand4.py"]) and not os.path.exists(
            params["prepare_ligand4.py"].replace('"', "")
        ):
            printout = (
                printout
                + "\nERROR: Could not find prepare_ligand4.py at "
                + params["prepare_ligand4.py"]
                + "\n"
            )
            print(printout)
            raise NotImplementedError(printout)
        if not os.path.exists(params["prepare_receptor4.py"]) and not os.path.exists(
            params["prepare_receptor4.py"].replace('"', "")
        ):
            printout = (
                printout
                + "\nERROR: Could not find prepare_receptor4.py at "
                + params["prepare_receptor4.py"]
                + "\n"
            )
            print(printout)
            raise NotImplementedError(printout)

    if not os.path.exists(params["filename_of_receptor"]):
        printout = (
            printout
            + '\nERROR: There receptor file does not exist: "'
            + params["filename_of_receptor"]
            + '".'
            + "\n"
        )
        print(printout)
        raise NotImplementedError(printout)

    # Check if the user wants to continue a run or start a new run.
    # Make new run directory if necessary. return the Run folder path
    # The run folder path will be where we place our generations and output files
    params["output_directory"] = set_run_directory(
        params["root_output_folder"], params["start_a_new_run"]
    )

    # Save variables in params dict to a .json file for later usage and reference
    # It saves the file to the output_directory + "vars.json"
    # -If AutoGrow has been run multiple times for the same directory it
    # will save the new vars file as append a number to the file name
    # starting with 2. The util scripts will only look at the original "vars.json"
    #     ie) output_directory + "vars_2.json"
    save_vars_as_json(params)

    return params, printout


def make_complete_children_dict(purpose_of_object: str) -> Dict[str, Any]:
    """
    This will retrieve all the names of every child class of the parent class
    This can be either filter, parent_pdbqt_converter, ParentDocking,
    or ParentScoring

    Inputs:
    :param str purpose_of_object: either filter, parent_pdbqt_converter,
        ParentDocking, or ParentScoring
    Returns:
    :returns: dict child_dict: Dictionary of all the class objects for either
        Filtering, docking, Dockingfile conversion or scoring
    """
    parent_object = None
    get_all_subclasses = None
    if purpose_of_object == "parent_pdbqt_converter":
        import autogrow.docking.docking_class.docking_file_conversion
        from autogrow.docking.docking_class.parent_pdbqt_converter import (
            ParentPDBQTConverter as parent_object,
        )
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses
    elif purpose_of_object == "ParentScoring":
        import autogrow.docking.scoring.scoring_classes.scoring_functions
        from autogrow.docking.scoring.scoring_classes.parent_scoring_class import (
            ParentScoring as parent_object,
        )
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses

    assert parent_object is not None, "parent_object is not recognized"
    assert get_all_subclasses is not None, "get_all_subclasses is not recognized"

    children = get_all_subclasses(parent_object)
    child_dict = {}
    for child in children:
        child_object = child()
        child_name = child_object.get_name()
        child_dict[child_name] = child_object

    return child_dict


############################################
######## Filter Handlining Settings ########
############################################
def filter_choice_handling(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    This function handles selecting the user defined Ligand filters.

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: dict params: Dictionary of User variables with the
        chosen_ligand_filters added
    """
    if "No_Filters" in list(params.keys()):
        if params["No_Filters"] is True:
            chosen_ligand_filters = None
        else:
            chosen_ligand_filters, params = picked_filters(params)
    else:
        chosen_ligand_filters, params = picked_filters(params)
    params["chosen_ligand_filters"] = chosen_ligand_filters

    return params


#
def picked_filters(params: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
    """
    This will take the user vars and return a list of the filters
    which a molecule must pass to move into the next generation.

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow
        the same nomenclature and enter that filter in the user
        params["alternative_filters"] as the name of that class and place
        that file in the same folder as the other filter classes.
    """
    import pdb

    pdb.set_trace()
    filter_list = []
    vars_keys = list(params.keys())

    if "LipinskiStrictFilter" in vars_keys:
        if params["LipinskiStrictFilter"] is True:
            filter_list.append("LipinskiStrictFilter")
    else:
        params["LipinskiStrictFilter"] = False

    if "LipinskiLenientFilter" in vars_keys:
        if params["LipinskiLenientFilter"] is True:
            filter_list.append("LipinskiLenientFilter")
    else:
        params["LipinskiLenientFilter"] = False

    if "GhoseFilter" in vars_keys:
        if params["GhoseFilter"] is True:
            filter_list.append("GhoseFilter")
    else:
        params["GhoseFilter"] = False

    if "GhoseModifiedFilter" in vars_keys:
        if params["GhoseModifiedFilter"] is True:
            filter_list.append("GhoseModifiedFilter")
    else:
        params["GhoseModifiedFilter"] = False

    if "MozziconacciFilter" in vars_keys:
        if params["MozziconacciFilter"] is True:
            filter_list.append("MozziconacciFilter")
    else:
        params["MozziconacciFilter"] = False

    if "VandeWaterbeemdFilter" in vars_keys:
        if params["VandeWaterbeemdFilter"] is True:
            filter_list.append("VandeWaterbeemdFilter")
    else:
        params["VandeWaterbeemdFilter"] = False

    if "PAINSFilter" in vars_keys:
        if params["PAINSFilter"] is True:
            filter_list.append("PAINSFilter")
    else:
        params["PAINSFilter"] = False

    if "NIHFilter" in vars_keys:
        if params["NIHFilter"] is True:
            filter_list.append("NIHFilter")
    else:
        params["NIHFilter"] = False

    if "BRENKFilter" in vars_keys:
        if params["BRENKFilter"] is True:
            filter_list.append("BRENKFilter")
    else:
        params["BRENKFilter"] = False

    if "alternative_filter" in vars_keys:
        filter_list = handle_alternative_filters(params, filter_list)
    else:
        params["alternative_filter"] = None

    # if there is no user specified ligand filters but they haven't set
    # filters to None ---> set filter to default of LipinskiLenientFilter.
    # TODO: I want the default to be none. Revisit this.
    if len(filter_list) == 0:
        params["LipinskiLenientFilter"] = True
        filter_list.append("LipinskiLenientFilter")

    return filter_list, params
