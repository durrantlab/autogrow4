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
from typing import Any, Dict, List, Tuple

from autogrow.config.config_run_directory import set_run_directory
from autogrow.utils.shellcmds import determine_bash_timeout_vs_gtimeout

# TODO: NOTE USED. CAN DELETE?
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
    from autogrow.utils.parallelizer import Parallelizer

    params["parallelizer"] = Parallelizer(
        params["multithread_mode"], params["number_of_processors"]
    )

    return params


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
        "receptor_path",
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
    input_params["receptor_path"] = os.path.abspath(input_params["receptor_path"])
    input_params["root_output_folder"] = os.path.abspath(
        input_params["root_output_folder"]
    )
    input_params["source_compound_file"] = os.path.abspath(
        input_params["source_compound_file"]
    )

    # Check receptor_path exists
    if os.path.isfile(input_params["receptor_path"]) is False:
        raise NotImplementedError(
            "Receptor file can not be found. File must be a .PDB file."
        )
    if ".pdb" not in input_params["receptor_path"]:
        raise NotImplementedError("receptor_path must be a .PDB file.")

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


# TODO: Rename this here and in `check_for_required_inputs`
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



# TODO: NEVER USED. DELETE?
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


############################################
######## Input Handlining Settings #########
############################################
# TODO: Never used?
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
            # Examples may be things like receptor_path or
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

