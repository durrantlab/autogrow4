"""Parameter validation utilities for AutoGrow.

This module handles validation of required parameters and input files, ensuring
that all necessary configuration is present and properly formatted before AutoGrow
execution begins.
"""
import os
from typing import Any, Dict, List


def validate_params(params: Dict[str, Any]) -> None:
    """
    Validate required input parameters and files.

    Checks for presence of all required parameters and verifies that specified
    input files exist with correct formats. Required parameters include receptor
    configuration, grid specifications, output folder, and source compounds.

    Args:
        params (Dict[str, Any]): Dictionary of parameter names and values

    Raises:
        NotImplementedError: If required parameters are missing, input files
            don't exist, or files have incorrect formats
    """
    keys_from_input = list(params.keys())

    list_of_required_params = [
        "receptor_path",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "output_directory",
        "source_compound_file",
    ]

    missing_params = [
        param for param in list_of_required_params if param not in keys_from_input
    ]

    if missing_params:
        _throw_missing_params_error(missing_params)

    # Check receptor_path and source_compound_file exist
    _check_file_exists("receptor_path", "pdb", ".PDB file", params)
    _check_file_exists(
        "source_compound_file", "smi", "tab delineated .smi file", params
    )

    if os.path.isdir(params["output_directory"]) is False:
        raise NotImplementedError(
            "output_directory is not a directory. \
            Check your input parameters."
        )


def _throw_missing_params_error(missing_params: List[str]) -> None:
    """
    Raise an error with details about missing required parameters.

    Generates a detailed error message listing all missing required parameters
    and instructions for getting help with parameter usage.

    Args:
        missing_params (List[str]): List of required parameter names missing
            from input

    Raises:
        NotImplementedError: Always raised with message detailing missing
            parameters and help instructions
    """
    printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunAutogrow -h"
    printout += "\nThe following required variables are missing: "
    for param in missing_params:
        printout += "\n\t" + param
    print("")
    print(printout)
    print("")
    raise NotImplementedError("\n" + printout + "\n")


def _check_file_exists(
    pname: str, ext: str, file_desc: str, params: Dict[str, Any]
) -> None:
    """
    Verify existence and format of required input files.

    Args:
        pname (str): Parameter name for file path in params dictionary
        ext (str): Expected file extension without dot (e.g. "pdb", "smi")
        file_desc (str): Description of expected file type 
        params (Dict[str, Any]): Dictionary of parameter names and values

    Raises:
        NotImplementedError: If file doesn't exist or has incorrect extension

    Example:
        >>> _check_file_exists("receptor_path", "pdb", ".PDB file", params)
        >>> _check_file_exists("source_compound_file", "smi", 
        ...                   "tab delineated .smi file", params)
    """
    if os.path.isfile(params[pname]) is False:
        raise NotImplementedError(
            f"{pname} can not be found. \
            File must be a {file_desc}."
        )
    if f".{ext}" not in params[pname]:
        raise NotImplementedError(f"{pname} must be a {file_desc}.")
