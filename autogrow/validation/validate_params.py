import os
from typing import Any, Dict, List


def validate_params(params: Dict[str, Any]) -> None:
    """
    Check for missing variables in the required inputs.

    Inputs:
    :param dict params: The parameters. A dictionary of {parameter name: value}.
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
        "root_output_folder",
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

    if os.path.isdir(params["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )


def _throw_missing_params_error(missing_params: List[str]) -> None:
    """
    Raises an error if any required parameters are missing.

    This function generates and prints a message listing all the missing parameters
    that are required for the program to run. It raises a NotImplementedError with 
    the list of missing parameters and instructions for how to obtain help.

    Inputs:
    :param List[str] missing_params: A list of the required parameter names that 
                                     are missing from the input.
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
    Verifies that a specified file exists and has the correct extension.

    This function checks whether the file specified by the parameter name exists 
    in the provided dictionary of parameters. It also ensures the file has the 
    expected extension (e.g., ".pdb" for PDB files). If the file does not exist or 
    has an incorrect extension, a NotImplementedError is raised.

    Inputs:
    :param str pname: The key in the params dictionary representing the file path.
    :param str ext: The expected file extension (e.g., "pdb", "smi").
    :param str file_desc: A description of the expected file type (e.g., ".PDB file").
    :param Dict[str, Any] params: A dictionary of parameter names and values.

    Raises:
    :raises NotImplementedError: If the file does not exist or does not have the expected extension.
    """
    if os.path.isfile(params[pname]) is False:
        raise NotImplementedError(
            f"{pname} can not be found. \
            File must be a {file_desc}."
        )
    if f".{ext}" not in params[pname]:
        raise NotImplementedError(f"{pname} must be a {file_desc}.")
