import os

def validate_params(params: dict):
    """
    Check for missing variables in the required inputs.

    Inputs:
    :param dict params: The parameters. A dictionary of {parameter name: value}.
    """
    keys_from_input = list(params.keys())

    list_of_required_params = [
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

    missing_params = [
        param
        for param in list_of_required_params
        if param not in keys_from_input
    ]

    if missing_params:
        _throw_missing_params_error(missing_params)

    # Check filename_of_receptor and source_compound_file exist
    _check_file_exists("filename_of_receptor", "pdb", ".PDB file", params)
    _check_file_exists(
        "source_compound_file", "smi", "tab delineated .smi file", params
    )

    if os.path.isdir(params["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )


def _throw_missing_params_error(missing_params):
    printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunAutogrow -h"
    printout += "\nThe following required variables are missing: "
    for param in missing_params:
        printout += "\n\t" + param
    print("")
    print(printout)
    print("")
    raise NotImplementedError("\n" + printout + "\n")

def _check_file_exists(pname: str, ext: str, file_desc: str, params: dict):
    if os.path.isfile(params[pname]) is False:
        raise NotImplementedError(
            f"{pname} can not be found. \
            File must be a {file_desc}."
        )
    if f".{ext}" not in params[pname]:
        raise NotImplementedError(f"{pname} must be a {file_desc}.")

