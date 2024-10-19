"""Utility functions for handling JSON parameters and saving variables."""
import json
import copy
import os
import sys
from typing import Dict, Union


def convert_json_params_from_unicode(
    params_unicode: Dict[str, Union[str, int]]
) -> Dict[str, Union[str, int]]:
    """
    Convert Unicode parameters to ASCII.

    Args:
        params_unicode (Dict[str, Union[str, int]]): The parameters dictionary
            with potentially Unicode keys and values.

    Returns:
        Dict[str, Union[str, int]]: A new dictionary with ASCII-encoded keys
        and values.

    Note:
        This function is particularly useful when working with JSON-loaded
        data, which may contain Unicode strings.
    """
    # Also, rdkit doesn't play nice with unicode, so convert to ascii

    # TODO: Don't param and val need to be encoded? This doesn't seem to
    # actually do anything...

    params = {}
    for param, val in params_unicode.items():
        params[param] = val
    return params


def save_vars_as_json(params: Dict[str, Union[str, int]]) -> None:
    """
    Save the parameters dictionary as a JSON file.

    This function saves the input parameters dictionary as a JSON file, which
    can be used to track experiments and is necessary for several utility
    scripts. It excludes the 'parallelizer' object from the saved data.

    Args:
        params (Dict[str, Union[str, int]]): Dictionary of user variables that
            govern how the program runs.

    Note:
        - The file is saved as 'vars.json' in the output directory specified in
          the params dictionary.
        - If 'vars.json' already exists, it appends a number to the filename
          (e.g., 'vars_2.json', 'vars_3.json', etc.).
        - Utility scripts will only look at the original 'vars.json' file.
    """
    output_directory = str(params["output_directory"])

    vars_file = output_directory + os.sep + "vars.json"
    if os.path.exists(vars_file):
        # vars.json already exists. lets make the next file.
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
