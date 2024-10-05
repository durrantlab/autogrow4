import os
from typing import Any, Dict


def config_params_under_windows(params: Dict[str, Any], printout: str) -> str:
    """
    Configures parameters specifically for Windows environments.

    If the operating system is Windows (`os.name` is either "nt" or "ce"), this function disables multiprocessing
    by setting `number_of_processors` to 1 and appends a warning message to the `printout` string. Additionally,
    it ensures that any path parameters containing spaces are enclosed in quotes to handle paths with spaces
    correctly.

    Inputs:
    :param params: Dictionary of user variables governing how the program runs.
    :param printout: A string message that accumulates warnings and information to be printed.

    Returns:
    :returns: The updated `printout` string with any added warnings or messages.
    """
    if os.name not in {"nt", "ce"}:
        # It's not windows
        return ""

    # so it's running under windows. multiprocessing disabled
    params["number_of_processors"] = 1
    printout += "\nWARNING: Multiprocessing is disabled on windows machines.\n"

    # convert path names with spaces if this is windows
    params_that_are_paths = [
        "filename_of_receptor",
        "root_output_folder",
    ]

    for pname in params_that_are_paths:
        if pname in params and " " in params[pname]:
            params[pname] = f'"{params[pname]}"'

    return printout
