import os
from typing import Any, Dict


def config_paths(params: Dict[str, Any]) -> None:
    """
    Configures paths by converting them to absolute paths, ensuring directories end with a separator,
    and creating the output folder if it does not exist.

    This function performs the following actions:
    1. Converts specified path parameters in `params` to absolute paths.
    2. Ensures that directory paths end with the OS-specific path separator.
    3. Creates the output folder if it does not already exist.

    Inputs:
    :param params: Dictionary of user variables governing how the program runs.

    Returns:
    :returns: None. The function modifies the `params` dictionary in place.
    """
    _make_paths_abs(params)
    _make_dirs_end_in_sep(params)
    _create_output_folder(params["root_output_folder"])


def _make_paths_abs(params: Dict[str, Any]) -> None:
    """
    Converts specified path parameters in `params` to absolute paths.

    This helper function iterates over a predefined list of parameter names and converts their
    values to absolute paths using `os.path.abspath`. It ensures that all relevant paths are
    absolute, which helps avoid issues related to relative paths during program execution.

    Inputs:
    :param params: Dictionary of user variables.

    Returns:
    :returns: None. The function modifies the `params` dictionary in place.
    """
    # convert paths to abspath, in case necessary
    for pname in [
        "filename_of_receptor",
        "root_output_folder",
        "source_compound_file"
    ]:
        if pname not in list(params.keys()):
            continue
        params[pname] = os.path.abspath(params[pname])


def _make_dirs_end_in_sep(params: Dict[str, Any]) -> None:
    """
    Ensures that specified directory paths in `params` end with the OS-specific path separator.

    This helper function iterates over a predefined list of directory parameter names and appends
    the OS-specific path separator (`os.sep`) to the end of each directory path if it does not
    already end with one. This standardizes directory paths, making file operations more predictable.

    Inputs:
    :param params: Dictionary of user variables.

    Returns:
    :returns: None. The function modifies the `params` dictionary in place.
    """
    dir_params = [
        "root_output_folder"
    ]

    for dir_param in dir_params:
        if dir_param not in params:
            continue
        if params[dir_param][-1] != os.sep:
            params[dir_param] = params[dir_param] + os.sep


def _create_output_folder(out_folder: str) -> None:
    """
    Creates the output folder if it does not already exist.

    This helper function checks whether the specified `out_folder` exists. If it does not, the
    function attempts to create the directory using `os.makedirs`. If directory creation fails,
    it raises a `NotImplementedError` with an informative message.

    Inputs:
    :param out_folder: The path to the output folder.

    Returns:
    :returns: None. The function may raise a `NotImplementedError` if the folder cannot be created.
    """
    # Check root_output_folder exists
    if os.path.exists(out_folder) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(out_folder)
        except Exception as e:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            ) from e
