"""
Path Configuration Module for AutoGrow

This module handles the configuration of file and directory paths for the
AutoGrow program, ensuring proper path formatting and directory creation.
"""
import os
from typing import Any, Dict


def config_paths(params: Dict[str, Any]) -> None:
    """
    Configure paths for the AutoGrow program.

    This function performs the following actions:
    1. Converts specified path parameters to absolute paths.
    2. Ensures directory paths end with the OS-specific path separator.
    3. Creates the output folder if it doesn't exist.

    Args:
        params (Dict[str, Any]): Dictionary of user variables governing how the
            program runs.

    The function modifies the `params` dictionary in place.
    """
    _make_paths_abs(params)
    _make_dirs_end_in_sep(params)
    _create_output_folder(params["root_output_folder"])


def _make_paths_abs(params: Dict[str, Any]) -> None:
    """
    Convert specified path parameters to absolute paths.

    Args:
        params (Dict[str, Any]): Dictionary of user variables.

    The function modifies the `params` dictionary in place, converting paths
    for 'receptor_path', 'root_output_folder', and 'source_compound_file' to
    absolute paths.
    """
    # convert paths to abspath, in case necessary
    for pname in ["receptor_path", "root_output_folder", "source_compound_file"]:
        if pname not in list(params.keys()):
            continue
        params[pname] = os.path.abspath(params[pname])


def _make_dirs_end_in_sep(params: Dict[str, Any]) -> None:
    """
    Ensure specified directory paths end with the OS-specific path separator.

    Args:
        params (Dict[str, Any]): Dictionary of user variables.

    The function modifies the `params` dictionary in place, appending the
    OS-specific path separator to 'root_output_folder' if necessary.
    """
    dir_params = ["root_output_folder"]

    for dir_param in dir_params:
        if dir_param not in params:
            continue
        if params[dir_param][-1] != os.sep:
            params[dir_param] = params[dir_param] + os.sep


def _create_output_folder(out_folder: str) -> None:
    """
    Create the output folder if it doesn't exist.

    Args:
        out_folder (str): The path to the output folder.

    Raises:
        NotImplementedError: If the folder cannot be created.
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
