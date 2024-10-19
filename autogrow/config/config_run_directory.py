"""
Run Directory Configuration Module for AutoGrow

This module handles the creation and configuration of the run directory for
AutoGrow output files.
"""
import os
from typing import Optional

from autogrow.utils.logging import LogLevel, log_info


def set_run_directory(root_folder_path: str) -> str:
    """
    Determine and create the folder for the run directory.

    This function checks if the specified root folder exists. If it doesn't,
    it creates the directory and logs the action. It then returns the path
    to this directory.

    Args:
        root_folder_path (str): The path of the root output folder where the
            run directory will be created.

    Returns:
        str: The path to the created or existing run directory.

    Note:
        If the directory already exists, this function will not create a new
        one but will return the path to the existing directory.
    """
    # TODO: This can be removed, and code merged into original function
    if not os.path.exists(root_folder_path):
        os.makedirs(root_folder_path)
        log_info(f"Making the output folder path: {root_folder_path}")

    return root_folder_path
