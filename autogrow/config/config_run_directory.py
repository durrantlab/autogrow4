import os
from typing import Optional

from autogrow.utils.logging import LogLevel, log_info


def set_run_directory(root_folder_path: str) -> str:
    """
    Determine and make the folder for the run directory.

    Inputs:
    :param str root_folder_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files

    Returns:
    :returns: str folder_path: the string of the newly created directory for
        puting output folders
    """

    if not os.path.exists(root_folder_path):
        os.makedirs(root_folder_path)
        log_info(f"Making the output folder path: {root_folder_path}")

    return root_folder_path
