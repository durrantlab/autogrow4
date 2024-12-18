"""Utilities for file format conversion using OpenBabel.

This module provides functions to construct and execute OpenBabel commands for
converting between different molecular file formats.
"""
import os

from autogrow.utils.logging import log_warning


def obabel_convert_cmd(
    in_file: str, out_file: str, obabel_path: str, extra_params: str = ""
) -> str:
    """
    Construct an OpenBabel conversion command string.

    Creates a shell command for OpenBabel to convert between molecular file
    formats. Handles special cases like 'vina' extension and supports additional
    parameters.

    Args:
        in_file (str): Path to input file
        out_file (str): Path for output file
        obabel_path (str): Path to OpenBabel executable
        extra_params (str, optional): Additional OpenBabel parameters. Defaults
            to empty string

    Returns:
        str: Complete OpenBabel conversion command with error output suppressed

    Note:
        - Treats '.vina' extension as '.pdbqt'
        - Redirects stdout and stderr to /dev/null
        - Uses '-e' flag for error checking
    """
    # prot_and_3d: bool = False
    in_ext = in_file.split(".")[-1]
    out_ext = out_file.split(".")[-1]

    if in_ext == "vina":
        in_ext = "pdbqt"

    cmd = f'{obabel_path} -i{in_ext} "{in_file}" -o{out_ext}'
    if extra_params != "":
        # cmd += " --gen3d --p 7.4"
        cmd += f" {extra_params}"
    cmd += f' -e -O "{out_file}"'
    cmd += "  > /dev/null 2>&1"

    return cmd


def obabel_convert(
    in_file: str, out_file: str, obabel_path: str, extra_params: str = ""
) -> bool:
    """Execute an OpenBabel conversion between file formats.

    Attempts to convert a molecular file from one format to another using
    OpenBabel. Validates the conversion by checking the output file exists and
    is non-empty.

    Args:
        in_file (str): Path to input file
        out_file (str): Path for output file
        obabel_path (str): Path to OpenBabel executable
        extra_params (str, optional): Additional OpenBabel parameters. Defaults
            to empty string

    Returns:
        bool: True if conversion succeeds, False otherwise

    Note:
        Logs warnings via log_warning() if conversion fails for any reason:
        - Command execution fails
        - Output file is not created
        - Output file is empty
    """
    cmd = obabel_convert_cmd(in_file, out_file, obabel_path, extra_params)
    try:
        os.system(cmd)
    except Exception as e:
        log_warning(f"Could not convert with obabel: {in_file}")
        return False

    if not os.path.exists(out_file):
        log_warning(f"Could not convert with obabel: {in_file}")
        return False

    with open(out_file, "r") as f:
        content = f.read().strip()
        if content == "":
            log_warning(f"Could not convert with obabel: {in_file}")
            return False

    return True
