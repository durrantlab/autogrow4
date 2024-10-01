"""
This script will convert a pdbqt file into a .pdb file.

This is done by removing a column of the PDB file.

# Run example:

# output example:
# python convert_ligands_pdb_to_smi.py \
#   -source_folder $PATH/OF/PDBS/ \
#   -output_folder $PATH/TO/OUTPUT/ \
#   -number_of_processors -1


"""


import __future__
import contextlib
import os
import argparse
from typing import Dict


def convert_pdbqt_to_pdb(pdbqt_file_in: str, pdb_file_out: str) -> None:
    """
    Converts a pdbqt file to a pdb file by removing the 3rd to last column.
    Inputs:
    :param str pdbqt_file_in: the string of .pdbqt to be formatted
    :param str pdb_file_out: the string of the output .pdb
    """
    printout = ""
    line_index_range = list(range(61)) + list(range(70, 80))

    with open(pdbqt_file_in) as f:
        for line in f:
            if "ATOM" in line:
                short_line = ""
                for i in line_index_range:
                    # print(i)
                    if i >= len(line):
                        continue

                    short_line = short_line + line[i]

                printout = printout + short_line
            elif (
                "REMARK                            x       y       z     vdW  Elec"
                + "       q    Type"
                in line
                or "REMARK                         _______ _______ _______ _____ _____"
                + "    ______ ____"
                in line
            ):
                short_line = ""
                for i in line_index_range:
                    # print(i)
                    if i >= len(line):
                        continue

                    short_line = short_line + line[i]

                printout = printout + short_line + "\n"
            else:
                printout = printout + line
    with open(pdb_file_out, "w") as f:
        f.write(printout)


def get_arguments_from_argparse(args_dict: Dict[str, str]) -> Dict[str, str]:
    """
    This function handles the arg parser arguments for the script.

    Inputs:
    :param dict args_dict: dictionary of parameters
    Returns:
    :returns: dict args_dict: dictionary of parameters
    """

    # Argument handling
    if type(args_dict["pdbqt_file"]) != str:
        raise Exception("provided pdbqt_file must be a .pdbqt file.")

    if not os.path.exists(args_dict["pdbqt_file"]):
        raise Exception("provided pdbqt_file must be a .pdbqt file.")
    if args_dict["pdbqt_file"].split(".")[-1] not in ["pdbqt", "PDBQT"]:

        raise Exception("provided pdbqt_file must be a .pdbqt file.")

    if args_dict["output_file"] is not None:
        if args_dict["output_file"].split(".")[-1] not in ["pdb", "PDB"]:
            raise Exception("provided output_file must be a .pdb file.")

        if not os.path.exists(os.path.dirname(args_dict["output_file"])):
            with contextlib.suppress(Exception):
                os.mkdir(os.path.dirname(args_dict["output_file"]))
            if not os.path.exists(os.path.dirname(args_dict["output_file"])):
                raise Exception(
                    "directory to output the file could not be made or found."
                )
    else:
        args_dict["output_file"] = os.path.dirname(args_dict["pdbqt_file"]) + args_dict[
            "pdbqt_file"
        ].replace(".pdbqt", ".pdb").replace(".PDBQT", ".pdb")

    return args_dict


# Argument parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument(
    "--pdbqt_file",
    "-f",
    required=True,
    default=None,
    help="Path to .pdbqt file to convert to a .pdb file. This must be a single \
    ligand and must end with .pdbqt",
)
PARSER.add_argument(
    "--output_file",
    "-o",
    type=str,
    default=None,
    help="Path to file where we will output .pdb file. \
    If not provide the output .pdb will be the same as the input \
    pdbqt_file but ending with .pdb instead of .pdbqt.",
)

ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)


# Run Converter
convert_pdbqt_to_pdb(ARGS_DICT["pdbqt_file"], ARGS_DICT["output_file"])
print("finished")
