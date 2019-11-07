"""
This script will convert a pdbqt file into a .pdb file.

This is done by removing a column of the PDB file.

# Run example:
# 
# output example:
# python convert_ligands_pdb_to_smi.py \
#   -source_folder $PATH/OF/PDBS/ \
#   -output_folder $PATH/TO/OUTPUT/ \
#   -number_of_processors -1


"""

import __future__
import os
import sys
import argparse


def convert_pdbqt_to_pdb(pdbqt_file_in, pdb_file_out):
    """
    Converts a pdbqt file to a pdb file by removing the 3rd to last column.
    Input:
    :param str pdbqt_file_in: the string of .pdbqt to be formated
    :param str pdb_file_out: the string of the output .pdb
    """
    printout = ""
    line_index_range = [x for x in range(0,61)] + [x for x in range(70,80)] 
    print(line_index_range)
    with open(pdbqt_file_in) as f:
        for line in f.readlines():
            if "ATOM" in line:
                print(line)
                short_line = ""
                for i in line_index_range:
                    # print(i)
                    if i >= len(line):continue
                        
                    else:
                        short_line = short_line + line[i]

                print(short_line)
                printout = printout + short_line 
            else:
                printout = printout + line
    with open(pdb_file_out,'w') as f:
        f.write(printout)


def get_arguments_from_argparse(ARGS_DICT):
    """
    This function handles the arg parser arguments for the script.
    
    Inputs:
    :param dict ARGS_DICT: dictionary of parameters
    Returns:
    :returns: dict ARGS_DICT: dictionary of parameters
    """


    # Argument handling
    if type(ARGS_DICT["pdbqt_file"]) != str:
        raise Exception("provided pdbqt_file must be a .pdbqt file.")

    #  argument_handling
    if os.path.exists(ARGS_DICT["pdbqt_file"]) == False:
        raise Exception("provided pdbqt_file must be a .pdbqt file.")
    else:
        if ARGS_DICT["pdbqt_file"].split(".")[-1] != "pdbqt" and ARGS_DICT["pdbqt_file"].split(".")[-1] != "PDBQT":

            raise Exception("provided pdbqt_file must be a .pdbqt file.")
        
    if ARGS_DICT["output_file"] != None:
        if ARGS_DICT["output_file"].split(".")[-1] != "pdb" and ARGS_DICT["output_file"].split(".")[-1] != "PDB":
            raise Exception("provided output_file must be a .pdb file.") 

        if os.path.exists(os.path.dirname(ARGS_DICT["output_file"])) == False:
            try:
                os.mkdir(os.path.dirname(ARGS_DICT["output_file"]))
            except:
                pass
            if os.path.exists(os.path.dirname(ARGS_DICT["output_file"])) == False:
                raise Exception("directory to output the file could not be made or found.")
    else:
        ARGS_DICT["output_file"] = os.path.dirname(ARGS_DICT["pdbqt_file"]) + \
            ARGS_DICT["pdbqt_file"].replace(".pdbqt", ".pdb").replace(".PDBQT", ".pdb")

    return ARGS_DICT



# Argment parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument('--pdbqt_file', '-f', required = False, default=None,
    help='Path to .pdbqt file to convert to a .pdb file. This must be a single \
    ligand and must end with .pdbqt')
PARSER.add_argument('--output_file', '-o', type=str, default=None,
    help='Path to file where we will output .pdb file. \
    If not provide the output .pdb will be the same as the input \
    pdbqt_file but ending with .pdb instead of .pdbqt.')

ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)


# Run Converter
convert_pdbqt_to_pdb(ARGS_DICT["pdbqt_file"], ARGS_DICT["output_file"])
print("finished")