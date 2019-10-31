"""
This script will convert a pdbqt file into a .pdb file.

This is done by removing a column of the PDB file.
"""

import os
import sys

def convert_pdbqt_to_pdb(pdbqt_file_in, pdb_file_out):
    """
    Converts a pdbqt file to a pdb file by removing the 3rd to last column.
    Input:
    :param str pdbqt_file_in: the string of .pdbqt to be formated
    :param str pdb_file_out: the string of the output .pdb
    """
    os.system("cut -c 1-60,70-79 {} > {}".format(pdbqt_file_in, pdb_file_out))



#######
if __name__ == "__main__":
    try:
        pdbqt_file = sys.argv[1]
    except:
        print("Need a path to a pdbqt to be converted to pdb")
    if os.path.exists(pdbqt_file) == False:
        PRINTOUT = "pdbqt file does not exists: {}".format(pdbqt_file)
        print(PRINTOUT)
        raise Exception(PRINTOUT)

    if ".pdbqt" not in pdbqt_file and ".PDBQT" not in pdbqt_file:
        PRINTOUT = "pdbqt file must be a .pdbqt file.: {}".format(pdbqt_file)
        print(PRINTOUT)
        raise Exception(PRINTOUT)
    if pdbqt_file[-2].lower() != "q" or pdbqt_file[-1].lower() != "t":
        PRINTOUT = "pdbqt file must be a .pdbqt file.: {}".format(pdbqt_file)
        print(PRINTOUT)
        raise Exception(PRINTOUT)

    try:
        pdb_file = sys.argv[2]
        print("converting {} to .pdb".format(pdbqt_file))
        print("\tThe new file will be {}".format(pdb_file))
    except:
        print("PDB outputfile will be the same as pdbqt file minus qt")
        pdb_file = pdbqt_file.replace(".pdbqt", ".pdb")
        pdb_file = pdbqt_file.replace(".PDBQT", ".PDB")
        print("\t pdbqt_file = {} \t pdb_file = {}".format(pdbqt_file, pdb_file))
    

    convert_pdbqt_to_pdb(pdbqt_file, pdb_file)
    print("finished")