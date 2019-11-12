# convert pdbs into smiles

# This script will take a folder and convert all pdb files into a single texted file.
# The text file will contain smiles strings of the respective pdb and the name of the file.

# Run example:
# 
# output example:
# python convert_ligands_pdb_to_smi.py \
#   -source_folder $PATH/OF/PDBS/ \
#   -output_folder $PATH/TO/OUTPUT/ \
#   -number_of_processors -1
# 
# This will convert all .pdb files within $PATH/OF/PDBS/ into a single 
# .smi file at $PATH/TO/OUTPUT/PDBS.smi
# using all available processors
# 
# CC1COC(=O)OC(=O)O1    ZINC60039447
# O = C1OC(=O)N2CCC12    ZINC59199492
# O = C1CC(C(=O)O)C(=O)O1    ZINC59901386
import __future__

import glob
import os
import sys

import argparse

from rdkit import Chem
       
import support_scripts.mol_object_handling as MOH
import support_scripts.Multiprocess as mp

def run_convert_on_single_pdb(pdb):

    """
    This function converts a ligand into SMILES
    and returns the list with the smiles with a name. 
    The names are the basename of the file minus the .pdb

    Inputs:
    :param str pdb: path to the folder to a pdb file
    Returns:
    :returns: list output_data: A list containing all SMILES info from the file
    """
    
    try:
        mol = Chem.MolFromPDBFile(pdb)

        mol_sanitized = MOH.check_sanitization(mol)
        if mol_sanitized is not None:
            smiles = Chem.MolToSmiles(mol_sanitized, isomericSmiles=True)
            fileName = os.path.basename(pdb)
            fileStripped = fileName.replace(".pdb","")
            output_data = smiles + "\t" + fileStripped
    except:
        pass
    return output_data
# 

def make_smile_list(sub_folder, number_of_processors):
    """
    This function converts every ligand within a folder into SMILES
    and returns the list of smiles with a name. 
    The names are the basename of the file minus the .pdb

    Inputs:
    :param str sub_folder: path to the folder to search for pdb files
    Returns:
    :returns: list smilesList: A list of lists containing all SMILES from
        the .pdb files and their respective name
    """
    sub_folder = sub_folder+os.sep
    smilesList = []
    pdb_list = glob.glob(os.sep + sub_folder+"*.pdb")
    pdb_list.extend(glob.glob(os.sep + sub_folder+"*.PDB"))
    pdb_list = tuple([tuple([pdb]) for pdb in pdb_list])

    # run convert in multithread
    smilesList = mp.multi_threading(pdb_list, -1,  run_convert_on_single_pdb)

    return smilesList
# 

def get_arguments_from_argparse(ARGS_DICT):
    """
    This function handles the arg parser arguments for the script.
    
    Inputs:
    :param dict ARGS_DICT: dictionary of parameters
    Returns:
    :returns: dict ARGS_DICT: dictionary of parameters
    """


    # Argument handling
    if type(ARGS_DICT["source_folder"]) != str:
        raise Exception("provided source folder must be a directory.")
    if type(ARGS_DICT["output_folder"]) != str:
        raise Exception("provided output_folder must be a directory.")

    #  argument_handling
    if os.path.exists(ARGS_DICT["source_folder"]) == False or os.path.isdir(ARGS_DICT["source_folder"])==False:
        raise Exception("provided source folder can not be found or is not a directory.")
    else:
        ARGS_DICT["source_folder"] = os.path.abspath(ARGS_DICT["source_folder"]) + os.sep
        
    if os.path.exists(ARGS_DICT["output_folder"]) == False:
        try:
            os.mkdir(ARGS_DICT["output_folder"])
        except:
            pass
        if os.path.exists(ARGS_DICT["output_folder"]) == False:
            raise Exception("output_folder could not be made or found.")
    else:
        if os.path.isdir(ARGS_DICT["output_folder"]) == False:
            raise Exception("output_folder needs to be a directory.")
        else:
            ARGS_DICT["output_folder"] = os.path.abspath(ARGS_DICT["output_folder"]) + os.sep

    return ARGS_DICT
# 

# Argment parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument('--source_folder', '-s', required = False, default=None,
    help='Path to folder containing .pdb files to convert. \
    File must contain a single small molecules. Without protein. \
    Files must end with either .pdb or .PDB')
PARSER.add_argument('--output_folder', '-o', required = False, default=None,
    help='Path to folder where we will output .smi file of converted .pdbs.')

# processors and multithread mode
PARSER.add_argument('--number_of_processors', '-p', type = int, metavar='N', default = 1,
    help='Number of processors to use for parallel calculations. \
    This script is not MPI enable but is able to multithread using SMP architecture. \
    Set to -1 for all availble CPUs.')

ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)

# Running converter
smilesList = make_smile_list(ARGS_DICT["source_folder"], ARGS_DICT["number_of_processors"])
name = [x for x in ARGS_DICT["source_folder"].split(os.sep)if x!=""][-1] 
outputFile = ARGS_DICT["output_folder"] + os.sep + name + ".smi"
with open(outputFile,"w") as f:
    f.write("\n".join(smilesList))

print("Converted ligands to .smi located:\n\t{}".format(outputFile))
        