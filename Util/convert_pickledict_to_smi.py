#This script will filter out anything which fails to be imported into rdkit
# sanitized, protanated/deprotanated... and will make a new .smi with the same 
# SMILES and ZincID. It will also conicalize the SMILES strings.
# It will correct for Nitrogens with improper valence charges.

# This also Removes STEREO INFORMATION!!!!!
# import __future__

import rdkit
from rdkit import Chem
import copy
import json
import sys
import os
import pickle
import random

# import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp
# import run_filter as rf

import glob


import rdkit
import rdkit.Chem as Chem


def get_mols_dict_from_pickle(file_path):

    with open(file_path, 'rb') as handle:
        mols_dict = pickle.load(handle)
    return mols_dict

def write_pickle_to_file(file_path, obj):
    with open(file_path, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
        

def convert_to_smi(pickle_infile, smi_outfile):
    mols_dict = get_mols_dict_from_pickle(pickle_infile)
    
    if type(mols_dict) != dict:
        raise Exception("FILE IS NOT A DICT: {}".format(mols_dict))

    list_keys = list(mols_dict.keys())
    list_keys.sort()


    with open(smi_outfile, 'w') as f:
        for ZincID in list_keys:
            printout = ""
            mol_list = mols_dict[ZincID]
            printout = mol_list[0] + "\t" + mol_list[1] + "\n"
            f.write(printout)
    

if __name__ == "__main__":

    infolder = sys.argv[1]
    outfolder = sys.argv[2]

    jobs = []

    if os.path.exists(infolder) == False:
        raise Exception("infolder doesn't exist: {}".format(infolder))
    if os.path.exists(outfolder) == False:
        raise Exception("outfolder doesn't exist: {}".format(outfolder))

    if os.path.isdir(infolder) == False:
        raise Exception("infolder isn't a directory: {}".format(infolder))
    if os.path.isdir(outfolder) == False:
        raise Exception("outfolder isn't a directory: {}".format(outfolder))

    for infile in glob.glob(infolder + "*_pickle"):
        base = os.path.basename(infile)
        smi_outfile = outfolder + base.split("_pickle")[0] + ".smi"
        if os.path.exists(smi_outfile)== True:
            print("")
            print("File already exists for {}".format(smi_outfile))
            print("skipping {}".format(base))
            continue

        elif os.path.exists(smi_outfile)== False:
            
            jobs.append((infile,smi_outfile))

    print("Number of jobs", len(jobs))
    # convert_to_smi(infile, smi_outfile)
    output = mp.MultiThreading(jobs, -1, convert_to_smi)    

