import glob
import os
import sys
import pickle

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

import support_scripts.Multiprocess as mp

def count_smiles(infile):
    with open(infile, 'rb') as handle:
        mols_dict = pickle.load(handle)

    printout = "{}:   {}".format(infile, len(list(mols_dict.keys())))
    return printout
    

folder = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/modified_pickled_lib/"
file_list = glob.glob(folder+"*")

job_input = [[file_list[i]] for i in range(0,len(file_list))]
print("     Multithreading mol_1 as variable, mol_2 as control")
output = mp.MultiThreading(job_input, -1,  count_smiles)
print("")
for x in output:
    print(x)



