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

import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp
import run_filter as rf

import glob


import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp

import rdkit
import rdkit.Chem as Chem



def get_mols_dict_from_pickle(file_path):

    with open(file_path, 'rb') as handle:
        mols_dict = pickle.load(handle)
    return mols_dict

def write_pickle_to_file(file_path, obj):
    with open(file_path, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
        

def trim_a_dictionary_from_pickle_dict(base, end_size, folder_dir, outfolder_dir):
    
    input_pickle_file = folder_dir + base
    modified_pickle_file = outfolder_dir + base

    if os.path.exists(modified_pickle_file) == True:
        print("Already filtered:  {}".format(modified_pickle_file))
        return "Already filtered{}".format(base)


    print("")
    print(input_pickle_file)
    mols_dict = get_mols_dict_from_pickle(input_pickle_file)

    keys_mols_dict = list(set(mols_dict.keys()))
    len_mols_dict_start = len(keys_mols_dict)

    # shuffle the list of keys
    random.shuffle(keys_mols_dict)
    random.shuffle(keys_mols_dict)
    random.shuffle(keys_mols_dict)
    random.shuffle(keys_mols_dict)
    random.shuffle(keys_mols_dict)
    random.shuffle(keys_mols_dict)
    
    final_key_list = []
    for i in range(0,end_size):
        key = keys_mols_dict[i]
        final_key_list.append(key)
        key = None
    keys_mols_dict = None


    final_key_list = list(set(final_key_list))
    
    if len(final_key_list) != end_size:
        print(len(final_key_list))
        raise Exception("WRONG SIZE OF KEY LIST!!!!!!")


    new_dict = {}
    # Take the desired number of ligands
    for key in final_key_list:  
        new_dict[key] = mols_dict[key]

    mols_dict = None

    keys_new_dict = list(set(new_dict.keys()))
    len_new_dict = len(keys_new_dict)
    if len_new_dict != end_size:
        raise Exception("WRONG SIZE OF KEY LIST!!!!!!")

    # pickle the new dictionary

    write_pickle_to_file(modified_pickle_file, new_dict)
    new_dict = None

    printout = "{} sucessfully trimmer {}  ligands to {} ligands".format(base,len_mols_dict_start,len_new_dict)

    return printout

    
if __name__ == "__main__":

    end_size = 50000

    # Laptop
    over_100k =['aldehyde_or_ketone_pickle', 'alkyl_halogen_or_alcohol_pickle', 'amide_pickle', 'aryl_halide_pickle', 'carboxylic_acid_pickle', 'carboxylic_acid_or_ester_pickle', 'halide_pickle', 'hydrazine_pickle', 'ketone_pickle', 'nitrile_pickle', 'phenole_pickle', 'primary_amine_pickle', 'pyridine_pyrimidine_triazine_pickle', 'terminal_alkyne_pickle']

    folder_dir = "/mnt/data_1/DataB/jspiegel/FILTER_FOR_AUTO/source_pickled_lib/"
    folder_dir = "/mnt/data_1/DataB/jspiegel/FILTER_FOR_AUTO/Large_filter_one/"
    outfolder_dir = "/mnt/data_1/DataB/jspiegel/FILTER_FOR_AUTO/Larger_trimmed_50k/"
   
    if os.path.exists(outfolder_dir)== True:
        if os.path.isdir(outfolder_dir) == False:
            raise Exception("outfolder_dir is not a directory")
    else:
        os.mkdir(outfolder_dir)

    list_pickle_files = glob.glob(folder_dir + "*")
    to_trim_list_base_pickle_files= []
    for pickle_file in list_pickle_files:
        base = os.path.basename(pickle_file)
        
        if os.path.exists('/mnt/data_1/DataB/jspiegel/FILTER_FOR_AUTO/FinishedPICKLESRanThroughBothFilters/' + base) == True:
            continue
        
        # if os.path.exists('/mnt/data_1/DataB/jspiegel/FILTER_FOR_AUTO/Large_filter_one/' + base) == True:
        #     continue
        

        to_trim_list_base_pickle_files.append(base)



    jobs = [(x, end_size, folder_dir, outfolder_dir) for x in to_trim_list_base_pickle_files]
    
    print("list of jobs: ", jobs)
    print("Start multithread to_trim_list_base_pickle_files")
    # output = mp.MultiThreading(jobs, number_of_processors, run_filter_on_mol)   

    output = mp.MultiThreading(jobs, 5, trim_a_dictionary_from_pickle_dict)   

    print("")
    print("################")
    print("Finished Multithreading")
    print("################")
    print("")
    print("")
    for x in output:
        print(x)

    print("")

