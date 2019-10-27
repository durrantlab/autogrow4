# This script should convert a list of ligand identifiers, this will compile the information about the ligand
# This will take the original info from the unranked smi file which created the ligand. This has no docking information. 

import glob
import sys
import os

file_path_ligands = "/mnt/data_1/DataB/jspiegel/projects/autogrow_testing/2018-05-07/timeout_list.txt"
run_folder = "/home/jspiegel/projects/AUTOGROW_PARP/Run_1"
#all_folders_list = [f for f in sorted(os.listdir(infolder)) if os.path.isdir(infolder+f)]


def get_list_of_ligands(file_path_ligands):
    list_of_ligands = []
    with open(file_path_ligands, 'r') as f:
        for line in f.readlines():
            line = line.replace("\t","")
            line = line.replace("\n","")
            parts = line.split(', ',)      # split line into parts seperated by ", "
            parts = [x for x in parts if x!=None]
            
            list_of_ligands.extend(parts)
    return list_of_ligands

def find_generation_of_ligand(ligand_shorthand,run_folder):
    """
    ligand_shorthand is a string
    """

    gen_num = ligand_shorthand.split("_")[1]
    print(gen_num)




if __name__ == "__main__":
    list_of_ligands = get_list_of_ligands(file_path_ligands)
    find_generation_of_ligand(list_of_ligands[0],run_folder)

