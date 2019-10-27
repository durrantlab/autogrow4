import rdkit
import glob
import os
import sys
import random

folder = "/home/jacob/Documents/autogrow4/autogrow/Operators/Mutation/SmileClickChem/Reaction_libraries/Robust_Rxns/complimentary_mol_dir/"

folder2 = "/home/jacob/Documents/autogrow4/autogrow/Operators/Mutation/SmileClickChem/Reaction_libraries/ClickChem/complimentary_mol_dir/"


output_file = "/home/jacob/Documents/autogrow/source_compounds/ZINC_fragments.smi"


def handle_one_file(filename):
    line_count = 0
    temp_list = []
    final_list = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            line_count = line_count+1
            temp_list.append(line)


    if line_count < 100:
        for line in temp_list:
            final_list.append(line)
    
    else:
        random.shuffle(temp_list)
        random.shuffle(temp_list)

        for x in range(0,100):
            final_list.append(temp_list[x])
    
    return final_list


file_list = glob.glob(folder + "*.smi") + glob.glob(folder2 + "*.smi")
print(file_list)

full_list = []
for filename in file_list:
    temp_list = handle_one_file(filename)
    full_list.extend(temp_list)

random.shuffle(full_list)
random.shuffle(full_list)
random.shuffle(full_list)
random.shuffle(full_list)
print("")
print("FINISHED")
print("LEN is: ", len(full_list))

printout=""
for line in full_list:
    printout = printout + line

with open(output_file,"w") as f:
    f.write(printout)



