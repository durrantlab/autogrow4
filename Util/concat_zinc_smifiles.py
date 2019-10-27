# Concatanate all .smis in zinc folder with sub folder


import glob
import os
import sys
import support_scripts.Multiprocess as mp

def get_all_mols(rank_file):
    counter = 0
    data_dict = {}
    num_redundancy_removed = 0

    with open(rank_file, 'r') as rf:
        rf.readline()
        for line in rf.readlines():
            if line == "\n":
                continue
                
            line = line.replace("\n","")
            line = line.replace("\r","")
            parts = line.split(" ")      # split line into parts seperated by 4-spaces
            
            try:
                smile = parts[0]
                name = parts[1]
            
            except:
                continue

            counter = counter + 1
            data_dict[name] = smile
    num_redundancy_removed = len(list(data_dict.keys())) - counter
    if num_redundancy_removed != 0:
        print("counter =", counter)
        print("len(list(data_dict.keys())) =", len(list(data_dict.keys())))
        print("num_redundancy_removed =", num_redundancy_removed)
    return data_dict, num_redundancy_removed

#########################
def get_sub_folders(folder):

    sub_folders_list = []
    sub_sub_folders_list = []

    file_list = []
    for sub_folders in glob.glob(folder + '*'):
        if os.path.isdir(sub_folders):
            sub_folders_list.append(sub_folders)

    for sub_folders in sub_folders_list:
            for sub_sub in glob.glob(sub_folders + '/*'):

                if os.path.isfile(sub_sub):
                    file_list.append(sub_sub)
                elif os.path.isdir(sub_sub):
                    sub_sub_folders_list.append(sub_folders)
                    raise Exception("Multiple layers!!!!!")

    return file_list

def compile_and_Remove_Redundancy(list_of_dicts):


    full_dic = {}
    num_redundancy_removed = 0
    counter = 0
    for results in list_of_dicts:
        dictionary = results[0]
        redundant_removed = results[1]
        num_redundancy_removed = num_redundancy_removed + redundant_removed

        counter =counter + len(list(dictionary.keys()))
        for zinc_id in list(dictionary.keys()):
            
            full_dic[zinc_id] = dictionary[zinc_id]



    start_num_lig = counter + num_redundancy_removed
    num_redundancy_removed = start_num_lig -  len(list(full_dic.keys()))
    print("The starting number of files was: {}".format(len(list_of_dicts)))
    print("The starting number of ligands was: {}".format(start_num_lig))
    print("The Number of ligands removed for redundant Zinc IDs was: {}".format(num_redundancy_removed))

    reduced_dic = {}
    for zinc_id in list(full_dic.keys()):
        smile_str = full_dic[zinc_id]
        reduced_dic[smile_str] = zinc_id

    removed_from_redundant_SMILES = len(list(full_dic.keys())) - len(list(reduced_dic.keys())) 

    print("The Number of ligands removed for redundant SMILES was: {}".format(removed_from_redundant_SMILES))
    print("#############################################################")
    redundant_removed_count = removed_from_redundant_SMILES + num_redundancy_removed
    if start_num_lig-len(list(reduced_dic.keys())) != redundant_removed_count:
        print("FAIL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("Total Number of Ligands removed for redundancy is: {}".format(redundant_removed_count))
    print("Total Number of unique Ligands is: {}".format(len(list(reduced_dic.keys()))))

    return reduced_dic


#############################################################

folder = str(sys.argv[1])

output_file = str(sys.argv[2])

file_list = get_sub_folders(folder)

print("")
print("Finished getting files!")
print("")

job_input = [[i] for i in file_list]
output_list = mp.MultiThreading(job_input, -1,  get_all_mols)

print("")
print("finished getting Mols!")
print("")

# Remove redundant SMILES Strings
reduced_dic = compile_and_Remove_Redundancy(output_list)


print("")
print("Finished Compiling Mol dictionary!")
print("")



with open(output_file,'w') as f:
    
    for smiles in list(reduced_dic.keys()):
        tmp_str = "{}\t{}".format(smiles, reduced_dic[smiles])

        f.write(tmp_str + "\n")



print("")
print("Finished writing!")
print("")
