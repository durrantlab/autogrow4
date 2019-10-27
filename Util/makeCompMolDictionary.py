# script for converting large smi file into files for complimentary dictionary


import __future__

import rdkit
from rdkit import Chem
import copy
import json
import sys
import os
import pickle

import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp
    
def create_blank_fun_files(fun_group_dict):
    file_name_dict = {}
    for fun_group_name in list(fun_group_dict.keys()):
        file_name = outfolder + fun_group_name + ".smi"

        file_name_dict[fun_group_name] = file_name

        if os.path.exists(file_name):

            print("Deleting the following files: {}".format(file_name))
            os.remove(file_name)

    return file_name_dict

def write_pickle_to_file(file_path, obj):
    with open(file_path, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def get_mols_list_from_pickle(file_path):

    with open(file_path, 'rb') as handle:
        mols_list = pickle.load(handle)
    return mols_list

def make_mol_list_dataset(infile):

    num_lines = 0
    with open(infile,'r') as smiles_file:
        for line in smiles_file.readlines():
            if line != '\n':

                num_lines = num_lines + 1


    # Run in Batches
    counter = 0
    batch_num = 0

    with open(infile,'r') as smiles_file:
        mols_list_full = []
        while counter < num_lines:
            job_batch = []
            print("")
            print("counter is :",counter)
            print("batch_num is :",batch_num)
            for x in range(counter, counter+100000):
                if counter > num_lines-1:
                    continue

                line = smiles_file.readline()
                if line == "\n":
                    break
                line = line.replace("\n","")
                parts = line.split('\t')      # split line into parts seperated by 4-spaces
                if len(parts) != 2:

                    parts = line.split('    ')      # split line into parts seperated by 4-spaces                
                    if len(parts) != 2:

                        continue

                tmp = [(parts[0], parts[1])]
                job_batch.append((tmp))
                tmp=[]

            print("Start multithread protanated_mols")
            protanated_mols = mp.MultiThreading(job_batch, number_of_processors, create_mol)    
            job_batch = []
            protanated_mols = [x for x in protanated_mols if x is not None]
            print("Fin multithread protanated_mols")
            mols_list_full.append(protanated_mols)
            # Handle the counters
            counter = counter + 100000
            batch_num = batch_num +1
        #create a pickled version 
        pickle_file = infile + "_pickle"
        write_pickle_to_file(pickle_file, mols_list_full)
        return mols_list_full

def test_smiles_from_file(infile, fun_group_dict, number_of_processors):
    """
    given a file path, retrieve a list of smiles
    return a set of lists contain each molecules information
    ([SMILES_1,ZINCID_1],[SMILES_2,ZINCID_2],[SMILES_3,ZINCID_3])
    Also return the start count
    """
    
    file_name_dict = create_blank_fun_files(fun_group_dict)
    
    #Check if file has already been made and pickled
    # This is a list of molecules which have been sanitized and protanated
    pickle_file = infile + "_pickle"

    if os.path.exists(pickle_file)==True:
        mols_list_full = get_mols_list_from_pickle(pickle_file)
    else:
        mols_list_full = make_mol_list_dataset(infile)
 
    print("")
    print("OBTAINED MOLECULES TO USE")
    print("")

    # Run in Batches
    counter = 0
    batch_num = 0
    number_of_mols_w_fun_groups = 0
    num_mols = len(mols_list_full)
    with open(infile,'r') as smiles_file:
        
        while counter < 100000:
            job_batch = []
            print("")
            print("counter is :",counter)
            print("batch_num is :",batch_num)
            protanated_mols = []
            for x in range(counter, counter+100000):
                if x > num_mols-1:
                    continue

                # find_fun_groups(SMILES_string, ZINCID, protanated_mol, fun_group_dict)
                temp = (mols_list_full[x][2][0],mols_list_full[x][1],mols_list_full[x][2][2], fun_group_dict)
                job_batch.append(temp)
               
            print("Start multithread test_fun_groups")
            mol_fun_group_results = mp.MultiThreading(job_batch, number_of_processors, find_fun_groups)    
            job_batch = [] # Clear for memory reasons
            mol_fun_group_results = [x for x in mol_fun_group_results if x is not None]
            print("Fin multithread test_fun_groups")
            
            # Handle the counters
            number_of_mols_w_fun_groups = number_of_mols_w_fun_groups + len(mol_fun_group_results)
            counter = counter + 100000
            batch_num = batch_num +1

            fun_mol_dictionary, zinc_id_dictionary = convert_fungroup_results_to_dict(mol_fun_group_results)
            mol_fun_group_results = [] # Clear for memory reasons

            # append the ligands to the fun_group files
            job_input = []
            for fun_group_name in list(fun_mol_dictionary.keys()):
                filename = file_name_dict[fun_group_name]
                list_of_mols = fun_mol_dictionary[fun_group_name]
                tmp = (fun_group_name, filename, list_of_mols, zinc_id_dictionary)
                job_input.append(tmp)
                tmp = () # Clear for memory reasons

            fun_mol_dictionary = {} # Clear for memory reasons
            zinc_id_dictionary = {} # Clear for memory reasons
            
            # Append to lists
            print("Start multithread append_smiles_to_fun_group_file")
            mp.MultiThreading(job_input, number_of_processors, append_smiles_to_fun_group_file) 
            print("Fin multithread append_smiles_to_fun_group_file")
            
            job_input = [] # Clear for memory reasons

    
    return number_of_mols_w_fun_groups, num_mols, file_name_dict


def append_smiles_to_fun_group_file(fun_group_name, filename, list_of_mols, zinc_id_dictionary):

    with open(filename, 'a') as f:    
        for zinc_id in list_of_mols:
            
            smiles_str = zinc_id_dictionary[zinc_id]

            save_line = "".join([smiles_str,"\t",zinc_id,"\n"])

            f.write(save_line)
            

def find_fun_groups(SMILES_string, ZINCID, protanated_mol, fun_group_dict):
    """
    mol_info:[(SMILES_string, ZINCID, protanated rdkit.Chem.rdchem.Mol)]
    """ 
    fun_groups_in_mol = []

    for fun_group_name in list(fun_group_dict.keys()):

        substructure = fun_group_dict[fun_group_name]

        protanated_mol = MOH.try_reprotanation(protanated_mol)
        if protanated_mol == None:
            return None


        if protanated_mol.HasSubstructMatch(substructure) == True:
            fun_groups_in_mol.append(fun_group_name)
        else:
            protanated_mol = MOH.try_deprotanation(protanated_mol)
            if protanated_mol == None:
                return None

            if protanated_mol.HasSubstructMatch(substructure) == True:
                fun_groups_in_mol.append(fun_group_name)
             
        substructure = None


    if len(fun_groups_in_mol) == 0:
        #print("No subgroups found within {}".format(ZINCID))
        return None
    else:
        return [SMILES_string, ZINCID], fun_groups_in_mol


def create_mol(mol_info):
    """
    mol_info: a list containing the smiles_string and the zinc_id
                (smile_str, zinc_ID)
    returns
    mol_info[0][0], mol_info[0][1], mol
    returns None if it Fails
    """

    try:
        mol = Chem.MolFromSmiles(mol_info[0][0])
    except:
        raise Exception("A MOL FAILED TO SANITIZE!!!!")
    # Confirm the SMILES of the molecule can be imported into RDKIT Sanitized and protanate/deprotanate 
    if type(mol) != rdkit.Chem.rdchem.Mol:
        raise Exception("A MOL FAILED TO SANITIZE!!!!")

    mol = MOH.try_reprotanation(mol)

    mol = Chem.AddHs(mol)

    if type(mol) != rdkit.Chem.rdchem.Mol:
        raise Exception("A MOL FAILED TO SANITIZE!!!!")


    mol = Chem.AddHs(mol)

    return mol_info[0][0], mol_info[0][1], mol


def retrieve_reaction_dict(functional_group_json_file):


    with open(functional_group_json_file, 'r') as rxn_file:
        reaction_dict = json.load(rxn_file)

    reaction_dict = rxn_lib_format_json_dict_of_dict(reaction_dict)

    new_dict = {}
    for group in reaction_dict:
        smart = reaction_dict[group]
        substructure = Chem.MolFromSmarts(smart)
        if substructure == None:
            print("substructure FAILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print(group)
            print(smart)
            raise Exception("SMART failed to import into rdkit. check the smart {} : {}".format(group,smart))
        new_dict[group] = substructure


    return new_dict


def rxn_lib_format_json_dict_of_dict(old_dict):
    """
    json dictionaries  import as type unicode. This script converts all the keys and items to strings, with a few specific exceptions.
    It takes both the functional group dictionary and the reaction library. 
    
    The reaction library is a dictionary of dictionary and has a few exceptions which are not intended to be strings.
        ie. the num_reactants which converts to interger and functional_groups which convert to a list of strings. 

    The functional_group_dictionary is simply a dictionary with all items and keys needing to be strings.

    :param dic old_dict: a dictionary of the the reaction library or functional groups. This is what is importated from the .json file.
    returns:
    :returns: dic new_dict: a dictionary of the the reaction library or functional groups where the unicode type items have been replaced with
                    the proper python data types.
    """
    new_dict = {}
    
    for rxn_key in list(old_dict.keys()):
        rxn_dic_old = old_dict[rxn_key]
        key_str = str(rxn_key)
        
        # For reaction libraries
        if type(rxn_dic_old) == dict:
            new_sub_dict = {}
            for key in rxn_dic_old.keys():
                sub_key_str = str(key)
                item = rxn_dic_old[key]

                if sub_key_str == "num_reactants":
                    item = int(item)
                elif  sub_key_str == "functional_groups":
                    new_list = []
                    for i in item:
                        i_str = str(i)
                        new_list.append(i_str)

                    item = new_list
                else:
                    item = str(item)
                
                new_sub_dict[sub_key_str] = item
            new_dict[key_str] = new_sub_dict

        # For functional groups
        else:
            item = old_dict[rxn_key]
            new_dict[key_str] = str(item)
    
    return new_dict
#
def convert_fungroup_results_to_dict(mol_fun_group_results):

    zinc_id_dict = {}

    new_dict = {}
    for mol in mol_fun_group_results:
        smile = mol[0][0]
        zinc_Id = mol[0][1]

        zinc_id_dict[zinc_Id] = smile
        new_dict [zinc_Id] = mol[1]

    # invert new_dict
    fun_dict = invert_dictionary(new_dict)

    return fun_dict, zinc_id_dict


def invert_dictionary(old_dic):
    """
    This will invert any dictionary so that the keys are the values and the values are the keys.
    Inputs: 
    :param dict old_dic: a dictionary to invert 
    Return: 
    :returns: dict inverted_dic: old_dict dict inverted so the keys are the items and the items are the keys
    """
    
    # inverted_dic = {}
    # for k, v in old_dic.iteritems():
    # keys = inverted_dic.setdefault(v, [])
    # keys.append(k)
    values = set([a for b in list(old_dic.values()) for a in b])
    values = list(values)
    inverted_dic = dict((new_key, [key for key, value in list(old_dic.items()) if new_key in value]) for new_key in values)
     
    return inverted_dic
#
  

if __name__ == "__main__":
    # infile = sys.argv[1]
    # outfolder = sys.argv[2]
    # number_of_processors = -1
    # functional_group_json_file = sys.argv[3]

    # functional_group_json_file = "/home/jspiegel/DataB/jspiegel/projects/autogrow4/autogrow/Operators/Mutation/SmileClickChem/Reaction_libraries/Robust_Rxns/Robust_Rxns_functional_groups.json"
    # infile = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_available/pass_sanitize_all.smi_no_dup.smi"
    # outfolder = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_available/compound_lib/"
    # number_of_processors = -1
    infile = "/home/jacob/Downloads/zinc15_available/pass_sanitize_all.smi_no_dup.smi"
    outfolder = "/home/jacob/Downloads/zinc15_available/compound_lib/"
    number_of_processors = -1
    #LAptop
    functional_group_json_file = "/home/jacob/Documents/autogrow4/autogrow/Operators/Mutation/SmileClickChem/Reaction_libraries/Robust_Rxns/Robust_Rxns_functional_groups.json"
    infile = "/home/jacob/Downloads/zinc15_available/concatinated_all.smi_no_dup.smi"
    outfolder = "/home/jacob/Downloads/zinc15_available/"
    number_of_processors = -1

    #Bob
    functional_group_json_file = "/home/jspiegel/DataB/jspiegel/projects/autogrow4/autogrow/Operators/Mutation/SmileClickChem/Reaction_libraries/Robust_Rxns/Robust_Rxns_functional_groups.json"
    infile = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_data/zinc15_data/concatinated_all.smi"
    outfolder = "/home/jspiegel/DataB/jspiegel/zinc15_data/compound_lib/"
    number_of_processors = 22
    #Bob
    functional_group_json_file = "/mnt/data_1/DataB/jspiegel/projects/autogrow/smileConversionTools/temp/reaction.json"
    infile = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_data/pass_filters.smi"
    outfolder = "/mnt/data_1/DataB/jspiegel/projects/autogrow/smileConversionTools/temp/"
    number_of_processors = -1
    

    functional_group_json_file = "/home/jspiegel/Desktop/temp/fungroup.json"
    infile =  "/home/jspiegel/DataB/jspiegel/compound_library_robust_rxns/small_Zinc_subset/zinc15_data/pass_filters.smi"
    outfolder = "/home/jspiegel/Desktop/temp/"
    number_of_processors = -1
    


    # # Run on Barry 15million
    # number_of_processors = 21
    # functional_group_json_file = "/mnt/c/Users/kcc35/Desktop/JAKE/autogrow4/autogrow/Operators/Mutation/SmileClickChem/Reaction_libraries/Robust_Rxns/Robust_Rxns_functional_groups.json"
    # infile = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/pass_filters.smi"
    # outfolder = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/compound_lib/"


    # functional_group_json_file = "/home/jspiegel/DataB/jspiegel/projects/autogrow4/autogrow/Operators/Mutation/SmileClickChem/Reaction_libraries/Robust_Rxns/Robust_Rxns_functional_groups.json"
    # infile = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_available/pass_sanitize_all.smi_no_dup.smi"
    # outfolder = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_available/compound_lib/"
    # number_of_processors = -1
    
    
    # functional_group_json_file = "/home/jspiegel/Desktop/tmp/fake_fun_group_lib.json"
    # outfolder = "/home/jspiegel/Desktop/tmp/"

    fun_group_dict = retrieve_reaction_dict(functional_group_json_file)

    print("TEST_SMILES_FROM_FILE")

    number_of_mols_w_fun_groups, num_lines, file_name_dict = test_smiles_from_file(infile, fun_group_dict, number_of_processors)

    if num_lines == 0:
        raise Exception("There were no ligands in the infile. Please check that it is a tab deliniated .smi file.")
    
    if number_of_mols_w_fun_groups == 0:
        raise Exception("None of the ligands in the files contained the searched for functional groups. Please check the functional groups and the infile of smiles.")
    


    mols_without_fun_groups =  num_lines - number_of_mols_w_fun_groups
    list_of_fun_groups = list(file_name_dict.keys())

    # Determine if any fun_group were empty. If no molecules were found with a given functional group
    #       The file in the file_name_dict will not exist so check os.path.exists(filename)
    list_of_empty_fun_groups = []
    for fun_group_name in list_of_fun_groups:
        file_path = file_name_dict[fun_group_name]
        if os.path.exists(file_path) == False:
            list_of_empty_fun_groups.append(fun_group_name)

    num_of_populated_fungroup_files = len(list_of_fun_groups)- len(list_of_empty_fun_groups)
    print("")
    print("")
    print("#######################################################################################")
    print("The number of entries in the initial file was: {}".format(num_lines))
    print("The number of Ligands with no found functional groups is: {}".format(mols_without_fun_groups))
    print("The number of Ligands with functional groups is: {}".format(number_of_mols_w_fun_groups))
    print("#######################################################################################")
    print("The number of Functional groups within the dictionary is: {}".format(len(list_of_fun_groups)))
    print("The number of Functional groups which were populated is: {}".format(num_of_populated_fungroup_files))
    print("The number of EMPTY Functional groups is: {}".format(len(list_of_empty_fun_groups)))
    print("#######################################################################################")
    print("")
    print("")

    if len(list_of_empty_fun_groups) != 0:
        printout = "We failed to find any ligands containing the following functional groups:\n"
        printout =  "{}(Autogrow will not run without these functional group librarys)\n".format(printout)
        for fun_group_name in list_of_empty_fun_groups:
            printout =  "{}\t{}\n".format(printout, fun_group_name)

        print(printout)
        print("")
        print("")
        raise Exception(printout)
    
print("")
print("Ready to run in Autogrow!!!!!")
print("FINISHED")
