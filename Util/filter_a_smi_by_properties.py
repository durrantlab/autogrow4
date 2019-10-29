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

import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp
import run_filter as rf

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
        





## 
def test_smiles_from_file(infile,new_file,number_of_processors,verbose):
    """
    given a file path, retrieve a list of smiles
    return a set of lists contain each molecules information
    ([SMILES_1,ZINCID_1],[SMILES_2,ZINCID_2],[SMILES_3,ZINCID_3])
    Also return the start count
    """
    
    # Run in Batches
    num_lines = 0
    with open(infile,'r') as smiles_file:
        for line in smiles_file.readlines():
            if line != '\n':

                num_lines = num_lines + 1


    counter = 0
    batch_num = 0
    number_of_passed_mols = 0
    none=0
    one=0
    two=0
    three=0
    four =0
    five=0
    six=0
    seven=0
    eight=0
    with open(new_file,'w') as nf:
        with open(infile,'r') as smiles_file:

            while counter < num_lines:
                job_batch = []
                output=[]
                print("")
                print("counter is :",counter)
                print("batch_num is :",batch_num)
                for x in range(counter, counter+10000):
                    if counter > num_lines-1:
                        continue

                    line = smiles_file.readline()
                    if line == "\n":
                        break
                    line = line.replace("\n","")
                    parts = line.split('\t')      # split line into parts seperated by 4-spaces
                    if len(parts) != 2:
                        continue
                    
                    job_batch.append((parts[0],parts[1]))


                print("Start multithread")
                output = mp.MultiThreading(job_batch, number_of_processors, run_filter)    
                job_batch = []

                nones = [x for x in output if x is None]
                ones = [x for x in output if x is 1]
                twos =  [x for x in output if x is 2]
                threes=  [x for x in output if x is 3]
                fours= [x for x in output if x is 4]
                fives= [x for x in output if x is 5]
                sixes= [x for x in output if x is 6]
                sevens= [x for x in output if x is 7]
                eights= [x for x in output if x is 8]
                nines= [x for x in output if x is 9]



                #################################
                none = none+ len(nones)
                one=one+ len(ones)
                two= two+ len(twos)
                three= three + len(threes)
                four= four + len(fours)
                five= five + len(fives)
                six= six + len(sixes)
                seven= seven + len(sevens)
                eight= eight + len(eights)
                #################################
   

                output = [(x) for x in output if x is not None and type(x) is not int]
  
                counter = counter + 10000
                batch_num = batch_num +1
                
                number_of_passed_mols = number_of_passed_mols + len(output)



                #write to a file and delete unnecessary objects to minimize mem
                for x in output:
                    smile = str(x[0])
                    zincId = str(x[1])
                    save_line = "".join([smile,"\t",zincId,"\n"])
                    nf.write(save_line)

                output=[]
    
    print("#################################")

    print("Nones: ", none)
    print("ones: ", one)
    print("twos: ", two)
    print("threes: ", three)
    print("fours: ", four)
    print("fives: ", five)
    print("sixes: ", six)
    print("sevens: ", seven)
    print("eights: ", eight)
    print("#################################")
    return number_of_passed_mols, num_lines

def run_on_pickled_dictionarys(pickle_infile,modified_pickle_file):

    mols_dict = get_mols_dict_from_pickle(pickle_infile)
    mol_dict_keys = list(mols_dict.keys())
    start_len = len(mol_dict_keys)
    counter = 0
    batch_num = 0
    number_of_passed_mols = 0
    passed_zinc_ids = []
    use_lienent = False
    while counter < start_len:
        job_batch = []
        output=[]
        print("")
        print("counter is :",counter)
        print("batch_num is :",batch_num)
        print("Total Number to do is : ", start_len)
        for x in range(counter, counter+10000):
            if x > start_len-1:
                continue
            zincID = mol_dict_keys[x]
            job_batch.append(mols_dict[zincID])
        # print(job_batch)

        print("Start multithread run_filter_on_mol")
        # output = mp.MultiThreading(job_batch, number_of_processors, run_filter_on_mol)   

        output = mp.MultiThreading(job_batch, number_of_processors, rf.run_normal_filter)   

        output = [x for x in output if x is not None]

        pass_list = [x[1] for x in output if x[0]==0]
        if len(pass_list) < 100 and len(job_batch) > 1000:

            pass_list = [x[1] for x in output if x[0]==0 or x[0]==1]
            if len(pass_list) < 300 and start_len > 1000:
                print("LESS THAN 300 in 1st Filter. May be some intrisic property.")
        output = None
        job_batch=[]
        for x in pass_list:
            job_batch.append(mols_dict[x])
        
        if use_lienent == False:
            output = mp.MultiThreading(job_batch, number_of_processors, rf.run_filter_on_mol_substructurefilters_strict)   

            output = [x for x in output if x is not None]
            if len(output)==0:
                
                use_lienent = True
              
        if use_lienent == True:
            output = mp.MultiThreading(job_batch, number_of_processors, rf.run_filter_on_mol_substructurefilters_lienent)   

        output = [x for x in output if x is not None]
        if len(output)==0:
            print("LIENENT FILTER REMOVED ALL LIGANDS!!!!")

        pass_list = [x for x in output if x is not None]
        if len(pass_list) < 100 and len(job_batch) > 1000:
            print("LESS THAN 300 in second Filter. May be some intrisic property.")
        


        counter = counter + 10000
        batch_num = batch_num +1
        number_of_passed_mols = number_of_passed_mols + len(output)
        passed_zinc_ids.extend(output)
    
    modified_pickle = pickle_infile +"filtered_pickle"
    if number_of_passed_mols < 300:
        print("LESS THAN 300 for {}".format(modified_pickle))
        print("Total count was {}".format(number_of_passed_mols))

    new_dict = {}
    for x in passed_zinc_ids:
        new_dict[x] = mols_dict[x]
    filtered_len = start_len - len(passed_zinc_ids)
    mols_dict= None
    
    print("Pickling modified dict")
    write_pickle_to_file(modified_pickle_file, new_dict)
    print("")
    print("")
    print("###########################")
    print("Started with :   ",start_len)
    print("Filtered   :   ",filtered_len)
    print("Ended with   :   ",len(passed_zinc_ids))
    print("###########################")
    print("")



def test_mol_sanitization(mol_info, verbose=False):
    """
    mol_info: a list containing the smiles_string and the zinc_id
                (smile_str, zinc_ID)
    returns
    mol_info:  mol_info     SMILES May be modified when Sanitizing or fixing Nitrogen valences
    returns None if it Fails
    """
    
    mol_info_copy = copy.deepcopy(mol_info)

    mol = Chem.MolFromSmiles(mol_info_copy[0], sanitize=False)

    # Confirm the SMILES of the molecule can be imported into RDKIT Sanitized and protanate/deprotanate 
    if type(mol) != rdkit.Chem.rdchem.Mol:
        if verbose == True:
            print("mol Failed to import into rdkit {}".format(mol_info_copy))
        return None
    mol = MOH.check_sanitization(mol)
    if type(mol) != rdkit.Chem.rdchem.Mol:
        if verbose == True:
            print("mol Failed to sanitize into rdkit {}".format(mol_info_copy))
        return None

    # Confirm it can protanate and Deprotanate (this is important in Crossover and mutation.)
    mol = MOH.handleHs(mol, True)
    if type(mol) != rdkit.Chem.rdchem.Mol:
        if verbose == True:
            print("mol Failed to handle protanation into rdkit {}".format(mol_info_copy))
        return None

    mol = MOH.try_deprotanation(mol)
    if type(mol) != rdkit.Chem.rdchem.Mol:
        if verbose == True:
            print("mol Failed to handle protanation into rdkit {}".format(mol_info_copy))
        return None

    # This is a Good Mol !!!

    # Test if the SMILES string needed correcting
    SMILES_string = Chem.MolToSmiles(mol,isomericSmiles=True,canonical=True)
    if SMILES_string != mol_info[0]:
        if verbose == True:
            printout = "SMILES string is modified for sanitization purpose for the Ligand: {}\n".format(mol_info[1])
            printout = printout + "{} Changed to ==> {}".format(mol_info[0], SMILES_string)
            print(printout)

    mol_info = (SMILES_string, mol_info[1])
    

    return mol_info

  
if __name__ == "__main__":


    # Laptop
    infile = "/home/jacob/Downloads/zinc15_available/pass_sanitize_all.smi_no_dup.smi"
    new_file = "/home/jacob/Downloads/zinc15_available/pass_filters.smi"
    outfolder = "/home/jacob/Downloads/zinc15_available/"
    number_of_processors = -1
    print("the number of processors is:  ",number_of_processors)
    

    # # Run on Barry 15million
    # number_of_processors = -1
    # infile = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/sanitize_no_dup.smi"
    # new_file = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/pass_filters.smi"
    # outfolder = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/"


    # For Bob
    folder_dir = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/Larger_trimmed_50k/"
    folder_dir ="/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/filtered_props/"
    outfolder_dir = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/Large_filter_one/"
    # For CRC
    # folder_dir = "/ihome/jdurrant/jspiegel/FILTER_FOR_AUTO/"

    list_pickle_files = glob.glob(folder_dir + "*")

    # for base in over_100k:
    for pickle_infile in list_pickle_files:
        # pickle_infile = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/Filtered_largest_files_pickled/" + base
        # # pickle_infile = folder_dir + base
        # if "aldehyde_or_ketone_pickle" in pickle_infile:
        #     print("aldehyde_or_ketone_pickle")
        #     continue

        base = os.path.basename(pickle_infile)
        modified_pickle_file = outfolder_dir + base
        if os.path.exists(modified_pickle_file) == True:
            print("Already filtered:  {}".format(modified_pickle_file))
            continue
        print("")
        print(pickle_infile)
        run_on_pickled_dictionarys(pickle_infile,modified_pickle_file)
    
        print(base)



