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

import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp


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

                    job_batch.append(((parts[0],parts[1]),verbose))


                print("Start multithread")
                output = mp.MultiThreading(job_batch, number_of_processors, test_mol_sanitization)    
                job_batch = []
                output = [(x) for x in output if x is not None]
                
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
    
    return number_of_passed_mols, num_lines
    

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
    # time python filter_smi_for_sanitization.py \
    # /home/jacob/Downloads/zinc15_available/concatinate_all_pass_sanitize.smi \
    # /home/jacob/Downloads/zinc15_available/compound_lib -1
    # infile = sys.argv[1]
    # outfolder = sys.argv[2]
    # number_of_processors = sys.argv[3]
    verbose = False

    #LAptop
    infile = "/home/jacob/Downloads/zinc15_available/concatinated_all.smi_no_dup.smi"
    new_file = "/home/jacob/Downloads/zinc15_available/pass_sanitize_all_no_dup2.smi"
    outfolder = "/home/jacob/Downloads/zinc15_available/"
    number_of_processors = -1

    #Bob
    infile = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_data/zinc15_data/concatinated_all.smi"
    new_file = "/home/jspiegel/DataB/jspiegel/small_Zinc_subset/zinc15_data/zinc15_data/pass_sanitize_all_no_dup.smi"
    outfolder = "/home/jacob/Downloads/zinc15_available/"
    number_of_processors = 22



    # Run on Barry 15million
    # number_of_processors = 21
    # infile = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/concatinated_all.smi_no_dup.smi "
    # new_file = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/pass_sanitize_all_no_dup.smi"
    # outfolder = "/mnt/c/Users/kcc35/Desktop/JAKE/zinc15_data/"


    # new_file = outfolder + 'concatinate_all_pass_sanitize_all2' + ".smi"


    number_smiles_pass_sanitize, number_of_smiles_originally = test_smiles_from_file(infile,new_file,number_of_processors,verbose)

    print("Finished multithread")
    print("")
    number_smiles_fail_sanitize = number_of_smiles_originally - number_smiles_pass_sanitize
    print("#######################################################################################")
    print("The number of entries in the initial file was: {}".format(number_of_smiles_originally))
    print("The number of Ligands which failed to Sanitize/protonate: {}".format(number_smiles_fail_sanitize))
    print("The number of Ligands which have substructures and passed Sanitize: {}".format(number_smiles_pass_sanitize))
    print("#######################################################################################")
    

    print("FINISHED")
