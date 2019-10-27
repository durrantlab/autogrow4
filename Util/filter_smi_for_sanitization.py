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

def get_usable_list_smiles(infile, verbose=True):
    """
    given a file path, retrieve a list of smiles
    return a set of lists contain each molecules information
    ([SMILES_1,ZINCID_1],[SMILES_2,ZINCID_2],[SMILES_3,ZINCID_3])
    Also return the start count
    """

    zinc_name_list = []
    smile_list = []
    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    number_of_smiles_originally = 0
    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n","")
            parts = line.split('\t')      # split line into parts seperated by 4-spaces
            try:
                smile_list.append(parts[0])
                zinc_name_list.append(parts[1])              # store the 1st part as that

                number_of_smiles_originally = number_of_smiles_originally + 1 
            except:
                continue

    job_input = []
    for smile, name in zip(smile_list, zinc_name_list):
        job_input.append(([smile, name], verbose))

    return job_input, number_of_smiles_originally

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
    infile = sys.argv[1]
    outfolder = sys.argv[2]
    number_of_processors = sys.argv[3]
    verbose = False

    
    job_input, number_of_smiles_originally = get_usable_list_smiles(infile,verbose)

    print("Start multithread")

    # Run in Batches
    counter = 0
    output_all = {}
    batch_num = 0
    while counter < len(job_input):
        print("Counter is :",counter)
        print("batch_num is :",batch_num)
        job_batch = []
        tmp = []
        for x in range(counter, 1000+counter):
            if x > len(job_input)-1:
                continue
            job_batch.append(job_input[x])
        
        counter = counter +1000
        batch_num = batch_num +1
        # Run in Multithread for speed
        output = mp.MultiThreading(job_batch, number_of_processors, test_mol_sanitization)
        for x in output:
            if x is None:
                continue    
            output_all[x[0]] = x[1]
        
    print("Finished multithread")
    
    
    zinc_id_dict = {}
    for result in output:
        if result == None:
            continue

        else:
            mol_info = result[0]
            zinc_ID = mol_info[1]
            SMILES_Str = mol_info[0]

            zinc_id_dict[zinc_ID] = SMILES_Str

    # Get Count of how many SMILES failed to be imported and how many total input and how many pass
    number_smiles_pass_sanitize = len(list(zinc_id_dict.keys()))
    number_smiles_fail_sanitize = number_of_smiles_originally - number_smiles_pass_sanitize

    print("#######################################################################################")
    print("The number of entries in the initial file was: {}".format(number_of_smiles_originally))
    print("The number of Ligands which failed to Sanitize/protonate or redundant: {}".format(number_smiles_fail_sanitize))
    print("The number of Ligands which have substructures and passed Sanitize: {}".format(number_smiles_pass_sanitize))
    print("#######################################################################################")
    

    # Make a single files with all ligands
    output=""
    file_name = outfolder + 'concatinate_all_pass_sanitize' + ".smi"
    for zinc_ID in list(zinc_id_dict.keys()):
        SMILES_string = zinc_id_dict[zinc_ID]
        output = output + "{}\t{}\n".format(SMILES_string,zinc_ID)
    output = output + "\n"
        
    with open(file_name, 'w') as f:

        f.write(output)


    print("FINISHED")
