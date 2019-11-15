import __future__

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #

import support_scripts.Multiprocess as mp

def get_usable_fomat(infile):
    """
    This code takes a string for an file which is formated as an .smi file. It
    opens the file and reads in the components into a usable list.

    The .smi must follow the following format for each line:
        MANDATORY INFO
            part 1 is the SMILES string
            part 2 is the SMILES name/ID

        Optional info
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating
                - For default setting this is the Docking score
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)

    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_list_of_smiles: list of SMILES and their associated
        information formated into a list which is usable by the rest of Autogrow
    """

    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_list_of_smiles = []

    if os.path.exists(infile) is False:
        print("\nFile of Source compounds does not exist: {}\n".format(infile))
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts seperated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts seperated by 4-spaces

            choice_list = []
            for i in range(0, len(parts)):
                choice_list.append(parts[i])
            usable_list_of_smiles.append(choice_list)

    return usable_list_of_smiles


def add_mol_to_list(usable_list_line):

    try:
        mol = Chem.MolFromSmiles(usable_list_line[0])
        usable_list_line.append(mol)
    except:
        usable_list_line = None
    return usable_list_line

def mcs_mols(two_mols_list, list_for_mol2):

    mcs_results = rdFMCS.FindMCS(two_mols_list, matchValences=False, ringMatchesRingOnly=True, completeRingsOnly=True)
    num_mcs = str(mcs_results.numAtoms)
    list_for_mol2.append(num_mcs)

    return list_for_mol2

def get_PARPis():
    olaparib = Chem.MolFromSmiles("C1CC1C(=O)N2CCN(CC2)C(=O)C3=C(C=CC(=C3)CC4=NNC(=O)C5=CC=CC=C54)F")
    rucaparib = Chem.MolFromSmiles("CNCC1=CC=C(C=C1)C2=C3CCNC(=O)C4=CC(=CC(=C34)N2)F")
    niraparib = Chem.MolFromSmiles("C1CC(CNC1)C2=CC=C(C=C2)N3C=C4C=CC=C(C4=N3)C(=O)N")
    veliparib = Chem.MolFromSmiles("CC1(CCCN1)C2=NC3=C(C=CC=C3N2)C(=O)N")
    iniparib = Chem.MolFromSmiles("C1=CC(=C(C=C1C(=O)N)[N+](=O)[O-])I")
    talazoparib = Chem.MolFromSmiles("CN1C(=NC=N1)C2C(N=C3C=C(C=C4C3=C2NNC4=O)F)C5=CC=C(C=C5)F")
    parpis= [olaparib, rucaparib, niraparib, veliparib, iniparib, talazoparib]
    parpis_names = ["olaparib", "rucaparib", "niraparib", "veliparib", "iniparib", "talazoparib"]

    #Change the index number in this list to visualize the mol
    return parpis,parpis_names


if __name__ == "__main__":

    input_smi = "/home/jacob/Desktop/Outputfolder/Trial_dock/LGMNB_PAINS.smi"
    usable_list = get_usable_fomat(input_smi)
    new_file = "/home/jacob/Desktop/Outputfolder/Trial_dock/ranked_PARPi_comparison.smi"

    # Convert Strings to RDKIT mol objects and append to end of each ligands list
    job_input = [[line] for line in usable_list]
    mol_usable_list = mp.multi_threading(job_input, -1,  add_mol_to_list)
    mol_usable_list = [x for x in mol_usable_list if x is not None]
    # This will be used to check rdkit types later
    rdkit_obj_type = type(mol_usable_list[0][-1])

    parpis, parpis_names = get_PARPis()
    counter = 0
    for inhib,inhib_name in zip(parpis,parpis_names):
        counter = counter-1
        job_input = [[[inhib, mol_list[counter]], mol_list] for mol_list in mol_usable_list]

        mol_usable_list = mp.multi_threading(job_input, -1,  mcs_mols)
        print(inhib_name)
        print(-7 - counter)
        print("")


    save_list = []
    #remove rdkit mol object from list
    for sublist in mol_usable_list:
        temp_list = [x for x in sublist if type(x)!=rdkit_obj_type]
        save_list.append(temp_list)


    save_list.sort(key = lambda x: float(x[-1]),reverse = True)


    with open(new_file, "w") as NF:
        for line in save_list:

            output_line = "\t".join(line) + "\n"
            NF.write(output_line)


