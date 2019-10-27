# convert pdbs into smiles

# This script will take a folder and convert all pdb files into a single texted file.
# The text file will contain smiles strings of the respective pdb and the name of the file.

# This script won't work on window as directory format is different than Linux.

# Run example:
# 
# output example:
# python convert_pdb_to_smiles.py /home/jspiegel/DataB/jspiegel/projects/autogrow/3.1.2/autogrow/fragments/MW_150 /home/jspiegel/DataB/jspiegel/projects/autogrow/smileConversionTools/output_smiles
# CC1COC(=O)OC(=O)O1    ZINC60039447
# O = C1OC(=O)N2CCC12    ZINC59199492
# O = C1CC(C(=O)O)C(=O)O1    ZINC59901386
import __future__

import glob
import os
import sys


from rdkit import Chem
       
import support_scripts.MolObjectHandling as MOH

sourceFolder = sys.argv[1]
outputFolder = sys.argv[2]


def make_smile_list(sub_folder):
    sub_folder = sub_folder+"/"
    smilesList = []
    for pdb in glob.glob("/" + sub_folder+"*.pdb"):
        print("PDB")
        print(pdb)
        try:
            mol = Chem.MolFromPDBFile(pdb)

            mol_sanitized = MOH.check_sanitization(mol)
            if mol_sanitized is not None:
                smiles = Chem.MolToSmiles(mol_sanitized, isomericSmiles=True)
                fileName = os.path.basename(pdb)
                fileStripped = fileName.replace(".pdb","")
                output_data = smiles + "\t" + fileStripped

        except:
            pass
        smilesList.append(output_data)
        print(smilesList)
    return smilesList

if __name__ == "__main__":
    master_list = []
    # for folder in glob.glob(sourceFolder + "/*"):
    #     print(folder)
    #     for sub_folder in glob.glob(folder + "/*"):
    #         print(sub_folder)
    #         print("")

    print("")
    for folder in glob.glob(sourceFolder + "/*"):
        print("FOLDER")
        folder_name_short = os.path.basename(folder)
        for sub_folder in glob.glob(folder + "/*"):
            print("SUBFOLDER")
            smilesList = make_smile_list(sub_folder)
            print(sub_folder)
            master_list.extend(smilesList)


            sub_folder_name_short = os.path.basename(sub_folder)
            outputFile = outputFolder+"/" + folder_name_short + "_"+ sub_folder_name_short + ".smi"


            with open(outputFile,"w") as f:
                f.write("\n".join(smilesList))
            
            



        