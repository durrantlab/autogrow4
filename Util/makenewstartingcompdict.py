import glob


import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.Lipinski as Lipinski
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.MolSurf as MolSurf



def run_filter(mol):
    """
    ghose
    mw 160-500
    # atoms=20-70
    logP -0,4 to +5,6
    MolMR 40 to 130
    """
    sanitize_string =  Chem.SanitizeMol(mol, sanitizeOps = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors = True)
    if sanitize_string.name != "SANITIZE_NONE":
        return False

    ExactMWt = Descriptors.ExactMolWt(mol)
    if 160 > ExactMWt > 500:
        return False
    
    num_atoms = mol.GetNumAtoms()
    if 20 > num_atoms > 70:
        return False

    num_H_bond_donors = Lipinski.NumHDonors(mol)
    if num_H_bond_donors > 5:
        return False
    
    num_H_bond_acceptors = Lipinski.NumHAcceptors(mol)
    if num_H_bond_acceptors > 10:
        return False
    
    # molar Refractivity
    MolMR = Crippen.MolMR(mol)
    if 40 > MolMR > 130:
        return False

    # molar LogP
    mol_logP = Crippen.MolLogP(mol)

    if -0.4 > mol_logP > 5:
        return False

    ExactMWt = Descriptors.ExactMolWt(mol)
    if ExactMWt > 500:
        return False
    
    num_H_bond_donors = Lipinski.NumHDonors(mol)
    if num_H_bond_donors > 5:
        return False
    
    num_H_bond_acceptors = Lipinski.NumHAcceptors(mol)
    if num_H_bond_acceptors > 10:
        return False
    
    mol_logP = Crippen.MolLogP(mol)
    if mol_logP > 5:
        return False
    
    ExactMWt = Descriptors.ExactMolWt(mol)
    if ExactMWt > 450:
        return False
    PSA = MolSurf.TPSA(mol)
    if PSA > 90:
        return False

    halogen = Chem.MolFromSmarts("[*;#9,#17,#35,#53,#85]")
    number_of_halogens = len(mol.GetSubstructMatches(halogen, maxMatches=8))
    if number_of_halogens > 8:
        return False

    oxygen_smarts = "[#8]"
    number_of_oxygens = len(mol.GetSubstructMatches(Chem.MolFromSmarts(oxygen_smarts),maxMatches=2))
    if number_of_oxygens < 1:
        return False
    
    num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    if num_rotatable_bonds > 15:
        return False
    else:
        return True
## 

if __name__ == "__main__":
    #################################### 
    folder_name = "/home/jspiegel/DataB/jspiegel/projects/autogrow/autogrow/Operators/Mutation/SmileClickChem/Complimentary_mol_dict/"


    files = glob.glob(folder_name + "*.smi")

    print("HI")
    print(files)
    list_of_mols = []
    for file in files:




        temp_list = []
        
        
        with open(file) as smiles_file:
            for line in smiles_file:
                line = line.replace("\n","")
                parts = line.split('\t')      # split line into parts seperated by 4-spaces



                temp_list = [parts[0], parts[1]]


                if temp_list not in list_of_mols:

                    mol = Chem.MolFromSmiles(temp_list[0], sanitize = False)
                    mol = MOH.handleHs(mol)
                    if mol is not None:
                        # filter mols by drug-likeliness
                        filter_result = run_filter(mol)
                        if filter_result is True:
                            list_of_mols.append(temp_list)

                else:
                    print("redundancy")


    output_file_name = "/home/jspiegel/DataB/jspiegel/projects/autogrow/cdc73/small_mol_library.smi"

    # write as a tab delineated .smi file
    with open(output_file_name, 'w') as f:
        for pair in list_of_mols:
            smile_string = pair[0]
            smile_id = pair[1]
            # smile_score_float = smile[2]
            x = str(smile_string + '\t' + smile_id + '\n')
            f.write(x)
            
