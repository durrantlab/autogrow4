import itertools
import copy
import random
import os
import sys

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem import rdFMCS
from rdkit import RDLogger
# Turn off warnings
RDLogger.DisableLog('rdApp.*')

import support_scripts.Multiprocess as mp
import support_scripts.mol_object_handling as MOH

def get_atom_w_iso_num(mol, iso_num):
    for atom in mol.GetAtoms():
        if atom.GetIsotope() == iso_num:
            return atom.GetIdx()
        else:
            continue
    return None
#
def label_iso_num_w_idx(mol):
    for atom in mol.GetAtoms():
        atom.SetIsotope(atom.GetIdx())
    return mol
#
def get_permutations_of_rotatable_bonds_to_cut(mol, c_c_bonds_off=False):
    RotatableBond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    Rotatable_Bonds_set = mol.GetSubstructMatches(RotatableBond)

    
    Rotatable_Bonds_to_frag = []
    for rot_bond in Rotatable_Bonds_set:
        atom1 = mol.GetAtomWithIdx(rot_bond[0])
        atom2 = mol.GetAtomWithIdx(rot_bond[1])
        atom_isos = [atom1.GetIsotope(), atom2.GetIsotope()]
        bond = mol.GetBondBetweenAtoms(rot_bond[0], rot_bond[1])
        if bond.GetIsAromatic() == True:
            continue
        # Remove any bonds including Hydrogen
        elif atom1.GetAtomicNum()==1 or atom2.GetAtomicNum()==1:
            continue
        # Remove any C-C single bonds
        elif atom1.GetAtomicNum()==6 and atom2.GetAtomicNum()==6:
            if c_c_bonds_off==True:continue
            else: 
                Rotatable_Bonds_to_frag.append(atom_isos)

        else:
            Rotatable_Bonds_to_frag.append(atom_isos)
    permutations_of_bonds_to_remove = []

    for i in range(1, len(Rotatable_Bonds_to_frag) + 1):  #  xrange will return the values 1,2,3,4 in this loop
        temp_perm_list = list(itertools.combinations(Rotatable_Bonds_to_frag, i))
        permutations_of_bonds_to_remove.extend(temp_perm_list)
    return permutations_of_bonds_to_remove
#
def remove_atoms(mol, list_of_idx_to_remove):
    """
    This function removes atoms from an rdkit mol based on
    a provided list. The RemoveAtom function in Rdkit requires
    converting the mol to an more editable version of the rdkit mol
    object (Chem.EditableMol).

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol
    :param list list_of_idx_to_remove: a list of idx values to remove
                                        from mol
    Returns:
    :returns: rdkit.Chem.rdchem.Mol new_mol: the rdkit mol as input but with
                                            the atoms from the list removed
    """

    if mol is None:
        return None

    try:
        atomsToRemove = list_of_idx_to_remove
        atomsToRemove.sort(reverse = True)
    except:
        return None

    try:
        em1 = Chem.EditableMol(mol)
        for atom in atomsToRemove:
            em1.RemoveAtom(atom)

        new_mol = em1.GetMol()

        return new_mol
    except:
        return None
#
def get_brics_permutations(mol,min_frag_size = 3):
    from rdkit.Chem.BRICS import BRICSDecompose
    res = list(BRICSDecompose(mol,returnMols=True,minFragmentSize=min_frag_size))
    smis = [Chem.MolToSmiles(x,True) for x in res]

    # Get larger pieces
    res = list(BRICSDecompose(mol,returnMols=True,keepNonLeafNodes=True,minFragmentSize=min_frag_size))
    smis.extend([Chem.MolToSmiles(x,True) for x in res])
    clean_frag_list = []
    for x in res:
        list_to_remove = []
        for i in x.GetAtoms():
            if i.GetAtomicNum()==0:
                list_to_remove.append(i.GetIdx())
        
        x = remove_atoms(x, list_to_remove)

        for atom in x.GetAtoms():
            atom.SetIsotope(0)

        clean_frag_list.append(Chem.MolToSmiles(x))
    list(set(list(clean_frag_list)))

    return clean_frag_list

def remove_bonds(mol, list_of_atomiso_bondsets_to_remove):
    """
    This function removes bond from an rdkit mol based on
    a provided list. This list is a list of sets, with each set containing
    two atoms with the isotope label of that atom. Using Isotopes is to ensure
    that atom Idx dont change.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol
    :param list list_of_idx_to_remove: a list of idx values to remove
                                        from mol
    Returns:
    :returns: rdkit.Chem.rdchem.Mol new_mol: the rdkit mol as input but with
                                            the atoms from the list removed
    """
    # None's often end up in a pipeline use of RDKit so we handle this data type as return None
    # instead of raise TypeError
    if mol is None:
        return None
    
    # If mol is wrong data type (excluding None) raise TypeError
    if type(mol) != rdkit.Chem.rdchem.Mol and type(mol) != rdkit.Chem.rdchem.RWMol:
        printout = "mol is the wrong data type. \n"
        printout = printout + "Input should be a rdkit.Chem.rdchem.Mol\n"
        printout = printout + "Input mol was {} type.".format(type(mol))
        raise TypeError(printout)
    new_mol = copy.deepcopy(mol)
    if len(list_of_atomiso_bondsets_to_remove) == 0:
        return None
    for atomiso_bondsets in list_of_atomiso_bondsets_to_remove:
        if len(atomiso_bondsets) == 0:
            continue
        elif len(atomiso_bondsets) != 2:
            printout = "list_of_atomiso_bondsets_to_remove needs to be 2 isolabels for the atoms"
            raise TypeError(printout)
            
        atom_1_idx = int(get_atom_w_iso_num(new_mol,atomiso_bondsets[0]))
        atom_2_idx = int(get_atom_w_iso_num(new_mol,atomiso_bondsets[1]))

        try:
            new_mol = Chem.FragmentOnBonds(new_mol, [atom_1_idx, atom_2_idx],addDummies=False)
        except:
            return None
        
        new_mol = MOH.check_sanitization(new_mol)
        if new_mol == None:
            return None
    new_mol = MOH.check_sanitization(new_mol)
    if new_mol== None:
        return None
    else:
        return new_mol
#
def make_list_of_all_unique_frags(fragment_list):
    """
    This function takes a list of all molecules after fragmentation and seperates the
    the fragments into individual rdkit mol objects, sanitizes each, removes isotopes
    and converts them into a SMILES string. The SMILES are compiled into a list,
    and then redudant strings are reduced to a single entry.
    
    It returns a list of all unique sanitized canonical SMILES for every fragment made 
    from all permutations of bond breaking. 
    
    Inputs:
    :param list fragment_list: list of fragmented rdkit mols which haven't been seperated
        yet
        
    Returns:
    :returns: list clean_frag_list: List of unique sanitized SMILES strings from all objects
                in fragment_list. Isotope labels are also removed here.
    """
    clean_frag_list = []
    for fragments in fragment_list:
        frags = Chem.GetMolFrags(fragments, asMols = True, sanitizeFrags = False)
        for frag in frags:
            frag = MOH.check_sanitization(frag)
            if frag == None:
                continue

            # Remove those under 2 atoms minumum
            list_mol_atoms = frag.GetAtoms()
            if len(list_mol_atoms) < 3:
                continue

            for atom in frag.GetAtoms():
                atom.SetIsotope(0)
            clean_frag_list.append(Chem.MolToSmiles(frag,isomericSmiles=True,canonical=True))
        list(set(list(clean_frag_list)))

    return clean_frag_list
#
def make_unique_lig_id(parent_lig_name, current_lig_list):
    """
    This will make a ligand name from the parent name. Keep start names simple.
    
    Format of names:
        - str(parent_lig_name) + "_Frag_" + str(random_int)
        
    Inputs:
    :param str parent_lig_name: str of the ligand Id for the parent mol
    :param list current_lig_list: the list of names already taken
        
    Returns:
    :returns: str unique_lig_id: A unique ID/name for the child ligand.
    """
    if type(parent_lig_name) != str:
        raise Exception("Ligand ID's to seed this must have Unique string IDs")
    parent_lig_name = parent_lig_name.replace(" ","")
    picked_name = False
    while picked_name == False:
        random_int = random.choice(range(100000,999999))
        unique_lig_id = str(parent_lig_name) + "_Frag_" + str(random_int)
        if unique_lig_id in current_lig_list:
            continue
        else:
            picked_name = True
            break
    return unique_lig_id
#      
def make_frag_list_for_one_mol(mol_info, frags_per_seed_lig, run_brics, run_frag, c_c_bonds_off=False):
    """
    This will take a ligand string and ID encased in the list mol_info.
    This will then be fragmented along all non Carbon-carbon rotatable bonds which
    are not aromatic.
    
    It will make all permutations of all potential bond breaks, reduce to only unique
    fragments and than pick the number of chosen fragments. Then it will create unique ID's
    for each and return a list of lists containing the chosen unique fragments.
    
    Inputs:
    :param list mol_info: list containing [mol_string, mol_id] 
                mol_info[0] = the SMILE string of the parent mol
                mol_info[1] = the Unique ID of the parent mol
    :param int frags_per_seed_lig: the number of fragments to generate from a given parent mol
        
    Returns:
    :returns: list final_frag_list: A list of lists containing the chosen unique fragments.
            final_frag_list[0] = [SMILE, mol_id]
    """
    mol_Smile = mol_info[0]
    lig_id = mol_info[1]

    mol = Chem.MolFromSmiles(mol_Smile, sanitize=False)
    mol = MOH.check_sanitization(mol)
    if mol == None:
        printout = "\nMolecule {} failed to sanitize. Could not make any fragments from it".format(mol_id)
        raise Exception(printout)
    mol_Smile = Chem.MolToSmiles(mol,isomericSmiles=True,canonical=True)
    
    
    mol = label_iso_num_w_idx(mol)
    mol_copy = copy.deepcopy(mol)
    permutations_of_bonds_to_remove = get_permutations_of_rotatable_bonds_to_cut(mol_copy,c_c_bonds_off)

    fragment_list = []
    for bond_set_to_del in permutations_of_bonds_to_remove:
        mol_copy = copy.deepcopy(mol)
        x = remove_bonds(mol_copy, bond_set_to_del)
        if x == None:
            continue
        fragment_list.append(x)
    
    clean_frag_list=[]
    if run_frag == True:
        clean_frag_list = make_list_of_all_unique_frags(fragment_list)
        clean_frag_list = list(set(clean_frag_list))
    
    if run_brics==True:
        mol_copy = copy.deepcopy(mol)
        bric_mols = get_brics_permutations(mol_copy, min_frag_size = 3)
        
        clean_frag_list.extend(bric_mols)
        clean_frag_list = list(set(clean_frag_list))
    
    if len(clean_frag_list) == 0:
        printout =  "\nNo fragments were made for {}.\n".format(lig_id)
        print(printout)
        return [[mol_Smile, lig_id]]

    # Pick the number of ligands to make
    final_frag_list = [[mol_Smile, lig_id]]

    if frags_per_seed_lig == -1:
        printout = "\nFor {}: {} fragmented were made.".format(lig_id,len(clean_frag_list))
        print(printout)
        for frag in clean_frag_list:
            unique_lig_id = make_unique_lig_id(lig_id, final_frag_list)
            temp_frag_info = [frag, unique_lig_id]
            final_frag_list.append(temp_frag_info)
    return final_frag_list
    
    if frags_per_seed_lig > len(clean_frag_list):
        printout = "\nFor {}: Not enough fragments made to make {} fragmented ligands.".format(lig_id,frags_per_seed_lig)
        printout =printout +"\n \t Only made {}. Will take all of these instead + 1 for the original ligand.".format(len(clean_frag_list))
        print(printout)

        for frag in clean_frag_list:
            unique_lig_id = make_unique_lig_id(lig_id, final_frag_list)
            temp_frag_info = [frag, unique_lig_id]
            final_frag_list.append(temp_frag_info)

    elif frags_per_seed_lig < len(clean_frag_list):
        printout = "\nFor {}: Made more fragments than requested. Made {} fragments.\n".format(lig_id,len(clean_frag_list))
        printout = printout + "\tWill Chose {} random unique fragments plus the original seed molecule. \n"
        print(printout)
        random_list_frags = random.sample(clean_frag_list, frags_per_seed_lig)
        for frag in random_list_frags:
            unique_lig_id = make_unique_lig_id(lig_id, final_frag_list)
            temp_frag_info = [frag, unique_lig_id]
            final_frag_list.append(temp_frag_info)
        
    else:
        printout = "\nFor {}: Made exact {} unique fragments..".format(frags_per_seed_lig)
        printout =printout +"\n \t Will take all of these instead + 1 for the original ligand."
        print(printout)

        for frag in clean_frag_list:
            unique_lig_id = make_unique_lig_id(lig_id, final_frag_list)
            temp_frag_info = [frag, unique_lig_id]
            final_frag_list.append(temp_frag_info)
            
    return final_frag_list
    



if __name__ == "__main__":
    try:
        smi_file_name = str(sys.argv[1])
        output_smi_file_name = str(sys.argv[2])
        frags_per_seed_lig = int(sys.argv[3])
        number_of_processors = int(sys.argv[4])
        if sys.argv[5] =="True" or sys.argv[5]=="true" or sys.argv[5]=="1":
            run_brics = True
        else:
            run_brics = False
        if sys.argv[6] =="True" or sys.argv[6]=="true" or sys.argv[6]=="1":
            run_frag = True
        else:
            run_frag = False

        try:
            c_c_bonds_off = sys.argv[7]
            if c_c_bonds_off==True or c_c_bonds_off=="True" or c_c_bonds_off=="true" or c_c_bonds_off==1 or c_c_bonds_off=="1":
                c_c_bonds_off=True
            else:
                c_c_bonds_off=False
            
        except:
            c_c_bonds_off = False

    except:
        printout = "to run the fragmentor: \n"

        printout = printout + "1st argv is the path to the.smi file of ligands to fragment\n"
        printout = printout + "2nd argv is the path to output .smi file of generated fragments\n"
        printout = printout + "3rd argv is number of fragments to generate per parent ligand. set to -1 if you want all output fragments.\n"
        printout = printout + "4th argv is number of processors\n"
        printout = printout + "5th argv is to add brics decomposition True or False \n"
        printout = printout + "6th argv is to add standard fragment along all rotatable bonds True or False \n"
        printout = printout + "7th argv is an optional arguement; if True it will ignore fragments on C-C bonds; if False it will fragment along C-C bonds; default is False\nn"
        print(printout+"\n\n\n")
        sys.exit(0)

    print("")
    print("STARTING FRAGMENTER")
    print("frags_per_seed_lig: ", frags_per_seed_lig)
    print("smi_file_name: ", smi_file_name)
    print("########")
    if os.path.isfile(smi_file_name) == False:
        raise Exception("\n.SMI file not found.\n") 
    print("Importing .smi file")

    list_of_ligands = []
    with open(smi_file_name, 'r') as smiles_file:
        line_counter = 0
        for line in smiles_file:
            line_counter = line_counter + 1
            line = line.replace("\n","")
            parts = line.split('\t')      # split line into parts seperated by 4-spaces
            if len(parts) == 1:
                parts = line.split('    ')      # split line into parts seperated by 4-spaces

            if len(parts) == 2 or len(parts) > 2:
                mol_string = parts[0]
                mol_id = parts[1]
                if type(mol_id) != str:
                    print("Miss Formatted within .SMI. Line number {}".format(str(line_counter)))
                    continue

                try:
                    mol = Chem.MolFromSmiles(mol_string, sanitize=False)
                except:
                    print("Miss Formatted within .SMI. Line number {}".format(str(line_counter)))
                    continue
                mol = MOH.check_sanitization(mol)
                if mol == None:
                    continue
                
                mol_Smile = Chem.MolToSmiles(mol,isomericSmiles=True,canonical=True)
                mol_info = [mol_Smile, mol_id]
                list_of_ligands.append(mol_info)

            else:
                continue
    print("Was able to import and sanitize {} ligands from the .smi.".format(len(list_of_ligands)))
    if line_counter != len(list_of_ligands):
        print("\t Failed to sanitize/import {} ligands from the .smi".format(line_counter - len(list_of_ligands)))
    original_ligands = copy.deepcopy(list_of_ligands)
    print("########")


    # create a set of jobs to multithread the fragmentation
    job_input = [tuple([mol_info, frags_per_seed_lig, run_brics, run_frag, c_c_bonds_off]) for mol_info in list_of_ligands]
    list_of_ligands = None
    output = mp.multi_threading(job_input, number_of_processors,  make_frag_list_for_one_mol)
    

    print("Finish multithread\n")
    #

    output = [x for x in output if x is not None]
    output = [x for x in output if x is not ""]
    
    initial_output_reduce = []

    for x in output:
        initial_output_reduce.extend(x)
    output = None
    
    initial_output_reduce = [x for x in initial_output_reduce if x[0] is not ""]
    initial_output_reduce = [x for x in initial_output_reduce if x[1] is not ""]
    
    # Reduce smile redundancies:
    smiles_list = []
    output_reduce = []
    for x in initial_output_reduce:
        if x[0] in smiles_list:continue
        else:
            output_reduce.append(x)
            smiles_list.append(x[0])
            
    final_mol_list = []
    master_smile_list = []
    master_id_list = []
    for x in output_reduce:
        temp_smile = x[0]
        temp_id =  x[1]
        if temp_smile in master_smile_list:
            continue
        if temp_id in master_id_list:
            continue

        # Append to master lists and final_mol_list 
        final_mol_list.append(x)
        master_smile_list.append(temp_smile)
        master_id_list.append(temp_id)


    # convert list of mols to a print statement
    printout = ""
    for x in final_mol_list:
        printout = printout + x[0] + "\t" + x[1] + "\n"
    
    print("####")
    print("\nSaving list to file")
    with open(output_smi_file_name, 'w') as f:
        f.write(printout)

    
    print("Number of parent ligands:         {}".format(len(job_input)))
    print("Number of new fragmented ligands: {}".format(len(final_mol_list)-len(job_input)))
    print("Total number ligs in output file: {}".format(len(final_mol_list)))
    print("Finished")



# Example submission
# python fragmenter_of_smi_mol.py ../source_compounds/PARPi.smi ../source_compounds/new_bric.smi -1 -1 True True 

# Frag and BRIC W C-C
# python fragmenter_of_smi_mol.py ../source_compounds/PARPi.smi ../source_compounds/new_bric.smi -1 -1 True True 

# Frag and BRIC only
# python fragmenter_of_smi_mol.py ../source_compounds/PARPi.smi ../source_compounds/only_bric.smi -1 -1 True False False

# Frag only no C-C
# python fragmenter_of_smi_mol.py ../source_compounds/PARPi.smi ../source_compounds/new_bric.smi -1 -1 False True False

# Frag only with C-C
# python fragmenter_of_smi_mol.py ../source_compounds/PARPi.smi ../source_compounds/new_bric.smi -1 -1 False True True


