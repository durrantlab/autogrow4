import __future__ 

import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS
#Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog('rdApp.*')

import autogrow.Operators.Crossover.SmilesMerge.MergeFunctions.Merge_w_core as MWC
import autogrow.Operators.Crossover.SmilesMerge.MergeFunctions.Dict_and_R_Groups as DnR
import autogrow.Operators.Crossover.SmilesMerge.MergeFunctions.Alignment_and_Breaks as AnB
import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


def run_main_SmilesMerge(vars, Lig_string_1, Lig_string_2):
    """
    This runs the main script for SmileMerge.

    Input:
    :param dict vars: User variables which will govern how the programs runs
    
    :param str Lig_string_1: smile string for lig 1
    :param str Lig_string_2: smile string for lig 2
        example Lig_string_1 = "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"
                Lig_string_2 = "C# CCOc1ccc2ccccc2c1CO"
    Returns:
    :returns: str Ligand_new_smiles: smile string for the child ligand
                              derived from lig_1 and lig_2   
                        Returns None if it failed at any point
    """
    # Lig_string_1 = "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"
    # Lig_string_2 = "C# CCOc1ccc2ccccc2c1CO"
    # Lig_string_1 = "C1 = CC = CC = C1"

    Lig_smile_1 = Chem.MolFromSmiles(Lig_string_1, sanitize = False)
    Lig_smile_2 = Chem.MolFromSmiles(Lig_string_2, sanitize = False)

    # Sanitize
    Lig_smile_1 = Chem.MolFromSmiles(Lig_string_1, sanitize = False)
    Lig_smile_2 = Chem.MolFromSmiles(Lig_string_2, sanitize = False)

    # Sanitize, deprotanate, and reprotanate both molecules    
    mol_1 = MOH.check_sanitization(Lig_smile_1)
    if mol_1 is None:
        return False
    mol_2 = MOH.check_sanitization(Lig_smile_2)
    if mol_2 is None:
        return False
    
    protanate_step = vars["protanate_step"]
    mol_1 = MOH.handleHs(Lig_smile_1, protanate_step)
    mol_2 = MOH.handleHs(Lig_smile_2, protanate_step)

    # check that handleHs() worked for both molecules
    # if fail move on to next pair of molecules
    if mol_1 is None or mol_2 is None:
        return None 

    # make a list of the two rdkit.Chem.rdchem.Mol objects
    mols = [mol_1,mol_2]
    
    # Use the below MCS_H function for Most Common Substructure searching. This will prevent broken rings.
    MCS_results = rdFMCS.FindMCS(mols, matchValences = False, ringMatchesRingOnly = True, completeRingsOnly = False, timeout = vars["max_time_MCS_thorough"])

    if MCS_results.canceled == True:
        return None

    # confirm that this meets the minimum number of matching atoms
    if MCS_results.numAtoms < vars["min_atom_match_MCS"]:
        return None

    ### Convert MCS_res from into usable and referable forms
    MCS_Mol = Chem.MolFromSmarts(MCS_results.smartsString)

    # handle_MCS_align_labeling_and_cyclicbreaks
    mol_1, mol_2, MCS_Mol = AnB.handle_MCS_align_labeling_and_cyclicbreaks(mol_1, mol_2, MCS_Mol)

    # confirm that this meets the minimum number of matching atoms
    if mol_1 is None or mol_2 is None or MCS_Mol is None:
        return None
    # This is a very big step. This is where all of the atoms are fragmented
    # to determine what are the R-groups options. R-groups are consolidated into
    # B-groups for each molecule individually, the 2 dictionaries of B-groups
    # are consolidated into a master B-group which is fed into the Mapping class
    # which will randomly select a set of non-clashing B-groups, such that no 2 chosen
    # B's will have connections to the same anchors (anchors are atoms in the common core).
    # It returns a list of smilestrings for the chosen R-groups 
    # It returns None if a step failed
    
    Rs_chosen_smiles = DnR.handle_dicts_and_select_B_groups(mol_1, mol_2, MCS_Mol)
    if Rs_chosen_smiles is None:
        return None

    # MERGE ALL CHOSEN R-GROUPS WITH THE CORE
    Ligand_new_mol = MWC.Merge_smiles_with_core(Rs_chosen_smiles, MCS_Mol)
    if Ligand_new_mol is None:
        return None

    Ligand_new_mol = MOH.check_sanitization(Ligand_new_mol)
    if Ligand_new_mol is None:
        return None

    # REMOVE ALL THE ISOTOPES IN THE NEW MOLECULE
    Ligand_new_mol_final = MWC.remove_all_isolabels(Ligand_new_mol)

    # Remove any fragments incase 1 made it through
    Ligand_new_mol_final = MOH.handle_frag_check(Ligand_new_mol_final)
    if Ligand_new_mol_final is None:
        return None
    
    # Make sure there are no unassigned atoms which made it through 
    # These are very unlikely but possible
    Ligand_new_mol_final = MOH.check_for_unassigned_atom(Ligand_new_mol_final)
    if Ligand_new_mol_final is None:
        return None

    Ligand_new_mol = MOH.check_sanitization(Ligand_new_mol_final)
    if Ligand_new_mol is None:
        return None

    Ligand_new_smiles = Chem.MolToSmiles(Ligand_new_mol, isomericSmiles = True)

    return Ligand_new_smiles

#
