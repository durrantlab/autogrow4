import rdkit
from rdkit import Chem
import copy
import json
import sys
import os
import pickle

import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp

import glob


import support_scripts.MolObjectHandling as MOH
import support_scripts.Multiprocess as mp

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.Lipinski as Lipinski
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.MolSurf as MolSurf
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

#BRENK
# Make a list of the BRENK filter.
params_BRENK = FilterCatalogParams()
params_BRENK.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
# # Make a list of the NIH filter.
params_NIH = FilterCatalogParams()
params_NIH.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
# Make a list of the PAINS filter.
params_PAINS = FilterCatalogParams()
params_PAINS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)


# list_of_params = [params_BRENK,params_NIH,params_PAINS_A,params_PAINS_B,params_PAINS_C]
list_of_params = [params_BRENK, params_NIH, params_PAINS]
# list_of_params = [params_PAINS_A,params_PAINS_B,params_PAINS_C]
stringent_filters_list = [FilterCatalog.FilterCatalog(x) for x in list_of_params]
lient_filters_list = [FilterCatalog.FilterCatalog(params_PAINS)]

list_of_params=None
params_BRENK = None
params_NIH = None
params_PAINS = None
params_PAINS_A = None
params_PAINS_B = None
params_PAINS_C = None



def run_normal_filter(list_of_mol_info):
    """
    ghose
    mw 160-500
    # atoms=20-70
    logP -0,4 to +5,6
    MolMR 40 to 130

    list_of_mol_info = [smiles_str,zinc_id, mol]
    """
    fail_counter = 0
    mol = list_of_mol_info[2]
    sanitize_string =  Chem.SanitizeMol(mol, sanitizeOps = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors = True)
    if sanitize_string.name != "SANITIZE_NONE":
        return None

    ExactMWt = Descriptors.ExactMolWt(mol)
    if ExactMWt > 500:
        fail_counter = fail_counter +1
    
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 70:
        fail_counter = fail_counter +1

    num_H_bond_donors = Lipinski.NumHDonors(mol)
    if num_H_bond_donors > 5:
        fail_counter = fail_counter +1
    
    num_H_bond_acceptors = Lipinski.NumHAcceptors(mol)
    if num_H_bond_acceptors > 10:
        fail_counter = fail_counter +1
    
    # molar Refractivity
    MolMR = Crippen.MolMR(mol)
    if 40 > MolMR > 130:
        fail_counter = fail_counter +1

    # molar LogP
    mol_logP = Crippen.MolLogP(mol)
    if -0.4 > mol_logP > 5:
        fail_counter = fail_counter +1

    
    # EXCLUDING THIS FILTER AS IT IS TOO PENALIZING!!!
    PSA = MolSurf.TPSA(mol)
    if PSA > 90:
        fail_counter = fail_counter +1
        
    halogen = Chem.MolFromSmarts("[*;#9,#17,#35,#53,#85]")
    number_of_halogens = len(mol.GetSubstructMatches(halogen, maxMatches=8))
    if number_of_halogens > 8:
        fail_counter = fail_counter +1

    # EXCLUDING THIS FILTER AS IT IS TOO PENALIZING!!!
    oxygen = Chem.MolFromSmarts("[#8]")
    number_of_oxygens = len(mol.GetSubstructMatches(oxygen, maxMatches=2))
    if number_of_oxygens < 1:
        fail_counter = fail_counter +1

    num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    if num_rotatable_bonds > 15:
        fail_counter = fail_counter +1
    
    if fail_counter > 3:
        return None
    else:

        zinc_id = list_of_mol_info[1]
        return fail_counter, zinc_id
   
def run_filter_on_mol_substructurefilters_lienent(list_of_mol_info):
    """
    ghose
    mw 160-500
    # atoms=20-70
    logP -0,4 to +5,6
    MolMR 40 to 130

    list_of_mol_info = [smiles_str,zinc_id, mol]
    """
    mol = list_of_mol_info[2]

    for i in range(0,len(lient_filters_list)):
        filters = lient_filters_list[i]
        
        if filters.HasMatch(mol) is True:
            return None
    
    else:
    
        zinc_id = list_of_mol_info[1]
        return zinc_id
    

def run_filter_on_mol_substructurefilters_strict(list_of_mol_info):
    """
    ghose
    mw 160-500
    # atoms=20-70
    logP -0,4 to +5,6
    MolMR 40 to 130

    list_of_mol_info = [smiles_str,zinc_id, mol]
    """
    mol = list_of_mol_info[2]

    for i in range(0,len(stringent_filters_list)):
        filters = stringent_filters_list[i]
        
        if filters.HasMatch(mol) is True:
            return None
    
    else:
    
        zinc_id = list_of_mol_info[1]
        return zinc_id
    
