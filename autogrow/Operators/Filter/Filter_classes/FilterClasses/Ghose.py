"""Ghose Filter
This runs a Ghose filter for drug-likeliness.
Ghose filter filters molecules by Molecular weight (MW),
the number of atoms, and the logP value.

To pass the filter a molecule must be
    MW between 160 and 500 dalton
    Number of Atoms: between 20 and 70
    logP  between -0,4 and +5,6

If you use the Ghose Filter please cite:
A.K. Ghose et al.
A knowledge-based approach in designing combinatorial or 
medicinal chemistry libraries for drug discovery. 1. A 
qualitative and quantitative characterization of known drug databases
Journal of Combinatorial Chemistry, 1 (1999), pp. 55-68
"""
import __future__

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.Lipinski as Lipinski
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem.Descriptors as Descriptors
#Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog('rdApp.*')

from autogrow.Operators.Filter.Filter_classes.ParentFilterClass import ParentFilter


class Ghose(ParentFilter):
    """
    This runs a Ghose filter for drug-likeliness.
    Ghose filter filters molecules by Molecular weight (MW),
    the number of atoms, and the logP value.
    
    To pass the filter a molecule must be
        MW between 160 and 500 dalton
        Number of Atoms: between 20 and 70
        logP  between -0,4 and +5,6

    If you use the Ghose Filter please cite:
    A.K. Ghose et al.
    A knowledge-based approach in designing combinatorial or 
    medicinal chemistry libraries for drug discovery. 1. A 
    qualitative and quantitative characterization of known drug databases
    Journal of Combinatorial Chemistry, 1 (1999), pp. 55-68

    Inputs:
    :param class ParentFilter: a parent class to initialize off
    """
    def run_filter(self, mol):
        """
        This runs a Ghose filter for drug-likeliness.
        Ghose filter filters molecules by Molecular weight (MW),
        the number of atoms, and the logP value.
        
        To pass the filter a molecule must be
            MW between 160 and 500 dalton
            Number of Atoms: between 20 and 70
            logP  between -0,4 and +5,6

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be tested if it passes the filters       
        Returns:
        :returns: bool bool: True if the mol passes the filter; False if it fails the filter
        """
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
        else:
            return True
    #
#
