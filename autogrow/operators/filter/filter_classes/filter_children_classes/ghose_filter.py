"""Ghose Filter
This runs a Ghose filter for drug-likeliness. Ghose filter filters molecules
by Molecular weight (MW), the number of atoms, and the logP value.

We protonate the mol in this filter because hydrogens affect
atom count. Our Ghose implementation counts hydrogens in against
the total number of atoms.

To pass the filter a molecule must be:
    MW between 160 and 480 dalton
    Number of Atoms: between 20 and 70
    logP  between -0,4 and +5,6

If you use the Ghose Filter please cite: A.K. Ghose et al. A knowledge-based
approach in designing combinatorial or medicinal chemistry libraries for drug
discovery. 1. A qualitative and quantitative characterization of known drug
databases Journal of Combinatorial Chemistry, 1 (1999), pp. 55-68
"""

import __future__

import copy

import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
import rdkit.Chem.Lipinski as Lipinski  # type: ignore
import rdkit.Chem.Crippen as Crippen  # type: ignore
import rdkit.Chem.Descriptors as Descriptors  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter


class GhoseFilter(ParentFilter):
    """
    This runs a Ghose filter for drug-likeliness. Ghose filter filters
    molecules by Molecular weight (MW), the number of atoms, and the logP
    value.

    We protonate the mol in this filter because hydrogens affect
    atom count. Our Ghose implementation counts hydrogens in against
    the total number of atoms.

    To pass the filter a molecule must be:
        MW between 160 and 480 dalton
        Number of Atoms: between 20 and 70
        logP  between -0,4 and +5,6

    If you use the Ghose Filter please cite: A.K. Ghose et al. A
    knowledge-based approach in designing combinatorial or medicinal chemistry
    libraries for drug discovery. 1. A qualitative and quantitative
    characterization of known drug databases Journal of Combinatorial
    Chemistry, 1 (1999), pp. 55-68

    Inputs:
    :param class ParentFilter: a parent class to initialize off
    """

    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        This runs a Ghose filter for drug-likeliness. Ghose filter filters
        molecules by Molecular weight (MW), the number of atoms, and the logP
        value.

        We protonate the mol in this filter because hydrogens affect
        atom count. Our Ghose implementation counts hydrogens in against
        the total number of atoms.

        To pass the filter a molecule must be:
            MW between 160 and 480 dalton
            Number of Atoms: between 20 and 70
            logP  between -0,4 and +5,6

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be
            tested if it passes the filters

        Returns:
        :returns: bool bool: True if the mol passes the filter; False if it
            fails the filter
        """
        # Make a copy of the mol so we can AddHs without affecting other filters
        # number of atoms is altered by the presence/absence of hydrogens.
        # Our Ghose filter counts hydrogenss towards atom count
        copy_mol = copy.deepcopy(mol)
        copy_mol = Chem.AddHs(copy_mol)
        exact_mwt = Descriptors.ExactMolWt(copy_mol)
        if (exact_mwt < 160) or (exact_mwt > 480):
            return False

        num_atoms = copy_mol.GetNumAtoms()
        if (num_atoms < 20) or (num_atoms > 70):
            return False

        # molar Refractivity
        MolMR = Crippen.MolMR(copy_mol)
        if (MolMR < 40) or (MolMR > 130):
            return False

        # molar LogP
        mol_log_p = Crippen.MolLogP(copy_mol)
        if (mol_log_p < -0.4) or (mol_log_p > 5.6):
            return False

        # passed all filters
        return True
