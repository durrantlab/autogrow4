"""Lipinski Strict
This runs a Strict Lipinski filter. Lipinski filter refines for orally
available drugs. It filters molecules by Molecular weight (MW), the number of
hydrogen donors, the number hydrogen acceptors, and the logP value.

To pass the Lipinski filter a molecule must be:
    MW: Max 500 dalton
    Number of H acceptors: Max 10
    Number of H donors: Max 5
    logP Max +5.0

If you use the Lipinski Filter please cite: C.A. Lipinski et al. Experimental
and computational approaches to estimate solubility and permeability in drug
discovery and development settings Advanced Drug Delivery Reviews, 46 (2001),
pp. 3-26
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
import rdkit.Chem.Lipinski as Lipinski  # type: ignore
import rdkit.Chem.Crippen as Crippen  # type: ignore
import rdkit.Chem.Descriptors as Descriptors  # type: ignore
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class LipinskiStrictFilter(SmilesFilterBase):
    """
    This runs a Strict Lipinski filter. Lipinski filter refines for orally
    available drugs. It filters molecules by Molecular weight (MW), the number
    of hydrogen donors, the number hydrogen acceptors, and the logP value.

    This is a strict Lipinski which means a ligand must pass all the
    requirements.

    If you use the Lipinski Filter please cite: C.A. Lipinski et al.
    Experimental and computational approaches to estimate solubility and
    permeability in drug discovery and development settings Advanced Drug
    Delivery Reviews, 46 (2001), pp. 3-26

    Inputs:
    :param class ParentFilter: a parent class to initialize off
    """

    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        This runs a Strict Lipinski filter. Lipinski filter refines for orally
        available drugs. It filters molecules by Molecular weight (MW), the
        number of hydrogen donors, the number hydrogen acceptors, and the logP
        value.

        This is a strict Lipinski which means a ligand must pass all the
        requirements.

        To pass the Lipinski filter a molecule must be:
            MW: Max 500 dalton
            Number of H acceptors: Max 10
            Number of H donors: Max 5
            logP Max +5.0

        If you use the Lipinski Filter please cite: C.A. Lipinski et al.
        Experimental and computational approaches to estimate solubility and
        permeability in drug discovery and development settings Advanced Drug
        Delivery Reviews, 46 (2001), pp. 3-26

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be
            tested if it passes the filters

        Returns:
        :returns: bool bool: True if the mol passes the filter; False if it
          fails the filter
        """

        exact_mwt = Descriptors.ExactMolWt(mol)
        if exact_mwt > 500:
            return False

        num_hydrogen_bond_donors = Lipinski.NumHDonors(mol)
        if num_hydrogen_bond_donors > 5:
            return False

        num_hydrogen_bond_acceptors = Lipinski.NumHAcceptors(mol)
        if num_hydrogen_bond_acceptors > 10:
            return False

        mol_log_p = Crippen.MolLogP(mol)
        if mol_log_p > 5:
            return False

        # Passed all filters
        return True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
        return (
            "SMILES Filters",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Lipinski filters for orally available drugs following Lipinski rule of fives. \
        Filters by molecular weight, logP and number of hydrogen bond donors and acceptors. \
        Strict implementation means a ligand must pass all requirements.",
                )
            ],
        )