"""
Strict Lipinski Filter for oral drug availability.

This module implements a Strict Lipinski filter to refine molecules for oral
drug availability. It filters molecules based on Molecular Weight (MW), number
of hydrogen donors and acceptors, and LogP value.

To pass the Lipinski filter, a molecule must meet all of the following criteria:
    - MW: Max 500 dalton
    - Number of H acceptors: Max 10
    - Number of H donors: Max 5
    - LogP: Max +5.0

If you use the Lipinski Filter, please cite: C.A. Lipinski et al. Experimental
and computational approaches to estimate solubility and permeability in drug
discovery and development settings Advanced Drug Delivery Reviews, 46 (2001),
pp. 3-26
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound
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
    Implements a Strict Lipinski filter for oral drug availability.

    This class runs a Strict Lipinski filter to refine molecules for oral drug
    availability. It filters molecules based on Molecular Weight (MW), number of
    hydrogen donors and acceptors, and LogP value.

    This strict implementation requires a ligand to pass all the requirements.

    To pass the Lipinski filter, a molecule must meet all of the following
    criteria:
        - MW: Max 500 dalton
        - Number of H acceptors: Max 10
        - Number of H donors: Max 5
        - LogP: Max +5.0

    If you use the Lipinski Filter, please cite: C.A. Lipinski et al.
    Experimental and computational approaches to estimate solubility and
    permeability in drug discovery and development settings Advanced Drug
    Delivery Reviews, 46 (2001), pp. 3-26
    """

    def run_filter(self, predock_cmpd: Compound) -> bool:
        """
        Run the Strict Lipinski filter on a given molecule.

        This method applies the Strict Lipinski filter criteria to determine if
        a molecule is likely to be orally available as a drug. It checks the
        molecular weight, number of hydrogen bond donors and acceptors, and LogP
        value.

        This strict implementation requires a molecule to pass all the
        requirements.

        Args:
            predock_cmpd (PostDockedCompound): A PostDockedCompound to be tested.

        Returns:
            bool: True if the molecule passes all filter criteria, False
                otherwise.

        If you use the Lipinski Filter, please cite: C.A. Lipinski et al.
        Experimental and computational approaches to estimate solubility and
        permeability in drug discovery and development settings Advanced Drug
        Delivery Reviews, 46 (2001), pp. 3-26
        """
        mol = self.predock_cmpd_to_rdkit_mol(predock_cmpd)
        if mol is None:
            return False

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
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the Strict
        Lipinski Filter. It allows users to enable the filter via command-line
        options.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of ArgumentVars objects defining the
                command-line arguments.
        """
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
