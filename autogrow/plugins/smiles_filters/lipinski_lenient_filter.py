"""
Lenient Lipinski Filter for oral drug availability.

This module implements a Lenient Lipinski filter to refine molecules for oral
drug availability. It filters molecules based on Molecular Weight (MW), number
of hydrogen donors and acceptors, and LogP value.

To pass the Lipinski filter, a molecule must meet the following criteria:
    - MW: Max 500 dalton
    - Number of H acceptors: Max 10
    - Number of H donors: Max 5
    - LogP: Max +5.0

This lenient implementation allows one violation of the Lipinski Rule of 5
constraints.

If you use the Lipinski Filter, please cite: C.A. Lipinski et al. Experimental
and computational approaches to estimate solubility and permeability in drug
discovery and development settings Advanced Drug Delivery Reviews, 46 (2001),
pp. 3-26
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import PreDockedCompound
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
import rdkit.Chem.Lipinski as Lipinski  # type: ignore
import rdkit.Chem.Crippen as Crippen  # type: ignore
import rdkit.Chem.Descriptors as Descriptors  # type: ignore
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class LipinskiLenientFilter(SmilesFilterBase):
    """
    Implements a Lenient Lipinski filter for oral drug availability.

    This class runs a Lenient Lipinski filter to refine molecules for oral drug
    availability. It filters molecules based on Molecular Weight (MW), number of
    hydrogen donors and acceptors, and LogP value.

    This lenient implementation allows one violation of the Lipinski Rule of 5
    constraints.

    To pass the Lipinski filter, a molecule must meet the following criteria:
        - MW: Max 500 dalton
        - Number of H acceptors: Max 10
        - Number of H donors: Max 5
        - LogP: Max +5.0

    If you use the Lipinski Filter, please cite: C.A. Lipinski et al.
    Experimental and computational approaches to estimate solubility and
    permeability in drug discovery and development settings Advanced Drug
    Delivery Reviews, 46 (2001), pp. 3-26
    """

    def run_filter(self, predock_cmpd: PreDockedCompound) -> bool:
        """
        Run the Lenient Lipinski filter on a given molecule.

        This method applies the Lenient Lipinski filter criteria to determine if
        a molecule is likely to be orally available as a drug. It checks the
        molecular weight, number of hydrogen bond donors and acceptors, and LogP
        value.

        This lenient implementation allows one violation of the Lipinski Rule of
        5 constraints.

        Args:
            predock_cmpd (PreDockedCompound): A PreDockedCompound to be tested.

        Returns:
            bool: True if the molecule passes the filter (allows up to one
                violation), False otherwise.
        """
        mol = self.predock_cmpd_to_rdkit_mol(predock_cmpd)
        if mol is None:
            return False

        violation_counter = 0

        exact_mwt = Descriptors.ExactMolWt(mol)
        if exact_mwt > 500:
            violation_counter += 1

        num_hydrogen_bond_donors = Lipinski.NumHDonors(mol)
        if num_hydrogen_bond_donors > 5:
            violation_counter += 1

        num_hydrogen_bond_acceptors = Lipinski.NumHAcceptors(mol)
        if num_hydrogen_bond_acceptors > 10:
            violation_counter += 1
        mol_log_p = Crippen.MolLogP(mol)
        if mol_log_p > 5:
            violation_counter += 1

        if violation_counter < 2:
            return True

        # Failed more than two filters
        return False

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the Lenient
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
        Lenient implementation means a ligand may fail all but one requirement and still passes.",
                )
            ],
        )
