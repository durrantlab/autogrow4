"""
VandeWaterbeemd Filter for blood-brain barrier permeability.

This module implements the VandeWaterbeemd filter to identify drugs that are
likely to be blood-brain barrier permeable. It filters molecules based on
Molecular Weight (MW) and Polar Surface Area (PSA).

To pass the VandeWaterbeemd filter, a ligand must have:
    - Molecular Weight: less than 450 dalton
    - Polar Surface Area: less than 90 A^2

If you use the Van de Waterbeemd Filter, please cite: Van de Waterbeemd, Han: et
al. Estimation of Blood-Brain Barrier Crossing of Drugs Using Molecular Size and
Shape, and H-Bonding Descriptors. Journal of Drug Targeting (1998), 6(2),
151-165.
"""
import __future__
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.smiles_filters import SmilesFilterBase
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
import rdkit.Chem.MolSurf as MolSurf  # type: ignore
import rdkit.Chem.Descriptors as Descriptors  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class VandeWaterbeemdFilter(SmilesFilterBase):
    """
    Implements the VandeWaterbeemd filter for blood-brain barrier permeability.

    This class filters molecules based on Molecular Weight (MW) and Polar
    Surface Area (PSA) to identify drugs that are likely to be blood-brain
    barrier permeable.

    To pass the VandeWaterbeemd filter, a ligand must have:
        - Molecular Weight: less than 450 dalton
        - Polar Surface Area: less than 90 A^2

    If you use the Van de Waterbeemd Filter, please cite: Van de Waterbeemd,
    Han: et al. Estimation of Blood-Brain Barrier Crossing of Drugs Using
    Molecular Size and Shape, and H-Bonding Descriptors. Journal of Drug
    Targeting (1998), 6(2), 151-165.
    """

    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        Run the VandeWaterbeemd filter on a given molecule.

        This method applies the VandeWaterbeemd filter criteria to determine if
        a molecule is likely to be blood-brain barrier permeable. It checks the
        molecular weight and polar surface area of the molecule.

        Args:
            mol (rdkit.Chem.rdchem.Mol): An RDKit mol object to be tested
                against the filter criteria.

        Returns:
            bool: True if the molecule passes all filter criteria (MW < 450
                dalton and PSA < 90 A^2), False otherwise.
        """
        exact_mwt = Descriptors.ExactMolWt(mol)
        if exact_mwt >= 450:
            return False
        psa = MolSurf.TPSA(mol)
        if psa >= 90:
            return False

        # passes everything
        return True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the
        VandeWaterbeemd Filter. It allows users to enable the filter via
        command-line options.

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
                    help="VandeWaterbeemd filters for drug likely to be blood brain barrier permeable. \
        Filters by the number of molecular weight and Polar Sureface Area (PSA).",
                )
            ],
        )
