"""VandeWaterbeemd Filter
This runs a VandeWaterbeemd filter for drugs which are likely to be blood
brain barrier permeable. VandeWaterbeemd filter filters molecules for
Molecular weight (MW), and Polar Sureface Area (PSA).

To pass the VandeWaterbeemd filter a ligand must have:
    Molecular Weight: less than 450 dalton
    Polar Sureface Area: less than 90 A^2

If you use the Van de Waterbeemd Filter please cite: Van de Waterbeemd, Han:
et al Estimation of Dlood-Brain Barrier Crossing of Drugs Using Molecular Size
and Shape, and H-Bonding Descriptors. Journal of Drug Targeting (1998), 6(2),
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
    This runs a VandeWaterbeemd filter for drugs which are likely to be blood
    brain barrier permeable. VandeWaterbeemd filter filters molecules for
    Molecular weight (MW), and Polar Sureface Area (PSA).

    To pass the VandeWaterbeemd filter a ligand must have:
        Molecular Weight: less than 450 dalton
        Polar Sureface Area: less than 90 A^2

    If you use the Van de Waterbeemd Filter please cite: Van de Waterbeemd,
    Han: et al Estimation of Dlood-Brain Barrier Crossing of Drugs Using
    Molecular Size and Shape, and H-Bonding Descriptors. Journal of Drug
    Targeting (1998), 6(2), 151-165.

    Inputs:
    :param class ParentFilter: a parent class to initialize off
    Returns:
    :returns: bool bool: True if the mol passes the filter; False if it fails the filter
    """

    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        This runs a VandeWaterbeemd filter for drugs which are likely to be
        blood brain barrier permeable. VandeWaterbeemd filter filters
        molecules for Molecular weight (MW), and Polar Sureface Area (PSA).

        To pass the VandeWaterbeemd filter a ligand must have:
            Molecular Weight: less than 450 dalton
            Polar Sureface Area: less than 90 A^2

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be
            tested if it passes the filters
        Returns:
        :returns: bool bool: True if the mol passes the filter; False if it
            fails the filter
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
        """Add command-line arguments required by the plugin."""
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
