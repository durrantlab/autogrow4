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
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound
from autogrow.plugins.registry_base import plugin_managers


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

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the VandeWaterbeemd filter on a given molecule.

        This method applies the VandeWaterbeemd filter criteria to determine if
        a molecule is likely to be blood-brain barrier permeable. It checks the
        molecular weight and polar surface area of the molecule.

        Args:
            cmpd (Compound): A Compound to be tested.

        Returns:
            bool: True if the molecule passes all filter criteria (MW < 450
                dalton and PSA < 90 A^2), False otherwise.
        """
        mol = self.cmpd_to_mol(cmpd)
        if mol is None:
            return False

        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        exact_mwt = chemtoolkit.descriptors_exact_mol_wt(mol)
        if exact_mwt >= 450:
            return False
        psa = chemtoolkit.molsurf_tpsa(mol)
        if psa >= 90:
            return False

        # passes everything
        return True

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the
        VandeWaterbeemd Filter. It allows users to enable the filter via
        command-line options.

        Returns:
            List[ArgumentVars]: A list of ArgumentVars objects defining the
            command-line arguments.
        """
        return [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable the Van de Waterbeemd filter to screen for blood-brain barrier (BBB) permeability. This filter assesses molecules based on their molecular weight and polar surface area (PSA) to predict their ability to cross the BBB.",
            )
        ]