"""
Mozziconacci Filter for drug-likeness.

This module implements a Mozziconacci filter to refine molecules for
drug-likeness. It filters molecules based on the number of rotatable bonds,
rings, oxygens, nitrogens, and halogens.

To pass the filter, a molecule must meet the following criteria:
    - Number of Rotatable bonds: Max 15
    - Number of Rings: Max 6
    - Number of Oxygens: Min 1
    - Number of Nitrogens: Min 1
    - Number of Halogens: Max 7

If you use the Mozziconacci Filter, please cite: Mozziconacci, J. C. et al.
Preparation of a Molecular Database from a Set of 2 Million Compounds for
Virtual Screening Applications: Gathering, Structural Analysis and Filtering.
9th Electronic Computational Chemistry Conference, World Wide Web, March (2003).
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.registry_base import plugin_managers


class MozziconacciFilter(SmilesFilterBase):
    """
    Implements a Mozziconacci filter for drug-likeness.

    This class runs a Mozziconacci filter to refine molecules for drug-likeness.
    It filters molecules based on the number of rotatable bonds, rings, oxygens,
    nitrogens, and halogens.

    To pass the filter, a molecule must meet the following criteria:
        - Number of Rotatable bonds: Max 15
        - Number of Rings: Max 6
        - Number of Oxygens: Min 1
        - Number of Nitrogens: Min 1
        - Number of Halogens: Max 7

    If you use the Mozziconacci Filter, please cite: Mozziconacci, J. C. et al.
    Preparation of a Molecular Database from a Set of 2 Million Compounds for
    Virtual Screening Applications: Gathering, Structural Analysis and
    Filtering. 9th Electronic Computational Chemistry Conference, World Wide
    Web, March (2003).
    """

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the Mozziconacci filter on a given molecule.

        This method applies the Mozziconacci filter criteria to determine if a
        molecule is likely to be drug-like. It checks the number of rotatable
        bonds, rings, oxygens, nitrogens, and halogens.

        Args:
            cmpd (Compound): A Compound to be tested.

        Returns:
            bool: True if the molecule passes all filter criteria, False
                otherwise.
        """
        mol = self.cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        halogen = chemtoolkit.mol_from_smarts("[*;#9,#17,#35,#53,#85]")
        number_of_halogens = len(
            chemtoolkit.get_substruct_matches(mol, halogen, max_matches=8)
        )
        if number_of_halogens > 7:
            return False

        oxygen = chemtoolkit.mol_from_smarts("[#8]")
        number_of_oxygens = len(
            chemtoolkit.get_substruct_matches(mol, oxygen, max_matches=2)
        )
        if number_of_oxygens < 1:
            return False

        nitrogen = chemtoolkit.mol_from_smarts("[#7]")
        number_of_nitrogen = len(
            chemtoolkit.get_substruct_matches(mol, nitrogen, max_matches=2)
        )
        if number_of_nitrogen < 1:
            return False

        num_rotatable_bonds = chemtoolkit.lipinski_num_rotatable_bonds(mol)
        if num_rotatable_bonds > 15:
            return False

        ring_count = chemtoolkit.rdmolops_get_sssr(mol)
        if ring_count > 6:
            return False

        # Passes everything
        return True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the
        Mozziconacci Filter. It allows users to enable the filter via
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
                    help="Mozziconacci filters for drug-likeliness; filters by the number of \
        rotatable bonds, rings, oxygens, and halogens.",
                )
            ],
        )
