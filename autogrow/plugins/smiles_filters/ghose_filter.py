"""
Implements the Ghose filter for drug-likeness screening in AutoGrow.

This module provides the GhoseFilter class that filters ligands using the Ghose
filter for drug-likeness. The filter screens molecules based on Molecular Weight
(MW), number of atoms, and LogP value.

To pass the filter, a molecule must meet the following criteria:
    - MW between 160 and 480 daltons
    - Number of Atoms: between 20 and 70
    - LogP between -0.4 and +5.6

Note: This implementation protonates the molecule before filtering because
hydrogens affect atom count. The Ghose filter implementation counts hydrogens
against the total number of atoms.

If using the Ghose Filter, please cite: A.K. Ghose et al. A knowledge-based
approach in designing combinatorial or medicinal chemistry libraries for drug
discovery. 1. A qualitative and quantitative characterization of known drug
databases Journal of Combinatorial Chemistry, 1 (1999), pp. 55-68
"""
import __future__

import copy

from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.registry_base import plugin_managers


class GhoseFilter(SmilesFilterBase):
    """
    A filter that implements the Ghose filter for drug-likeness screening.

    This class filters molecules based on Molecular Weight (MW), number of
    atoms, and LogP value according to the Ghose criteria for drug-likeness.

    Note:
        This implementation protonates the molecule before filtering because
        hydrogens affect atom count. The Ghose filter implementation counts
        hydrogens against the total number of atoms.
    """

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the Ghose filter on a given molecule.

        This method applies the Ghose filter criteria to determine if a molecule
        is drug-like. It checks the molecular weight, number of atoms, molar
        refractivity, and molar LogP.

        Args:
            cmpd (Compound): A Compound to be tested.

        Returns:
            bool: True if the molecule passes the filter; False if it fails.

        Note:
            This method creates a copy of the input molecule and adds hydrogens
            to it before applying the filter. This ensures that hydrogens are
            counted in the total atom count without affecting other filters.
        """
        mol = self.cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        # Make a copy of the mol so we can AddHs without affecting other filters
        # number of atoms is altered by the presence/absence of hydrogens.
        # Our Ghose filter counts hydrogenss towards atom count
        copy_mol = copy.deepcopy(mol)
        copy_mol = chemtoolkit.add_hs(copy_mol)
        exact_mwt = chemtoolkit.descriptors_exact_mol_wt(copy_mol)
        if (exact_mwt < 160) or (exact_mwt > 480):
            return False

        num_atoms = chemtoolkit.get_num_atoms(copy_mol)
        if (num_atoms < 20) or (num_atoms > 70):
            return False

        # molar Refractivity
        MolMR = chemtoolkit.crippen_mol_mr(copy_mol)
        if (MolMR < 40) or (MolMR > 130):
            return False

        # molar LogP
        mol_log_p = chemtoolkit.crippen_mol_log_p(copy_mol)
        if (mol_log_p < -0.4) or (mol_log_p > 5.6):
            return False

        # passed all filters
        return True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the Ghose filter.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("SMILES Filters")
                - A list with one ArgumentVars object defining the argument to
                  enable the Ghose filter
        """
        return (
            "SMILES Filters",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Ghose filters for drug-likeliness; filters by molecular weight,\
        logP and number of atoms.",
                )
            ],
        )
