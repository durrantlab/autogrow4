"""
Ghose Filter for drug-likeness.

This module implements a modified Ghose filter for drug-likeness. It filters
molecules based on Molecular Weight (MW), number of atoms, and LogP value. The
upper bound of MW is relaxed from 480Da to 500Da, making it less restrictive and
compatible with Lipinski's rule. This also matches AutoGrow 3.1.3's
implementation.

The filter protonates molecules as hydrogens affect atom count. This
implementation counts hydrogens against the total number of atoms.

To pass the filter, a molecule must meet the following criteria:
    - MW between 160 and 500 dalton
    - Number of Atoms: between 20 and 70
    - LogP between -0.4 and +5.6

If you use the Ghose Filter, please cite: A.K. Ghose et al. A knowledge-based
approach in designing combinatorial or medicinal chemistry libraries for drug
discovery. 1. A qualitative and quantitative characterization of known drug
databases Journal of Combinatorial Chemistry, 1 (1999), pp. 55-68
"""
import __future__

import copy

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


class GhoseModifiedFilter(SmilesFilterBase):
    """
    Implements a modified Ghose filter for drug-likeness.

    This class runs a Ghose filter that checks molecules for drug-likeness based
    on Molecular Weight (MW), number of atoms, and LogP value. The upper bound
    of MW is relaxed from 480Da to 500Da, making it less restrictive and
    compatible with Lipinski's rule. This also matches AutoGrow 3.1.3's
    implementation.

    The filter protonates molecules as hydrogens affect atom count. This
    implementation counts hydrogens against the total number of atoms.

    To pass the filter, a molecule must meet the following criteria:
        - MW between 160 and 500 dalton
        - Number of Atoms: between 20 and 70
        - LogP between -0.4 and +5.6

    If you use the Ghose Filter, please cite: A.K. Ghose et al. A
    knowledge-based approach in designing combinatorial or medicinal chemistry
    libraries for drug discovery. 1. A qualitative and quantitative
    characterization of known drug databases Journal of Combinatorial Chemistry,
    1 (1999), pp. 55-68
    """

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the modified Ghose filter on a given molecule.

        This method applies the Ghose filter criteria to determine if a molecule
        is drug-like. It checks the molecular weight, number of atoms, molar
        refractivity, and molar LogP.

        The molecule is protonated before analysis as hydrogens affect the atom
        count. This implementation counts hydrogens against the total number of
        atoms.

        Args:
            cmpd (Compound): A Compound to be tested.

        Returns:
            bool: True if the molecule passes all filter criteria, False
                otherwise.
        """
        mol = self.predock_cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        copy_mol = copy.deepcopy(mol)
        copy_mol = Chem.AddHs(copy_mol)
        exact_mwt = Descriptors.ExactMolWt(copy_mol)
        if (exact_mwt < 160) or (exact_mwt > 500):
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

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the Ghose
        Modified Filter. It allows users to enable the filter via command-line
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
                    help="Ghose filters for drug-likeliness; filters by molecular weight,\
        logP and number of atoms. This is the same as the GhoseFilter, but \
        the upper-bound of the molecular weight restrict is loosened from \
        480Da to 500Da. This is intended to be run with Lipinski Filter and \
        to match AutoGrow 3's Ghose Filter.",
                )
            ],
        )
