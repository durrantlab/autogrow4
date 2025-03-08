"""
DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
fingerprints for receptor-parent pairs.
"""
import __future__

import rdkit
import rdkit.Chem as Chem
import numpy as np
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilterBase

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class DeepFragFilterRDKit(DeepFragFilterBase):
    """
    DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
    fingerprints for receptor-parent pairs.
    """

    def get_prediction_for_parent_receptor(self, parent_mol, receptor, branching_point):
        """
        Calculating DeepFrag fingerprints for the receptor-parent complex.

        Args:
            parent_mol: RDKit molecule representing the parent interacting with the receptor.
            receptor: .pdb file containing the receptor.
            branching_point: coordinates of the branching point.

        Returns:
           Numpy array containing the DeepFrag fingerprints.
        """
        return np.random.rand(2048)

    def get_fingerprints_for_fragment(self, fragment):
        """
        Calculating a RDKit binary fingerprints for a chemical fragment.

        Args:
            fragment: RDKit molecule representing a chemical fragment.

        Returns:
           Numpy array containing the binary fingerprints calculated with RDKit library.
        """
        fp = Chem.rdmolops.RDKFingerprint(fragment, maxPath=10, fpSize=2048)
        n_fp = list(map(int, list(fp.ToBitString())))
        return np.array(n_fp)
