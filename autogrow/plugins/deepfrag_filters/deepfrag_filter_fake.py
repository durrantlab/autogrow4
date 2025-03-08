"""
DeepFrag plugin calculating fake fingerprints.
"""
import __future__

import numpy as np
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilterBase


class DeepFragFilterFake(DeepFragFilterBase):
    """
    DeepFrag plugin using fake fingerprints.
    """

    def get_prediction_for_parent_receptor(self, parent_mol, receptor, branching_point):
        """
        Calculating random fingerprints.

        Args:
            parent_mol: RDKit molecule representing the parent interacting with the receptor.
            receptor: .pdb file containing the receptor.
            branching_point: coordinates of the branching point.

        Returns:
           Numpy array containing fake fingerprints.
        """
        return np.random.rand(2048)

    def get_fingerprints_for_fragment(self, fragment):
        """
        Calculating a random binary fingerprints.

        Args:
            fragment: RDKit molecule representing a fragment.

        Returns:
           Numpy array containing fake binary fingerprints.
        """
        return np.random.randint(2, size=2048)
