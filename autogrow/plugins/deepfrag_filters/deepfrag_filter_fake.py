"""
DeepFrag plugin calculating fake fingerprints.
"""
import __future__

import numpy as np
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilterBase
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars

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

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments required by the plugin.
        This method defines the command-line arguments specific to the
        DeepFrag filter.
        It allows users to enable the filter via command-line options.
        Returns:
            List[ArgumentVars]: A list of ArgumentVars objects defining the
            command-line arguments.
        """
        parent_args = super().add_arguments()
        child_args = [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable the fake DeepFrag filter for testing purposes. This filter assigns random scores and does not perform any real analysis, useful for debugging the pipeline.",
            )
        ]
        return child_args + parent_args