"""
DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
fingerprints for receptor-parent pairs.
"""
import __future__

from autogrow.utils.logging import log_warning
import numpy as np
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilter
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.registry_base import plugin_managers


class DeepFragFilterRDKit(DeepFragFilter):
    """
    DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
    fingerprints for receptor-parent pairs.
    """

    def get_fingerprints_for_fragment(self, fragment):
        """
        Calculating RDKit binary fingerprints for a chemical fragment.

        Args:
            fragment: RDKit molecule representing a chemical fragment.

        Returns:
           Numpy array containing the binary fingerprints calculated with RDKit library.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        try:
            fp = chemtoolkit.get_rdk_fingerprint(fragment, max_path=10, fp_size=2048)
            n_fp = list(map(int, list(chemtoolkit.bit_vect_to_bit_string(fp))))
            return np.array(n_fp)
        except Exception as e:
            log_warning(f"Error calculating RDKit fingerprint: {e}")
            return np.zeros(2048)

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments required by the plugin.
        This method defines the command-line arguments specific to the
        DeepFrag RDKit filter.
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
                help="Enable the RDKit-based DeepFrag filter. This filter combines a deep learning model for analyzing the binding pocket with traditional RDKit-based fingerprinting for the appended fragment.",
            )
        ]
        return child_args + parent_args

