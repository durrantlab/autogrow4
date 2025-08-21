"""
DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
fingerprints for receptor-parent pairs.
"""
import __future__

from autogrow.utils.logging import log_warning
import rdkit
import rdkit.Chem as Chem
import numpy as np
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilter
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


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
        try:
            fp = Chem.rdmolops.RDKFingerprint(fragment, maxPath=10, fpSize=2048)
            n_fp = list(map(int, list(fp.ToBitString())))
            return np.array(n_fp)
        except Exception as e:
            log_warning(f"Error calculating RDKit fingerprint: {e}")
            return np.zeros(2048)

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.
        This method defines the command-line arguments specific to the
        DeepFrag RDKit filter.
        Returns:
        Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                        group name and a list of ArgumentVars objects defining the
                        command-line arguments.
        """
        group_name, parent_args = super().add_arguments()
        child_args = [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable the RDKit-based DeepFrag filter. This filter combines a deep learning model for analyzing the binding pocket with traditional RDKit-based fingerprinting for the appended fragment.",
            )
        ]
        return group_name, child_args + parent_args