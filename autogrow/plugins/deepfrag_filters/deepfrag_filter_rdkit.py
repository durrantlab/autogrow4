"""
DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
fingerprints for receptor-parent pairs.
"""
import __future__

import rdkit
import rdkit.Chem as Chem
import numpy as np
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilter

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
        except:
            return np.zeros(2048)
