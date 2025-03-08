"""
DeepFrag plugin using fake fingerprints.
"""
import __future__

import rdkit
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilterBase
from autogrow.deepfrag_integration.DeepFragIntegration import DeepFragFake

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class DeepFragFilterFake(DeepFragFilterBase):
    """
    DeepFrag plugin using fake fingerprints.
    """

    def get_concrete_deepfrag(self):
        return DeepFragFake()
