"""
Class for filtering out a molecule that does not have hydrophobic interactions against a given receptor at specific
residues.

This module uses the ProLIF library to get the hydrophobic interactions. To pass the filter, a molecule must have each
hydrophobic interaction specified as input.

ProLIF library is available at https://github.com/chemosim-lab/ProLIF
"""

from autogrow.plugins.pose_filters.prolif_specific_interaction_filter import SpecificInteractionFilter
import prolif
import rdkit  # type: ignore
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class HydrophobicInteractionFilter(SpecificInteractionFilter):
    """
    Class for filtering out a molecule that does not have hydrophobic interactions against a given receptor at specific
    residues.

    This module uses the ProLIF library to get the hydrophobic interactions. To pass the filter, a molecule must have
    each hydrophobic interaction specified as input.

    ProLIF library is available at https://github.com/chemosim-lab/ProLIF
    """
    fingerprints = prolif.Fingerprint(["Hydrophobic"])

    def get_interaction_type(self) -> str:
        return "Hydrophobic"

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the ProLIF Filter.
        It allows users to enable the filter via command-line options.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of ArgumentVars objects defining the
                command-line arguments.
        """
        return (
            "Pose Filters",
            [
                ArgumentVars(
                    name=self.name,
                    type=str,
                    default=False,
                    help="Comma-separated list of specific residues in a receptor that should have hydrophobic "
                         "interactions with respect to a docked molecule. For example, --HydrophobicInteractionFilter "
                         "HIS862,TYR896,TYR907",
                )
            ],
        )
