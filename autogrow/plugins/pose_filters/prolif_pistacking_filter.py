"""
Class for filtering out a molecule that does not have PiStacking interactions against a given receptor at specific
residues.

This module uses the ProLIF library to get the PiStacking interactions. To pass the filter, a molecule must have each
PiStacking interaction specified as input.

ProLIF library is available at https://github.com/chemosim-lab/ProLIF
"""

from autogrow.plugins.pose_filters.prolif_specific_interaction_filter import SpecificInteractionFilter
import prolif
import rdkit  # type: ignore
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class PiStackingInteractionFilter(SpecificInteractionFilter):
    """
    Class for filtering out a molecule that does not have PiStacking interactions against a given receptor at specific
    residues.

    This module uses the ProLIF library to get the PiStacking interactions. To pass the filter, a molecule must have
    each PiStacking interaction specified as input.

    ProLIF library is available at https://github.com/chemosim-lab/ProLIF
    """
    fingerprints = prolif.Fingerprint(["PiStacking"])

    def get_interaction_type(self) -> str:
        return "PiStacking"

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
                    help="Enable the Pi-Stacking Interaction Filter. Provide a comma-separated list of specific residues in the receptor (e.g., 'TYR907,PHE123') that must form pi-stacking interactions with the docked molecule. The molecule will be filtered out if it does not form all specified interactions. For example, --PiStackingInteractionFilter TYR907",
                )
            ],
        )
