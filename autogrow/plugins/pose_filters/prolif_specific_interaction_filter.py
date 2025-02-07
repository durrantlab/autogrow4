"""
Class for modeling filters of specific interactions.

This module uses the ProLIF library to get the interactions.

ProLIF library is available at https://github.com/chemosim-lab/ProLIF
"""
import __future__

from abc import abstractmethod
from autogrow.plugins.pose_filters.prolif_filter import ProLIFFilter
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
import rdkit.Chem.Descriptors as Descriptors  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class SpecificInteractionFilter(ProLIFFilter):
    """
    Class for modeling filters of specific interactions.

    This module uses the ProLIF library to get the interactions.

    ProLIF library is available at https://github.com/chemosim-lab/ProLIF
    """

    def run_filter(self, **kwargs) -> bool:
        """
        Run the ProLIF filter on a given molecule against a given receptor.

        This method is a generic implementation to filter out a docked compound that not contains specific interactions
        regarding a receptor.

        Args:
        **kwargs:a dictionary of arguments to pass to the plugin. It must contain the path to the
        receptor (receptor_path), a list containing Compound objects that represent docked molecules, and a
        dictionary containing the input parameters specified at the command line

        Returns:
            bool: True if the molecule passes the filter (allows up to one
                violation), False otherwise.
        """
        try:
            counter = 0
            _, df = self._compute_interaction_fingerprints(kwargs["receptor"], kwargs["docked_cmpd"])
            specific_sites = kwargs["docking_plugin_manager_params"][self.name].split(',')
            specific_interaction = self.get_interaction_type()
            for column_name in df.columns:
                if specific_interaction in column_name[2]:
                    for specific_site in specific_sites:
                        if specific_site in column_name[1]:
                            counter = counter + 1
                            break
            return counter == len(specific_sites)
        except:
            return False

    @abstractmethod
    def get_interaction_type(self) -> str:
        """
        Abstract method that should be implemented for child classes to return the interaction to be considered.

        Returns:
            string: the interaction.
        """
        pass
