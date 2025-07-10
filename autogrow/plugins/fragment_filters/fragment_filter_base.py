"""
Defines the base class for fragment filter plugins.
"""
from abc import abstractmethod
from typing import Any
from autogrow.plugins.plugin_base import PluginBase

class FragmentFilterBase(PluginBase):
    """
    Base class for fragment filters.
    This class defines the interface for all fragment filter plugins.
    """

    def run(self, **kwargs) -> bool:
        """
        Execute the filter plugin with the provided molecule.
        Args:
            **kwargs: Keyword arguments containing filter parameters, must
                include: mol (rdkit.Chem.rdchem.Mol): The molecule to be
                filtered.
        Returns:
            bool: True if the molecule passes the filter criteria, False
                otherwise.
        """
        return self.run_filter(kwargs["mol"])

    @abstractmethod
    def run_filter(self, mol: Any) -> bool:
        """
        Run the filter on a single RDKit Mol object.
        This method must be implemented by each subclass.
        Args:
            mol (rdkit.Chem.rdchem.Mol): A molecule to filter.
        Returns:
            bool: True if the molecule passes the filter, False if it fails.
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        PASS