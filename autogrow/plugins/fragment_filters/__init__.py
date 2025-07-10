"""
Manages fragment filter plugins for AutoGrow.
This module provides the plugin manager for fragment filtering operations.
"""
from typing import Any, List, cast
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.plugins.fragment_filters.fragment_filter_base import FragmentFilterBase
from autogrow.utils.logging import log_warning

class FragmentFilterPluginManager(PluginManagerBase):
    """Manages and executes fragment filter plugins."""

    def execute(self, **kwargs) -> List:
        """
        Run all selected filters on a list of molecules.
        Args:
            **kwargs: Keyword arguments containing:
                - mols (List[rdkit.Chem.rdchem.Mol]): A list of molecules to be filtered.
        Returns:
            List[rdkit.Chem.rdchem.Mol]: A list of molecules that passed all filters.
        """
        passed_mols: List[Any] = []
        passed_mols.extend(
            mol for mol in kwargs["mols"] if self._run_all_selected_filters(mol)
        )
        return passed_mols

    def _run_all_selected_filters(self, mol: Any) -> bool:
        """
        Determine if a molecule passes all selected filters.
        Args:
            mol: An RDKit mol object to be tested.
        Returns:
            bool: True if the mol passes all the filters. False if it fails any filters.
        """
        filters_failed = 0
        for plugin_name in self.plugins:
            plugin = cast(FragmentFilterBase, self.plugins[plugin_name])
            filter_function = plugin.run
            if not filter_function(mol=mol):
                filters_failed += 1
        if filters_failed == 0:
            return True
        return False