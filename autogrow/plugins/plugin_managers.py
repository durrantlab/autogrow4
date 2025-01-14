"""
Central management of AutoGrow's plugin system.

This module provides centralized access and coordination for all of AutoGrow's
plugin managers. It defines a PluginManagers class that maintains instances of
all plugin managers and provides utilities for their setup and coordination.

A single global instance of PluginManagers is created and maintained to provide
consistent access to all plugin managers throughout the application.
"""
from typing import Any, Dict

from autogrow.plugins.chem_toolkit import ChemToolkitBase, ChemToolkitPluginManager
from autogrow.plugins.crossover import CrossoverBase, CrossoverPluginManager
from autogrow.plugins.docking import DockingBase, DockingPluginManager
from autogrow.plugins.mutation import MutationBase, MutationPluginManager
from autogrow.plugins.pose_filters import PoseFilterBase, PoseFilterPluginManager
from autogrow.plugins.selectors import SelectorBase, SelectorPluginManager
from autogrow.plugins.shell_parallelizer import (
    ShellParallelizerBase,
    ShellParallelizerPluginManager,
)
from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase, SmiTo3DSdfPluginManager
from autogrow.plugins.smiles_filters import SmilesFilterBase, SmilesFilterPluginManager
from autogrow.plugins.plugin_manager_instances import PluginManagerRegistry, plugin_managers


def setup_plugin_managers(params: Dict[str, Any]):
    """Set up all plugin managers with the provided parameters.
    
    This function creates and initializes all plugin managers through the factory,
    then populates the global registry with references to them.
    
    Args:
        params (Dict[str, Any]): Configuration parameters for all plugin managers.
    
    Raises:
        AssertionError: If the factory and registry classes have mismatched attributes
    """
    # Create instances manually since order matters (ChemToolkit must be first)
    plugin_managers.ChemToolkit = ChemToolkitPluginManager(ChemToolkitBase)
    plugin_managers.ChemToolkit.setup_plugin_manager(params, plugin_managers)

    # Initialize others after ChemToolkit
    plugin_managers.SmilesFilter = SmilesFilterPluginManager(SmilesFilterBase)
    plugin_managers.Selector = SelectorPluginManager(SelectorBase)
    plugin_managers.Docking = DockingPluginManager(DockingBase)
    plugin_managers.Mutation = MutationPluginManager(MutationBase)
    plugin_managers.Crossover = CrossoverPluginManager(CrossoverBase)
    plugin_managers.SmiTo3DSdf = SmiTo3DSdfPluginManager(SmiTo3DSdfBase)
    plugin_managers.ShellParallelizer = ShellParallelizerPluginManager(
        ShellParallelizerBase
    )
    plugin_managers.PoseFilter = PoseFilterPluginManager(PoseFilterBase)

    # Setup all plugin managers
    for name in PluginManagerRegistry.__annotations__.keys():
        if name != "ChemToolkit":  # Skip since we already set it up
            plugin_manager = getattr(plugin_managers, name)
            if plugin_manager is not None:
                plugin_manager.setup_plugin_manager(params, plugin_managers)