"""
Central management of AutoGrow's plugin system.

This module provides centralized access and coordination for all of AutoGrow's
plugin managers. It defines a PluginManagers class that maintains instances of
all plugin managers and provides utilities for their setup and coordination.

A single global instance of PluginManagers is created and maintained to provide
consistent access to all plugin managers throughout the application.
"""
from dataclasses import dataclass
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


@dataclass
class PluginManagerFactory:
    """Factory class for creating and initializing AutoGrow plugin managers.
    
    This class serves as the factory that instantiates and configures all plugin 
    managers in the system. It creates fresh instances of each plugin manager type 
    (ChemToolkit, SmilesFilter, etc.) and handles their initialization and setup.

    The factory creates plugin managers but does not store global references to them - 
    that is handled by the PluginManagerRegistry. This separation avoids circular
    dependencies since plugins can access other plugins through the registry without
    importing the factory directly.
    """

    ChemToolkit: ChemToolkitPluginManager = ChemToolkitPluginManager(ChemToolkitBase)
    SmilesFilter: SmilesFilterPluginManager = SmilesFilterPluginManager(
        SmilesFilterBase
    )
    Selector: SelectorPluginManager = SelectorPluginManager(SelectorBase)
    Docking: DockingPluginManager = DockingPluginManager(DockingBase)
    Mutation: MutationPluginManager = MutationPluginManager(MutationBase)
    Crossover: CrossoverPluginManager = CrossoverPluginManager(CrossoverBase)
    SmiTo3DSdf: SmiTo3DSdfPluginManager = SmiTo3DSdfPluginManager(SmiTo3DSdfBase)
    ShellParallelizer: ShellParallelizerPluginManager = ShellParallelizerPluginManager(
        ShellParallelizerBase
    )
    PoseFilter: PoseFilterPluginManager = PoseFilterPluginManager(
        PoseFilterBase
    )


def setup_plugin_managers(params: Dict[str, Any]):
    """Set up all plugin managers with the provided parameters.
    
    This function creates and initializes all plugin managers through the factory,
    then populates the global registry with references to them.
    
    Args:
        params (Dict[str, Any]): Configuration parameters for all plugin managers.
    
    Raises:
        AssertionError: If the factory and registry classes have mismatched attributes
    """
    # First verify that factory and registry have matching attributes
    factory_attrs = set(PluginManagerFactory.__annotations__.keys())
    registry_attrs = set(PluginManagerRegistry.__annotations__.keys())
    
    if factory_attrs != registry_attrs:
        extra_in_factory = factory_attrs - registry_attrs
        extra_in_registry = registry_attrs - factory_attrs
        msg = "Plugin manager factory and registry have mismatched attributes.\n"
        if extra_in_factory:
            msg += f"Attributes in factory but not registry: {extra_in_factory}\n"
        if extra_in_registry:
            msg += f"Attributes in registry but not factory: {extra_in_registry}"
        raise AssertionError(msg)

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