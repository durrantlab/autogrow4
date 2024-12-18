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
from autogrow.plugins.selectors import SelectorBase, SelectorPluginManager
from autogrow.plugins.shell_parallelizer import (
    ShellParallelizerBase,
    ShellParallelizerPluginManager,
)
from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase, SmiTo3DSdfPluginManager
from autogrow.plugins.smiles_filters import SmilesFilterBase, SmilesFilterPluginManager


@dataclass
class PluginManagers:
    """
    Container class for all plugin managers in the AutoGrow system.

    This class serves as a central point of access for all plugin managers,
    making it easy to coordinate their setup and usage. It maintains one
    instance of each type of plugin manager.
    """

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
    ChemToolkit: ChemToolkitPluginManager = ChemToolkitPluginManager(ChemToolkitBase)


plugin_managers = PluginManagers()


def setup_plugin_managers(params: Dict[str, Any]):
    """
    Set up all plugin managers with the provided parameters.

    This function iterates through all plugin managers defined in the
    PluginManagers class and sets up each one with the provided parameters. Each
    manager is initialized with access to all other managers to enable
    inter-plugin coordination.

    Args:
        params (Dict[str, Any]): Configuration parameters for all plugin
            managers.
    """
    # Iterate through all plugin managers in the PluginManagers class and setup each one
    for name, plugin_manager in PluginManagers.__annotations__.items():
        getattr(plugin_managers, name).setup_plugin_manager(params, plugin_managers)
