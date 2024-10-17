from dataclasses import dataclass
from typing import TYPE_CHECKING

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


plugin_managers = PluginManagers()


def setup_plugin_managers(params):
    # Iterate through all plugin managers in the PluginManagers class and setup each one
    for name, plugin_manager in PluginManagers.__annotations__.items():
        getattr(plugin_managers, name).setup_plugin_manager(params, plugin_managers)
