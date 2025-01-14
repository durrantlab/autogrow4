from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING

from autogrow.plugins.pose_filters import PoseFilterPluginManager

if TYPE_CHECKING:
    from autogrow.plugins.chem_toolkit import ChemToolkitPluginManager
    from autogrow.plugins.crossover import CrossoverPluginManager
    from autogrow.plugins.docking import DockingPluginManager
    from autogrow.plugins.mutation import MutationPluginManager
    from autogrow.plugins.selectors import SelectorPluginManager
    from autogrow.plugins.shell_parallelizer import ShellParallelizerPluginManager
    from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfPluginManager
    from autogrow.plugins.smiles_filters import SmilesFilterPluginManager


@dataclass
class PluginManagerRegistry:
    """Global registry providing access to initialized AutoGrow plugin managers.
    
    This class serves as a centralized registry storing references to all initialized
    plugin managers. It provides a singleton-style access point that any part of the
    system can import to access plugin functionality without creating circular imports.
    
    The global instance of this registry is created as 'plugin_managers' in this
    module and should be imported by other modules needing access to plugins.
    """

    SmilesFilter: Optional["SmilesFilterPluginManager"] = None
    Selector: Optional["SelectorPluginManager"] = None
    Docking: Optional["DockingPluginManager"] = None
    Mutation: Optional["MutationPluginManager"] = None
    Crossover: Optional["CrossoverPluginManager"] = None
    SmiTo3DSdf: Optional["SmiTo3DSdfPluginManager"] = None
    ShellParallelizer: Optional["ShellParallelizerPluginManager"] = None
    ChemToolkit: Optional["ChemToolkitPluginManager"] = None
    PoseFilter: Optional["PoseFilterPluginManager"] = None

plugin_managers = PluginManagerRegistry()
