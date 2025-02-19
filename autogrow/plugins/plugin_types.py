"""
Central registry of all plugin types. This is the single place to register new plugins.
"""

# These imports are fine here as they don't import plugin_managers
from autogrow.plugins.chem_toolkit import ChemToolkitPluginManager, ChemToolkitBase
from autogrow.plugins.smiles_filters import SmilesFilterPluginManager, SmilesFilterBase
from autogrow.plugins.selectors import SelectorPluginManager, SelectorBase
from autogrow.plugins.docking import DockingPluginManager, DockingBase
from autogrow.plugins.mutation import MutationPluginManager, MutationBase
from autogrow.plugins.crossover import CrossoverPluginManager, CrossoverBase
from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfPluginManager, SmiTo3DSdfBase
from autogrow.plugins.shell_parallelizer import ShellParallelizerPluginManager, ShellParallelizerBase
from autogrow.plugins.pose_filters import PoseFilterPluginManager, PoseFilterBase
from autogrow.plugins.rescoring import RescoringPluginManager, RescoringBase

# TODO: I bet it would be possible to autogenerate this list through dynamic
# imports. Good to investigate.

# The definitive list of all plugin types
PLUGIN_TYPES = [
    # ChemToolkit must be first as others depend on it
    ("ChemToolkit", ChemToolkitPluginManager, ChemToolkitBase),

    ("SmilesFilter", SmilesFilterPluginManager, SmilesFilterBase),
    ("Selector", SelectorPluginManager, SelectorBase),
    ("Docking", DockingPluginManager, DockingBase),
    ("Mutation", MutationPluginManager, MutationBase),
    ("Crossover", CrossoverPluginManager, CrossoverBase),
    ("SmiTo3DSdf", SmiTo3DSdfPluginManager, SmiTo3DSdfBase),
    ("ShellParallelizer", ShellParallelizerPluginManager, ShellParallelizerBase),
    ("PoseFilter", PoseFilterPluginManager, PoseFilterBase),
    ("Rescoring", RescoringPluginManager, RescoringBase),
]