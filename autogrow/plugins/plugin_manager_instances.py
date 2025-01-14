# from dataclasses import dataclass
# from typing import Dict, List, Optional, Type, Tuple
# from autogrow.plugins.plugin_base import PluginBase
# from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.plugins.registry_base import plugin_managers


# @dataclass
# class PluginManagerRegistry:
#     """Global registry providing access to initialized AutoGrow plugin managers.
    
#     This class serves as a centralized registry storing references to all initialized
#     plugin managers. It provides a singleton-style access point that any part of the
#     system can import to access plugin functionality without creating circular imports.
    
#     The global instance of this registry is created as 'plugin_managers' in this
#     module and should be imported by other modules needing access to plugins.
#     """

#     def __init__(self):
#         # Initialize all plugin managers as None
#         self._managers: Dict[str, Optional[PluginManagerBase]] = {}
        
#         # Register all plugin types
#         self._register_plugin_types()

#     @classmethod
#     def get_plugin_types(cls) -> List[Tuple[str, Type[PluginManagerBase], Type[PluginBase]]]:
#         """
#         Returns the definitive list of all plugin types in the system.
#         This is the single source of truth for available plugins.
        
#         Returns:
#             List of tuples containing:
#                 - Name of the plugin type
#                 - Plugin manager class
#                 - Plugin base class
#         """
#         # Import here to avoid circular imports
#         from autogrow.plugins.chem_toolkit import ChemToolkitPluginManager, ChemToolkitBase
#         from autogrow.plugins.smiles_filters import SmilesFilterPluginManager, SmilesFilterBase
#         from autogrow.plugins.selectors import SelectorPluginManager, SelectorBase
#         from autogrow.plugins.docking import DockingPluginManager, DockingBase
#         from autogrow.plugins.mutation import MutationPluginManager, MutationBase
#         from autogrow.plugins.crossover import CrossoverPluginManager, CrossoverBase
#         from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfPluginManager, SmiTo3DSdfBase
#         from autogrow.plugins.shell_parallelizer import ShellParallelizerPluginManager, ShellParallelizerBase
#         from autogrow.plugins.pose_filters import PoseFilterPluginManager, PoseFilterBase

#         return [
#             # ChemToolkit must be first as others depend on it
#             ("ChemToolkit", ChemToolkitPluginManager, ChemToolkitBase),
#             ("SmilesFilter", SmilesFilterPluginManager, SmilesFilterBase),
#             ("Selector", SelectorPluginManager, SelectorBase),
#             ("Docking", DockingPluginManager, DockingBase),
#             ("Mutation", MutationPluginManager, MutationBase),
#             ("Crossover", CrossoverPluginManager, CrossoverBase),
#             ("SmiTo3DSdf", SmiTo3DSdfPluginManager, SmiTo3DSdfBase),
#             ("ShellParallelizer", ShellParallelizerPluginManager, ShellParallelizerBase),
#             ("PoseFilter", PoseFilterPluginManager, PoseFilterBase),
#         ]

#     def _register_plugin_types(self) -> None:
#         """Register all plugin types in the registry"""
#         for name, _, _ in self.get_plugin_types():
#             self._managers[name] = None
#             # Add property accessor
#             setattr(self.__class__, name, property(
#                 lambda self, n=name: self._managers[n]
#             ))

#     def setup_plugin_managers(self, params: dict) -> None:
#         """Initialize and setup all plugin managers"""
#         # First initialize all managers
#         for name, manager_class, base_class in self.get_plugin_types():
#             self._managers[name] = manager_class(base_class)

#         # Then setup ChemToolkit first
#         assert self._managers["ChemToolkit"] is not None, "ChemToolkit not set"
#         self._managers["ChemToolkit"].setup_plugin_manager(params, self)

#         # Then setup all others
#         for name, _, _ in self.get_plugin_types():
#             if name != "ChemToolkit":
#                 assert self._managers[name] is not None, f"{name} not set"
#                 self._managers[name].setup_plugin_manager(params, self)

#     @classmethod
#     def get_all_plugin_managers(cls) -> List[PluginManagerBase]:
#         """
#         Get list of all plugin manager instances for argument registration.
#         Used by argparser.py.
#         """
#         return [
#             manager_class(base_class)
#             for _, manager_class, base_class in cls.get_plugin_types()
#         ]

# plugin_managers = PluginManagerRegistry()
