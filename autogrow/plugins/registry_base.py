from dataclasses import dataclass
from typing import Dict, List, Optional, Type, Tuple, Any
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase

class PluginManagerRegistry:
    """Centralized registry for all plugin managers"""
    
    def __init__(self):
        self._managers: Dict[str, Optional[PluginManagerBase]] = {}
        self._initialized = False

    def initialize(self):
        """Initialize the registry if not already done"""
        if self._initialized:
            return
            
        self._register_plugin_types()
        self._initialized = True

    def _register_plugin_types(self) -> None:
        """Register all plugin types in the registry"""
        for name, _, _ in self.get_plugin_types():
            self._managers[name] = None
            # Add property accessor if it doesn't exist
            if not hasattr(self.__class__, name):
                setattr(self.__class__, name, property(
                    lambda self, n=name: self._managers[n]
                ))
    
    @classmethod
    def get_plugin_types(cls) -> List[Tuple[str, Type[PluginManagerBase], Type[PluginBase]]]:
        """Returns the definitive list of all plugin types in the system"""
        # Import here to avoid circular imports
        from autogrow.plugins.plugin_types import PLUGIN_TYPES
        return PLUGIN_TYPES

    def setup_plugin_managers(self, params: dict) -> None:
        """Initialize and setup all plugin managers"""
        self.initialize()  # Ensure registry is initialized
        
        # First initialize all managers
        for name, manager_class, base_class in self.get_plugin_types():
            self._managers[name] = manager_class(base_class)

        # Then setup ChemToolkit first
        self._managers["ChemToolkit"].setup_plugin_manager(params, self)

        # Then setup all others
        for name, _, _ in self.get_plugin_types():
            if name != "ChemToolkit":
                self._managers[name].setup_plugin_manager(params, self)

    @classmethod
    def get_all_plugin_managers(cls) -> List[PluginManagerBase]:
        """Get list of all plugin manager instances for argument registration"""
        return [
            manager_class(base_class)
            for _, manager_class, base_class in cls.get_plugin_types()
        ]

# Create the global instance
plugin_managers = PluginManagerRegistry()