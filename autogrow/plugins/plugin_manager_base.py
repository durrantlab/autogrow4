"""
Plugin management system.

This module provides the base classes for AutoGrow's plugin management system,
which enables extensible, modular functionality through plugins. It handles
plugin discovery, loading, validation, setup, execution, and caching.

The system consists of two main components:
    - PluginManagerBase: Abstract base class for plugin managers that handle
      specific types of plugins (e.g., filters, mutations).
    - PluginBase: Abstract base class that all plugins must inherit from,
      defining the required interface.

Plugin managers automatically discover and load plugins from their respective
directories. Plugins can be selectively enabled via command-line arguments,
which are automatically registered during plugin loading. Results can be cached
to enable quick restart of interrupted runs.
"""

import os
import importlib
import inspect
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Type, TYPE_CHECKING
from autogrow.config.argument_vars import register_argparse_group
from autogrow.plugins.plugin_base import PluginBase

from autogrow.utils.caching import CacheManager

if TYPE_CHECKING:
    from autogrow.plugins.plugin_manager_instances import PluginManagerRegistry


class PluginManagerBase(ABC):
    """
    Abstract base class for plugin managers in the AutoGrow system.

    This class handles plugin loading, setup, execution, and caching. It
    provides a framework for managing different types of plugins and their
    lifecycle.

    Attributes:
        plugin_base_class (Type[PluginBase]): The base class for plugins managed
            by this manager.
        plugins (Dict[str, PluginBase]): Dictionary of loaded plugin instances.
        params (dict): Parameters for configuring the plugin manager and its
            plugins.
    """

    plugin_base_class: Optional[Type[PluginBase]] = None  # Class variable for the base plugin class

    def __init__(self, plugin_base_class: Type[PluginBase]):
        """
        Initialize the plugin manager.

        Args:
            plugin_base_class (Type[PluginBase]): The base class for plugins
                managed by this manager.
        """
        self.plugin_base_class = plugin_base_class
        self.__class__.plugin_base_class = plugin_base_class
        # Load plugins but don't register arguments
        self.plugins = self.load_plugins()

    @classmethod
    def register_plugin_arguments(cls) -> None:
        """Register command-line arguments for all plugins of this type.
        
        This class method discovers and registers arguments for all plugins without 
        requiring full plugin manager instantiation. It should be called during
        argument parsing setup.
        
        Args:
            plugin_base_class: Base class for the type of plugins to discover
        """
        if cls.plugin_base_class is None:
            raise ValueError(f"plugin_base_class not set for {cls.__name__}")
        
        plugin_directory = os.path.dirname(os.path.abspath(__file__))
        # Get the package name by removing the last component of the module path
        package_name = cls.__module__.rsplit(".", 1)[0]
        
        processed_modules = set()  # Keep track of processed modules
        
        for root, dirs, files in os.walk(plugin_directory):
            for file in files:
                if file.endswith(".py") and file != "__init__.py":
                    # Calculate the relative path to build the module name
                    rel_path = os.path.relpath(root, plugin_directory)
                    # Build the full module name
                    if rel_path == ".":
                        module_name = f"{package_name}.{file[:-3]}"
                    else:
                        module_name = f"{package_name}.{rel_path.replace(os.path.sep, '.')}.{file[:-3]}"
                    
                    # Skip if we've already processed this module
                    if module_name in processed_modules:
                        continue
                    processed_modules.add(module_name)
                    
                    try:
                        module = importlib.import_module(module_name)
                        for name, obj in inspect.getmembers(module):
                            if (inspect.isclass(obj) 
                                and issubclass(obj, cls.plugin_base_class)
                                and obj is not cls.plugin_base_class):
                                # Create an instance and register its arguments
                                plugin = obj()
                                title, args = plugin.add_arguments()
                                register_argparse_group(title, args)
                    except ImportError as e:
                        print(f"Failed to import {module_name}: {e}")
                        
    def on_plugin_manager_setup_done(self):
        """
        Perform any initialization tasks for the plugin manager.

        This method is called once during initialization of the plugin manager.
        Children can overwrite it.
        """
        pass

    @property
    def name(self) -> str:
        """
        Get the name of the current class.

        Returns:
            str: The name of the current class.
        """
        return self.__class__.__name__

    def setup_plugin_manager(
        self, params: dict, plugin_managers: Optional["PluginManagerRegistry"] = None,
    ):
        """
        Set up the plugin manager with provided parameters.

        Note: This sets up the plugin manager, not individual plugins.

        Args:
            params (dict): Dictionary of parameters to set up the plugin
                manager.
            plugin_managers (Optional[PluginManagerRegistry]): PluginManagers object
                containing all plugin managers. Used to access plugins from
                within other plugins to avoid circular imports.
        """
        self.params = params

        names_of_plugins_to_load: Optional[List[str]] = None
        plugins_to_load = self.get_selected_plugins_from_params()
        if plugins_to_load is not None:
            names_of_plugins_to_load = [os.path.basename(p) for p in plugins_to_load]

        # Now we filter the plugins to only keep the ones we want to load
        self.plugins = {
            name: plugin
            for name, plugin in self.plugins.items()
            if names_of_plugins_to_load is None or name in names_of_plugins_to_load
        }

        for plugin in self.plugins.values():
            # This also sets params on the plugin.
            plugin._validate(params, plugin_managers)

        self.on_plugin_manager_setup_done()

    def setup_plugins(self, **kwargs):
        """
        Set up the plugins with provided arguments.

        This allows for scenarios where you want to setup a plugin once, then
        execute the run function multiple times.

        Note: This sets up the plugins, not the plugin manager.

        Args:
            **kwargs: Arguments to pass to the plugin setup functions.
        """
        for plugin in self.plugins.values():
            plugin.setup(**kwargs)

    def get_selected_plugins_from_params(self) -> Optional[List[str]]:
        """
        Extract the list of plugins to load from the provided parameters.

        Returns:
            Optional[List[str]]: List of plugins to load, or None to load all
            plugins. Child classes should override this method.
        """
        # For debugging
        print("Available plugins:", list(self.plugins.keys()))
        print("Parameters:", list(self.params.keys()))
        print("Common keys:", set(self.plugins.keys()) & set(self.params.keys()))
        
        keys_in_common = set(self.plugins.keys()) & set(self.params.keys())
        return [key for key in keys_in_common if self.params[key]]

    def load_plugins(self) -> Dict[str, PluginBase]:
        """
        Load all available plugins from the plugin directory.

        Returns:
            Dict[str, PluginBase]: Dictionary mapping plugin names to plugin
                instances.
        """
        plugins: Dict[str, PluginBase] = {}
        plugin_directory = os.path.dirname(os.path.abspath(__file__))
        package_name = self.__class__.__module__.rsplit(".", 1)[0]
        
        for root, dirs, files in os.walk(plugin_directory):
            for file in files:
                if file.endswith(".py") and file != "__init__.py":
                    rel_path = os.path.relpath(root, plugin_directory)
                    if rel_path == ".":
                        module_name = f"{package_name}.{file[:-3]}"
                    else:
                        module_name = f"{package_name}.{rel_path.replace(os.path.sep, '.')}.{file[:-3]}"
                    try:
                        module = importlib.import_module(module_name)
                        for name, obj in inspect.getmembers(module):
                            if (inspect.isclass(obj)
                                and issubclass(obj, self.plugin_base_class)
                                and obj is not self.plugin_base_class):
                                plugin = obj()
                                plugin.on_init()
                                plugins[name] = plugin
                    except ImportError as e:
                        print(f"Failed to import {module_name}: {e}")
        return plugins


    def run(self, cache_dir: Optional[str] = None, **kwargs) -> Any:
        """
        Run the selected plugin(s) with optional caching.

        Note: If your plugin runs only once per generation, you should cache it
        using the cache_dir parameter. If your plugin runs multiple times per
        generation (e.g., mutation or crossover plugins), you should cache it
        manually elsewhere.

        Args:
            cache_dir (Optional[str]): Directory for caching results. If
                provided (same as gen_dir), attempts to use caching.
            **kwargs: Additional arguments to pass to execute().

        Returns:
            Any: The result of executing the plugin(s).
        """
        # Selects which plugin(s) to run and runs them. If cache_dir is provided
        # (same as gen_dir), attempts to use caching.

        # NOTE: If your plugin runs only once per generation, you should cache
        # it using the cache_dir parameter. If your plugin runs multiple times
        # per generation (e.g., the mutation or crossover plugins), you should
        # cache it manually (not via plugin_manager).

        if cache_dir is not None:
            # Using the cache system
            with CacheManager(self.name, cache_dir) as cache:
                if cache.exists:
                    # Cached data exists. Use that instead.
                    return cache.data

                # No cached data, so need to generate
                resp = self.execute(**kwargs)

                # Save to cache
                cache.data = resp

                return resp

        # No cache system
        return self.execute(**kwargs)

    @abstractmethod
    def execute(self, **kwargs) -> Any:
        """
        Execute the selected plugin(s).

        This method should be implemented by child classes to define how
        plugins are selected and executed.

        Args:
            **kwargs: Arguments to pass to the plugin execution.

        Returns:
            Any: The result of executing the plugin(s).
        """
        # Selects which plugin(s) to run and runs them. Defiend on child
        # classes.
        pass

    # @abstractmethod
    # def load_from_cache(self, gen_dir: str) -> Optional[Any]:
    #     """Load the plugin from the cache. All plugins should cache their
    #     results somehow so when autogrow is restarted, it will quickly advance
    #     to the prevous stopping point. Return None if there is no cache."""
    #     pass
