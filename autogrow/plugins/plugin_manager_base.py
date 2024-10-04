import os
import importlib
import inspect
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Type, TypeVar
from autogrow.config.argparser import register_argparse_group
from autogrow.plugins.plugin_base import PluginBase


class PluginManagerBase(ABC):
    def __init__(self, plugin_base_class: Type[PluginBase]):
        self.plugin_base_class = plugin_base_class

        # Initially loads all plugins
        self.plugins = self.load_plugins()

        plugin_manager_name = self.__class__.__name__
        _pluginManagers[plugin_manager_name] = self

    def setup(self, params: dict):
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

    def get_selected_plugins_from_params(self) -> Optional[List[str]]:
        """
        Extract the list of plugins to load from the provided parameters
        (self.params). Children classes should override this method.

        Returns:
        :returns: list of str: the list of plugins to load, or None if the
            program should load all plugins.
        """
        # Get the keys taht self.plugins and self.params have in common.
        keys_in_common = set(self.plugins.keys()) & set(self.params.keys())
        return [key for key in keys_in_common if self.params[key]]

    def load_plugins(self) -> Dict[str, PluginBase]:
        plugins: Dict[str, PluginBase] = {}
        plugin_directory = os.path.dirname(os.path.abspath(__file__))
        # Get the package name
        package_name = self.__class__.__module__.rsplit(".", 1)[0]

        for root, dirs, files in os.walk(plugin_directory):
            for file in files:
                if file.endswith(".py") and file != "__init__.py":
                    # Calculate the relative path from the plugin directory
                    rel_path = os.path.relpath(root, plugin_directory)
                    # Construct the module name relative to the package
                    if rel_path == ".":
                        module_name = f"{package_name}.{file[:-3]}"
                    else:
                        module_name = f"{package_name}.{rel_path.replace(os.path.sep, '.')}.{file[:-3]}"

                    try:
                        module = importlib.import_module(module_name)

                        for name, obj in inspect.getmembers(module):
                            if (
                                inspect.isclass(obj)
                                and issubclass(obj, self.plugin_base_class)
                                and obj is not self.plugin_base_class
                            ):
                                plugins[name] = obj()
                                plugins[name].onInit()

                                title, args = plugins[name].add_arguments()
                                register_argparse_group(title, args)

                    except ImportError as e:
                        print(f"Failed to import {module_name}: {e}")

        return plugins

    @abstractmethod
    def run(self, **kwargs) -> Any:
        # Selects which plugin(s) to run and runs them.
        pass


_pluginManagers: Dict[str, PluginManagerBase] = {}


def get_plugin_manager(plugin_manager_name: str) -> PluginManagerBase:
    return _pluginManagers[plugin_manager_name]


def get_all_plugin_managers() -> Dict[str, PluginManagerBase]:
    return _pluginManagers
