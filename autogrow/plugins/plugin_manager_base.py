import os
import importlib
import inspect
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Type, TypeVar, TYPE_CHECKING
from autogrow.config.argparser import register_argparse_group
from autogrow.plugins.plugin_base import PluginBase
import pickle as pkl

from autogrow.utils.logging import LogLevel, log_warning

if TYPE_CHECKING:
    from autogrow.plugins.plugin_managers import PluginManagers


class PluginManagerBase(ABC):
    def __init__(self, plugin_base_class: Type[PluginBase]):
        self.plugin_base_class = plugin_base_class

        # Initially loads all plugins
        self.plugins = self.load_plugins()

    @property
    def name(self) -> str:
        return self.__class__.__name__

    def setup_plugin_manager(
        self, params: dict, plugin_managers: Optional["PluginManagers"] = None,
    ):
        """
        Sets up the plugin manager with the provided parameters. NOTE: This
        sets up the plugin manager, not individual plugins.
        
        Inputs:
        :param dict params: a dictionary of parameters to set up the plugin
            manager.
        :param PluginManagers plugin_managers: a PluginManagers object that
            contains all the plugin managers. This is used to access a plugin
            from within another plugin. Simple import would be circular.
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

    def setup_plugins(self, **kwargs):
        """
        Sets up the plugins with the provided arguments. This is required because
        one can imagine a scenario where you want to setup a plugin only once,
        then execute the run function multiple times.
        
        NOTE: This sets up the plugins, not the plugin manager.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin
            setup functions.
        """
        for plugin in self.plugins.values():
            plugin.setup(**kwargs)

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

    def run(self, cache_dir: Optional[str] = None, **kwargs) -> Any:
        # Selects which plugin(s) to run and runs them. If cache_dir is provided
        # (same as gen_dir), attempts to use caching.

        cache_filename = ""
        if cache_dir is not None:
            cache_filename = os.path.join(cache_dir, f"{self.name}_results_cache.pkl")

        if cache_dir is not None and os.path.exists(cache_filename):
            log_warning(f"Loading previous {self.name} results from cache: {cache_filename}")
            return pkl.load(open(cache_filename, "rb"))
        
        resp = self.execute(**kwargs)

        if cache_dir is not None:
            pkl.dump(resp, open(cache_filename, "wb"))
        
        return resp

    @abstractmethod
    def execute(self, **kwargs) -> Any:
        # Selects which plugin(s) to run and runs them. Defiend on child
        # classes.
        pass


    # @abstractmethod
    # def load_from_cache(self, gen_dir: str) -> Optional[Any]:
    #     """Load the plugin from the cache. All plugins should cache their
    #     results somehow so when autogrow is restarted, it will quickly advance
    #     to the prevous stopping point. Return None if there is no cache."""
    #     pass
