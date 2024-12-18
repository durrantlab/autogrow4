from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import TYPE_CHECKING, Any, List, Optional, Tuple

from autogrow.config.argparser import ArgumentVars

if TYPE_CHECKING:
    from autogrow.plugins.plugin_managers import PluginManagers


class PluginBase(ABC):
    """
    Abstract base class for plugins in the AutoGrow system.

    This class defines the basic structure and interface that all plugins should
    follow. It includes methods for initialization, argument handling,
    validation, setup, and execution.

    Attributes:
        params (dict): A dictionary of parameters for the plugin.
        plugin_managers (Optional[PluginManagers]): Reference to the plugin
            managers.
    """

    def on_init(self):
        """
        Initialization method that can be overwritten by child classes.
        """
        # children can overwrite
        pass

    @abstractmethod
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of ArgumentVars objects defining the
                command-line arguments.
        """
        pass

    def _validate(
        self, params: dict, plugin_managers: Optional["PluginManagers"] = None,
    ):
        """
        Validate the provided arguments. This is called by the plugin manager.

        Args:
            params (dict): A dictionary of parameters to validate.
            plugin_managers (Optional[PluginManagers]): Reference to the plugin
                managers.
        """
        self.params = params
        self.plugin_managers = plugin_managers
        self.validate(params)

    def setup(self, **kwargs):
        """
        Setup the plugin with provided arguments.

        This method can be implemented by child classes if needed. It's useful
        for scenarios where you want to setup a plugin only once, then execute
        the run function multiple times.

        Args:
            **kwargs: Arbitrary keyword arguments.

        Returns:
            Any: The result of the setup process, if any.
        """
        return

    @abstractmethod
    def validate(self, params: dict):
        """
        Validate the provided arguments.

        Args:
            params (dict): A dictionary of parameters to validate.
        """
        pass

    @abstractmethod
    def run(self, **kwargs) -> Any:
        """
        Run the plugin with provided arguments.

        Args:
            **kwargs: Arbitrary keyword arguments.

        Returns:
            Any: The result of running the plugin.
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
