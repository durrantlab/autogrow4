"""
Provides base classes and manager for rescoring plugins in AutoGrow.

This module defines the base classes for rescoring plugins and a plugin manager
for handling multiple rescoring plugins. It includes abstract methods that must
be implemented by specific rescoring plugins.
"""

from abc import abstractmethod
from typing import Optional

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound, Compound

from autogrow.plugins.plugin_base import PluginBase


class RescoringBase(PluginBase):
    """
    Base class for rescoring plugins.

    This abstract class defines the interface that all rescoring plugins must
    implement. It inherits from PluginBase and adds rescoring-specific methods.
    """

    def run(self, **kwargs) -> Optional[Compound]:
        """
        Run the rescoring plugin with provided arguments.

        This method serves as a wrapper for the run_rescoring method, extracting
        the 'mol' from the kwargs and passing it to run_rescoring.

        Args:
            **kwargs: Arbitrary keyword arguments. Must include 'mol'.

        Returns:
            Optional[Compound]: Result of the rescoring operation.
                Returns None if rescoring fails.
        """
        return self.run_rescoring(kwargs["mol"])

    @abstractmethod
    def run_rescoring(self, file_path: str) -> Optional[Compound]:
        """
        Abstract method to be implemented by each rescoring plugin.

        This method should contain the core logic for performing the rescoring
        operation.

        Args:
            file_path (str): The path to the file to be scored.

        Returns:
            Optional[Compound]: Result of the rescoring operation.
                Returns None if rescoring fails.
        """
        # raise NotImplementedError("run_scoring() not implemented")
        pass

    def validate(self, params: dict):
        """
        Validate the provided arguments for the rescoring plugin.

        This method should be implemented to check the validity of the
        parameters passed to the plugin.

        Args:
            params (dict): A dictionary of parameters to validate.

        Note:
            This method is currently a placeholder and needs to be implemented.
        """
        pass


class RescoringPluginManager(PluginManagerBase):
    """
    Manager class for handling multiple rescoring plugins.

    This class is responsible for executing rescoring plugins based on the
    provided arguments and managing the selection of plugins to run.
    """

    def execute(self, **kwargs) -> Optional[Compound]:
        """
        Execute the rescoring plugin with the provided arguments.

        Args:
            **kwargs: A dictionary of arguments to pass to the plugin. Must
                include 'mol'.

        Returns:
            Optional[Compound]: A Compound object containing 
                the rescored information. Returns None if rescoring fails.

        Note:
            This method is currently a placeholder and needs to be implemented.
        """
        # TODO: Implement here.
