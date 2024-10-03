from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import Any, List, Tuple

from autogrow.config.argparser import ArgumentVars


class PluginBase(ABC):
    def onInit(self):
        # children can overwrite
        pass

    @abstractmethod
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
        pass

    def _validate(self, params: dict):
        """Validate the provided arguments. This is called by the plugin manager."""
        self.params = params
        self.validate(params)

    @abstractmethod
    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass

    @abstractmethod
    def run(self, **kwargs) -> Any:
        """Run the plugin with provided arguments."""
        pass

    @property
    def name(self) -> str:
        """
        Returns the current class name.

        Returns:
        :returns: str self.__class__.__name__: the current class name.
        """
        return self.__class__.__name__
