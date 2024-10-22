"""Module for crossover operations in autogrow genetic algorithm."""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.utils.logging import log_warning
import copy
import random

from autogrow.plugins.plugin_base import PluginBase


class CrossoverBase(PluginBase):
    """
    Abstract base class for crossover operations in the autogrow framework.

    This class defines the interface for crossover plugins. Subclasses must
    implement the `run_crossover` method to define specific crossover behavior.
    """

    def run(self, **kwargs) -> Optional[str]:
        """
        Run the crossover plugin with provided arguments.

        Args:
            **kwargs: Arbitrary keyword arguments including:
                lig_string_1 (str): SMILES string of the first ligand.
                lig_string_2 (str): SMILES string of the second ligand.

        Returns:
            Optional[str]: The SMILES string of the crossed-over molecule, or None
                if crossover fails.
        """
        return self.run_crossover(kwargs["lig_string_1"], kwargs["lig_string_2"])

    @abstractmethod
    def run_crossover(self, lig_string_1: str, lig_string_2: str) -> Optional[str]:
        """
        Implement crossover operation for two ligands.

        This method needs to be implemented in each subclass.

        Args:
            lig_string_1 (str): SMILES string of the first ligand.
            lig_string_2 (str): SMILES string of the second ligand.

        Returns:
            Optional[str]: The SMILES string of the crossed-over molecule, or None
                if crossover fails.
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class CrossoverPluginManager(PluginManagerBase):
    """
    Manages and executes crossover plugins in the autogrow framework.

    This class is responsible for selecting and running crossover operations
    from the available plugins. It randomly chooses a crossover method from
    the selected plugins and executes it.
    """

    def execute(self, **kwargs) -> Optional[str]:
        """
        Execute a randomly selected crossover plugin with provided arguments.

        Args:
            **kwargs: Arbitrary keyword arguments to pass to the plugin.

        Returns:
            Optional[str]: The SMILES string of the crossed-over molecule, or None
                if crossover fails or no crossovers are selected.

        Notes:
            - If no crossovers plugins are selected, a warning is logged.
            - If only one crossover plugin is specified, it will always be 
                selected.
        """
        crossover_names = self.get_selected_plugins_from_params()

        if crossover_names is None or len(crossover_names) == 0:
            log_warning("No crossovers selected.")
            return None

        # Randomly select a mutation to run. If only one is specified, it will
        # always be picked.
        crossover_name = random.choice(crossover_names)

        crossover = cast(CrossoverBase, self.plugins[crossover_name])

        return crossover.run(**kwargs)
