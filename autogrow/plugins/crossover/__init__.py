"""Module for crossover operations in autogrow genetic algorithm."""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound
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
                predock_cmpd1 (str): PostDockedCompound of the first ligand.
                predock_cmpd2 (str): PostDockedCompound of the second ligand.

        Returns:
            Optional[str]: The SMILES string of the crossed-over molecule, or None
                if crossover fails.
        """
        return self.run_crossover(kwargs["predock_cmpd1"], kwargs["predock_cmpd2"])

    @abstractmethod
    def run_crossover(
        self, predock_cmpd1: Compound, predock_cmpd2: Compound
    ) -> Optional[str]:
        """
        Implement crossover operation for two ligands.

        This method needs to be implemented in each subclass.

        Args:
            predock_cmpd1 (PostDockedCompound): PostDockedCompound of the first
                ligand.
            predock_cmpd2 (PostDockedCompound):PostDockedCompound of the second
                ligand.

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

        resp = crossover.run(**kwargs)

        if resp in [kwargs["predock_cmpd1"].smiles, kwargs["predock_cmpd2"].smiles]:
            # log_warning("Crossover failed. Child molecule is the same as one of the parents. Skipping.")
            return None

        return resp
