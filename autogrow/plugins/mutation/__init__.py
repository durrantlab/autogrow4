"""
Provides base classes and manager for mutation plugins in AutoGrow.

This module defines the base classes for mutation plugins and a plugin manager
for handling multiple mutation plugins. It includes abstract methods that must
be implemented by specific mutation plugins.
"""
from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound
from autogrow.utils.logging import LogLevel, log_debug, log_warning
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy
import random

from autogrow.plugins.plugin_base import PluginBase


class MutationBase(PluginBase):
    """
    Base class for mutation plugins.

    This abstract class defines the interface that all mutation plugins must
    implement. It inherits from PluginBase and adds mutation-specific methods.
    """

    def run(self, **kwargs) -> Optional[Tuple[str, int, Union[str, None]]]:
        """
        Run the mutation plugin with provided arguments.

        Args:
            **kwargs: Arbitrary keyword arguments. Must include 'predock_cmpd'.

        Returns:
            Optional[Tuple[str, int, Union[str, None]]]: A tuple containing:
                - [0] str: SMILES string of the mutated molecule
                - [1] int: ID number of the reaction used
                - [2] Optional[str]: Name of the complementary molecule used (if
                                     any), or None for single-reactant reactions
            Returns None if the mutation fails.
        """
        return self.run_mutation(kwargs["predock_cmpd"])

    @abstractmethod
    def run_mutation(
        self, predock_cmpd: Compound
    ) -> Optional[Tuple[str, int, Union[str, None]]]:
        """
        Abstract method to be implemented by each mutation plugin.

        Args:
            predock_cmpd (PostDockedCompound): The PostDockedCompound of the
                parent molecule to be mutated.

        Returns:
            Optional[Tuple[str, int, Union[str, None]]]: A tuple containing:
                - [0] str: SMILES string of the mutated molecule
                - [1] int: ID number of the reaction used
                - [2] Optional[str]: Name of the complementary molecule used (if
                                     any), or None for single-reactant reactions
            Returns None if the mutation fails.

        Note:
            TODO: The 'existing_smiles' parameter mentioned in the original
            docstring is not present in the method signature. Subclasses should
            handle duplicate prevention if necessary.
        """
        pass


class MutationPluginManager(PluginManagerBase):
    """Manager class for handling multiple mutation plugins."""

    def execute(self, **kwargs) -> Optional[Tuple[str, int, Union[str, None]]]:
        """
        Execute a randomly selected mutation plugin with the provided arguments.

        Args:
            **kwargs: A dictionary of arguments to pass to the plugin. Must
                include 'predock_cmpd'.

        Returns:
            Optional[Tuple[str, int, Union[str, None]]]: A tuple containing:
                - [0] str: SMILES string of the mutated molecule
                - [1] int: ID number of the reaction used
                - [2] Optional[str]: Name of the complementary molecule used (if
                                     any), or None for single-reactant reactions
            Returns None if no mutations are selected or if the mutation fails.

        Note:
            The method logs a warning if no mutations are selected and logs
            debug information about the mutation performed.
        """
        mutation_names = self.get_selected_plugins_from_params()

        if mutation_names is None or len(mutation_names) == 0:
            # Throw warning. No mutations selected.
            log_warning("No mutations selected.")
            return None

        # Randomly select a mutation to run. If only one is specified, it will
        # always be picked.
        mutation_name = random.choice(mutation_names)

        mutation = cast(MutationBase, self.plugins[mutation_name])

        resp = mutation.run(**kwargs)

        if resp is not None:
            log_debug(f'{mutation.name}: {kwargs["predock_cmpd"].smiles} => {resp[0]}')

        return resp
