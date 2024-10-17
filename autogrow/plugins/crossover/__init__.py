from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.utils.logging import log_warning
import copy
import random

from autogrow.plugins.plugin_base import PluginBase


class CrossoverBase(PluginBase):
    def run(self, **kwargs) -> Optional[str]:
        """Run the plugin(s) with provided arguments."""
        return self.run_crossover(kwargs["lig_string_1"], kwargs["lig_string_2"])

    @abstractmethod
    def run_crossover(self, lig_string_1: str, lig_string_2: str) -> Optional[str]:
        """
        run_mutation is needs to be implemented in each class.

        Inputs:
        :param str parent_smiles: the SMILES string of the parent molecule
        :param list existing_smiles: a list of existing molecules, to prevent
            duplicates.

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        # TODO: Implement validation
        pass


class CrossoverPluginManager(PluginManagerBase):
    def execute(self, **kwargs) -> Optional[str]:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
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
