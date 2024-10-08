from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompoundInfo
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy
import random

from autogrow.plugins.plugin_base import PluginBase
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


class MutationBase(PluginBase):
    def run(self, **kwargs) -> Optional[List[Union[str, int, None]]]:
        """Run the plugin(s) with provided arguments."""
        return self.run_mutation(kwargs["parent_smiles"])

    @abstractmethod
    def run_mutation(self, parent_smiles: str) -> Optional[List[Union[str, int, None]]]:
        """
        run_mutation is needs to be implemented in each class.

        Inputs:
        :param str parent_smiles: the SMILES string of the parent molecule

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class MutationPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> Optional[List[Union[str, int, None]]]:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        mutation_names = self.get_selected_plugins_from_params()

        if mutation_names is None or len(mutation_names) == 0:
            # TODO: Throw warning? No mutations selected.
            return None

        # Randomly select a mutation to run. If only one is specified, it will
        # always be picked.
        mutation_name = random.choice(mutation_names)

        mutation = cast(MutationBase, self.plugins[mutation_name])

        mutation.run(**kwargs)
        