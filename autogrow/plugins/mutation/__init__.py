from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import Dict, List, Optional, Tuple, Union

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompoundInfo
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy

from autogrow.plugins.plugin_base import PluginBase
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


class MutationBase(PluginBase):
    def run(self, **kwargs) -> bool:
        """Run the plugin(s) with provided arguments."""
        return self.run_mutation(kwargs["mol"])

    @abstractmethod
    def run_mutation(self, mol: Chem.rdchem.Mol) -> bool:
        """
        run_filter is needs to be implemented in each class.

        Inputs:
        :param rdkit.Chem.rdchem.Mol mol: a molecule to filter

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class MutationPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> List:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        # TODO: Implement
        pass
