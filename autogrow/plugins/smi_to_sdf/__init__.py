from abc import ABC, abstractmethod
from argparse import ArgumentParser
import os
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompoundInfo
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy
import glob

from autogrow.plugins.plugin_base import PluginBase
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


class SmiToSdfBase(PluginBase):
    def run(self, **kwargs) -> None:
        """Run the plugin(s) with provided arguments."""
        self.run_smi_to_sdf_convertor(kwargs["smi_file"])

    @abstractmethod
    def run_smi_to_sdf_convertor(self, smi_file: str) -> None:
        """
        run_smi_to_sdf_convertor is needs to be implemented in each class.

        Inputs:
        :param str smi_file: The file path to the SMILES file. Note that it is
            a single with with many smi files in it.

        Returns:
        :returns: str: The file path to the 3D SDF file. This one file contains
            multiple molecules (all those that could be successfully converted).
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class SmiToSdfPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> None:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin
        """

        smi_to_sdf_convertors = self.get_selected_plugins_from_params()

        if smi_to_sdf_convertors is None or len(smi_to_sdf_convertors) == 0:
            raise Exception(
                f"You must specify an smi-to-3d-sdf convertor! Choose from {str(self.plugins.keys())}"
            )
        if len(smi_to_sdf_convertors) > 1:
            raise Exception(
                f"Only one smi-to-3d-sdf convertor can be selected at a time! You selected {smi_to_sdf_convertors}"
            )

        # Get the selector plugin to use
        smi_to_sdf_convertor = cast(
            SmiToSdfBase, self.plugins[smi_to_sdf_convertors[0]]
        )

        smi_to_sdf_convertor.run(**kwargs)

        if not glob.glob(kwargs["smi_file"] + "*.sdf"):
            raise Exception(f"Could not find 3D SDF files associated with input file {kwargs['smi_file']}. Conversion error?")

