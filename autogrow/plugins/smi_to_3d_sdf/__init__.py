from abc import ABC, abstractmethod
from argparse import ArgumentParser
import os
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompound
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy
import glob

from autogrow.plugins.plugin_base import PluginBase
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


class SmiTo3DSdfBase(PluginBase):
    def run(self, **kwargs) -> PreDockedCompound:
        """Run the plugin(s) with provided arguments."""
        pwd = kwargs["pwd"]
        if pwd[-1] != os.sep:
            pwd += os.sep
        return self.run_smi_to_3d_sdf_convertor(
            kwargs["predock_cmpd"], pwd, kwargs["cmpd_idx"]
        )

    @abstractmethod
    def run_smi_to_3d_sdf_convertor(
        self, predock_cmpd: PreDockedCompound, pwd: str, cmpd_idx: int
    ) -> PreDockedCompound:
        """
        run_smi_to_sdf_convertor is needs to be implemented in each class.

        Inputs:
        :param str predock_cmpd: A PreDockedCompound object. Conains a
            SMILES string, a name, etc.
        :param str pwd: The path to the working directory.
        :param int cmpd_idx: The index of the compound in the generation.

        Returns:
        :returns: PreDockedCompound: A PreDockedCompound,
            the same as the input, but with the sdf_3d_path field filled in.
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class SmiTo3DSdfPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> PreDockedCompound:
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
            SmiTo3DSdfBase, self.plugins[smi_to_sdf_convertors[0]]
        )

        # if not glob.glob(kwargs["smi_file"] + "*.sdf"):
        #     raise Exception(
        #         f"Could not find 3D SDF files associated with input file {kwargs['smi_file']}. Conversion error?"
        #     )

        return smi_to_sdf_convertor.run(**kwargs)