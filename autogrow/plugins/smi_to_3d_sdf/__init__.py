"""
Defines base classes and plugin manager for SMILES to 3D SDF conversion in
AutoGrow.

This module provides the SmiTo3DSdfBase abstract base class and
SmiTo3DSdfPluginManager for managing SMILES to 3D SDF conversion plugins.
"""

from abc import ABC, abstractmethod
from argparse import ArgumentParser
import os
from typing import Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound
from autogrow.utils.logging import LogLevel, log_info
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy
import glob

from autogrow.plugins.plugin_base import PluginBase


class SmiTo3DSdfBase(PluginBase):
    """
    An abstract base class for SMILES to 3D SDF converter plugins.

    This class defines the interface for SMILES to 3D SDF converter plugins and
    provides some common utility methods.
    """

    def run(self, **kwargs) -> List[Compound]:
        """
        Run the SMILES to 3D SDF conversion with provided arguments.

        Args:
            **kwargs: Arbitrary keyword arguments. Expected keys:
                - predock_cmpds (List[Compound]): Compounds to convert.
                - pwd (str): Working directory path.

        Returns:
            List[Compound]: Updated list of compounds with 3D SDF
                paths.
        """
        pwd = kwargs["pwd"]
        if pwd[-1] != os.sep:
            pwd += os.sep
        return self.run_smi_to_3d_sdf_converter(kwargs["predock_cmpds"], pwd)

    @abstractmethod
    def run_smi_to_3d_sdf_converter(
        self, predock_cmpds: List[Compound], pwd: str
    ) -> List[Compound]:
        """
        Convert SMILES to 3D SDF files.

        This method must be implemented by subclasses.

        Args:
            predock_cmpds (List[Compound]): List of compounds to
                convert.
            pwd (str): Working directory path.

        Returns:
            List[Compound]: Updated list of compounds with 3D SDF
                paths (must populate sdf_path properties of each
                Compound).
        """
        pass

    def validate(self, params: dict):
        """
        Validate the arguments provided for the SMILES to 3D SDF converter.

        This method is a placeholder and should be overridden by subclasses if
        specific validation is required.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.
        """
        pass


class SmiTo3DSdfPluginManager(PluginManagerBase):
    """
    A plugin manager for SMILES to 3D SDF converter plugins.

    This class manages the selection and execution of SMILES to 3D SDF converter
    plugins.
    """

    def execute(self, **kwargs) -> List[Compound]:
        """
        Execute the selected SMILES to 3D SDF converter plugin.

        This method selects and runs the appropriate SMILES to 3D SDF converter
        based on the provided parameters.

        Args:
            **kwargs: Arbitrary keyword arguments to pass to the selected
                plugin.

        Returns:
            List[Compound]: Updated list of compounds with 3D SDF
                paths.

        Raises:
            Exception: If no SMILES to 3D SDF converter is specified or if
                multiple are selected.

        Note:
            Only one SMILES to 3D SDF converter can be selected at a time.
        """
        smi_to_sdf_converters = self.get_selected_plugins_from_params()

        if smi_to_sdf_converters is None or len(smi_to_sdf_converters) == 0:
            raise Exception(
                f"You must specify an smi-to-3d-sdf Converter! Choose from {str(self.plugins.keys())}"
            )
        if len(smi_to_sdf_converters) > 1:
            raise Exception(
                f"Only one smi-to-3d-sdf Converter can be selected at a time! You selected {smi_to_sdf_converters}"
            )

        # Get the selector plugin to use
        smi_to_sdf_converter = cast(
            SmiTo3DSdfBase, self.plugins[smi_to_sdf_converters[0]]
        )

        # if not glob.glob(kwargs["smi_file"] + "*.sdf"):
        #     raise Exception(
        #         f"Could not find 3D SDF files associated with input file {kwargs['smi_file']}. Conversion error?"
        #     )

        resp = smi_to_sdf_converter.run(**kwargs)

        # Validate that the sdf_path property is populated
        for cmpd in resp:
            if cmpd.sdf_path is None:
                raise Exception(
                    "ERROR! Your SmiTo3DSdf plugin must populate the sdf_path property of each Compound."
                )
            cmpd.add_history("CONVERSION", f"Converted {cmpd.smiles} to 3D SDF file")

        return resp
