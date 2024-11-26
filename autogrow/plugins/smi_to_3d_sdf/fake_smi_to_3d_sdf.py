"""
Implements a fake SMILES to 3D SDF converter for AutoGrow.

This module provides the FakeSmiTo3DSDF class, which simulates the conversion of
SMILES representations to 3D SDF files for testing purposes.
"""
import __future__

from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
import os


class FakeSmiTo3DSDF(SmiTo3DSdfBase):
    """
    A plugin that simulates conversion of SMILES to 3D SDF files.

    This class extends SmiTo3DSdfBase to provide a fake SMILES to 3D SDF
    conversion functionality for testing purposes. It actually creates 2D
    representations for speed and should be used in conjunction with
    FakeDocking.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the Fake SMILES to 3D SDF
        converter.

        This method defines the command-line arguments that can be used to
        configure the Fake SMILES to 3D SDF converter.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("SMILES-to-3D-SDF Converter")
                - A list with one ArgumentVars object defining the argument to
                  enable the Fake SMILES to 3D SDF converter
        """
        return (
            "SMILES-to-3D-SDF Converter",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Fake conversion to 3D SDF. Converts only to 2D for speed. This is for testing purposes only. Use together with FakeDocking.",
                )
            ],
        )

    def validate(self, params: dict):
        """
        Validate the arguments provided for the Fake SMILES to 3D SDF converter.

        This method is a placeholder and currently performs no validation.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.
                Not used in the current implementation.
        """
        pass

    def run_smi_to_3d_sdf_converter(
        self, predock_cmpds: List[PreDockedCompound], pwd: str
    ) -> List[PreDockedCompound]:
        """
        Simulate the conversion of SMILES representations to 3D SDF files.

        This method takes a list of PreDockedCompound objects containing SMILES
        strings and creates fake 3D SDF files. Instead of actual 3D conversion,
        it creates placeholder files with a simple string content.

        Args:
            predock_cmpds (List[PreDockedCompound]): A list of PreDockedCompound
                objects, each containing a SMILES string and other compound
                information.
            pwd (str): The path to the working directory where fake SDF files
                will be created.

        Returns:
            List[PreDockedCompound]: The input list of PreDockedCompound
                objects, updated with the paths to the generated fake 3D SDF
                files.

        Note:
            - This method is intended for testing purposes only and should be
              used in conjunction with FakeDocking.
            - The generated SDF files contain only the string "fake 3D SDF
              file".
        """
        for cmpd_idx, predock_cmpd in enumerate(predock_cmpds):
            out_file = f"{pwd}compound{cmpd_idx}.sdf"
            predock_cmpd.sdf_path = out_file

            with open(out_file, "w") as f:
                f.write("fake 3D SDF file")

        return predock_cmpds
