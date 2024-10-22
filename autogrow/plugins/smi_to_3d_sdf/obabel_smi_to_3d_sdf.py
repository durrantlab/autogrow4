"""
Implements an OpenBabel-based SMILES to 3D SDF converter for AutoGrow.

This module provides the ObabelSmiTo3DSDF class, which uses OpenBabel to convert
SMILES representations of compounds to 3D SDF files.
"""
import __future__

# from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert, obabel_convert_cmd
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
import os


class ObabelSmiTo3DSDF(SmiTo3DSdfBase):
    """
    A plugin that uses OpenBabel to convert SMILES to 3D SDF files.

    This class extends SmiTo3DSdfBase to provide SMILES to 3D SDF conversion
    functionality using OpenBabel.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the OpenBabel SMILES to 3D SDF
        converter.

        This method defines the command-line arguments that can be used to
        configure the OpenBabel SMILES to 3D SDF converter.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("SMILES-to-3D-SDF Converter")
                - A list with one ArgumentVars object defining the argument to
                  enable the OpenBabel SMILES to 3D SDF converter
        """
        return (
            "SMILES-to-3D-SDF Converter",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Convert SMILES to 3D SDF using obabel.",
                )
            ],
        )

    def validate(self, params: dict):
        """
        Validate the arguments provided for the OpenBabel SMILES to 3D SDF
        converter.

        This method checks if the required 'obabel_path' parameter is provided.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.

        Raises:
            Exception: If the 'obabel_path' parameter is not provided.
        """
        if "obabel_path" not in params:
            raise Exception(
                "You must provide the path to obabel via the --obabel_path parameter."
            )

    def run_smi_to_3d_sdf_converter(
        self, predock_cmpds: List[PreDockedCompound], pwd: str
    ) -> List[PreDockedCompound]:
        """
        Convert a list of SMILES representations to 3D SDF files using
        OpenBabel.

        This method takes a list of PreDockedCompound objects containing SMILES
        strings and converts them to 3D SDF files using OpenBabel. The
        conversion is done in parallel using a shell parallelizer plugin.

        Args:
            predock_cmpds (List[PreDockedCompound]): A list of PreDockedCompound
                objects, each containing a SMILES string and other compound
                information.
            pwd (str): The path to the working directory where temporary files
            will be created.

        Returns:
            List[PreDockedCompound]: The input list of PreDockedCompound
                objects, updated with the paths to the generated 3D SDF files.

        Note:
            - This method uses OpenBabel with the '--gen3d' and '--p 7.4'
              options.
            - If the conversion fails for a compound, a warning is logged and
              the compound's sdf_3d_path is not set.
            - The method assumes that a ShellParallelizer plugin is available
              for parallel execution of OpenBabel commands.
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). If No matches are found to filter list this will
        # return a True as it Passed the filter.

        cmds = []
        out_files = []
        for cmpd_idx, predock_cmpd in enumerate(predock_cmpds):
            base_file = f"{pwd}compound{cmpd_idx}"
            in_file = f"{base_file}.smi"
            out_file = f"{base_file}.sdf"
            obabel_path = self.params["obabel_path"]

            with open(in_file, "w") as f:
                f.write(predock_cmpd.smiles)

            cmd = obabel_convert_cmd(
                in_file, out_file, obabel_path, extra_params="--gen3d --p 7.4"
            )

            cmds.append(cmd)
            out_files.append(out_file)

        # Get parallelizer plugin to use
        # TODO: Need to specify nprocs?
        assert self.plugin_managers is not None, "Plugin managers is None"
        assert (
            self.plugin_managers.ShellParallelizer is not None
        ), "Shell parallelizer is None"
        self.plugin_managers.ShellParallelizer.run(cmds=cmds)

        for cmpd_idx, predock_cmpd in enumerate(predock_cmpds):
            out_file = out_files[cmpd_idx]

            if not os.path.exists(out_file):
                log_warning(
                    f"Could not convert smiles to 3D SDF with obabel: {predock_cmpd.smiles}"
                )
                continue

            with open(out_file, "r") as f:
                content = f.read().strip()
                if content == "":
                    log_warning(
                        f"Could not convert smiles to 3D SDF with obabel: {predock_cmpd.smiles}"
                    )
                    continue

            predock_cmpd.sdf_3d_path = out_file

        return predock_cmpds
