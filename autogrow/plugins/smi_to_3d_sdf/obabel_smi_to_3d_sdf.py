"""
Implements an OpenBabel-based SMILES to 3D SDF converter for AutoGrow.

This module provides the ObabelSmiTo3DSDF class, which uses OpenBabel to convert
SMILES representations of compounds to 3D SDF files.
"""
import __future__

from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase
from autogrow.types import Compound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert_cmd
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
import os


class ObabelSmiTo3DSDF(SmiTo3DSdfBase):
    """
    A plugin that uses OpenBabel to convert SMILES to 3D SDF files.

    This class extends SmiTo3DSdfBase to provide SMILES to 3D SDF conversion
    functionality using OpenBabel.
    """

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line args specific to the Obabel SMILES to 3D SDF converter.

        This method defines the command-line arguments that can be used to
        configure the OpenBabel SMILES to 3D SDF converter.

        Returns:
            List[ArgumentVars]: A list with one ArgumentVars object defining
            the argument to enable the OpenBabel SMILES to 3D SDF converter.
        """
        return [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable SMILES to 3D SDF conversion using OpenBabel. This plugin generates 3D coordinates for molecules from their SMILES strings and saves them in SDF format.",
            ),
            ArgumentVars(
                name="ph",
                type=float,
                default=7.4,
                help="pH to use for protonation when generating 3D structures with OpenBabel. Default is 7.4.",
            ),
        ]

    def validate(self, params: dict):
        """
        Validate the args provided for the OpenBabel SMILES to 3D SDF converter.
        This method checks if the required 'obabel_path' parameter is provided and
        if the pH value is valid.
        Args:
            params (dict): A dictionary of parameters provided to the plugin.
        Raises:
            Exception: If the 'obabel_path' parameter is not provided.
            ValueError: If 'ph' is not a valid float.
        """
        if "obabel_path" not in params:
            raise Exception(
                "You must provide the path to obabel via the --obabel_path parameter."
            )
        try:
            ph_val = float(params["ph"])
            if not (0.0 <= ph_val <= 14.0):
                log_warning(
                    f"ph ({ph_val}) is outside the typical pH range of 0-14. This may be unintended."
                )
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid value for ph: {params['ph']}. Must be a number."
            ) from e

    def run_smi_to_3d_sdf_converter(
        self, predock_cmpds: List[Compound], pwd: str
    ) -> List[Compound]:
        """
        Convert a list of SMILES to 3D SDF files using OpenBabel.

        This method takes a list of Compound objects containing SMILES
        strings and converts them to 3D SDF files using OpenBabel. The
        conversion is done in parallel using a shell parallelizer plugin.

        Args:
            predock_cmpds (List[Compound]): A list of Compound
                objects, each containing a SMILES string and other compound
                information.
            pwd (str): The path to the working directory where temporary files
                will be created.
        Returns:
            List[Compound]: The input list of Compound
                objects, updated with the paths to the generated 3D SDF files.

        Note:
            - This method uses OpenBabel with the '--gen3d' and '--p' options.
            - If the conversion fails for a compound, a warning is logged and
                the compound's sdf_path is not set.
            - The method assumes that a ShellParallelizer plugin is available
                for parallel execution of OpenBabel commands.
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). If No matches are found to filter list this will
        # return a True as it Passed the filter.

        cmds = []
        out_files = []
        for cmpd_idx, cmpd in enumerate(predock_cmpds):
            base_file = f"{pwd}compound{cmpd_idx}"
            in_file = f"{base_file}.smi"
            out_file = f"{base_file}.sdf"
            obabel_path = self.params["obabel_path"]
            ph_value = self.params["ph"]
            with open(in_file, "w") as f:
                f.write(cmpd.smiles)
            cmd = obabel_convert_cmd(
                in_file, out_file, obabel_path, extra_params=f"--gen3d --p {ph_value}"
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

        for cmpd_idx, cmpd in enumerate(predock_cmpds):
            out_file = out_files[cmpd_idx]

            if not os.path.exists(out_file):
                log_warning(
                    f"Could not convert smiles to 3D SDF with obabel: {cmpd.smiles}"
                )
                continue

            with open(out_file, "r") as f:
                content = f.read().strip()
                if content == "":
                    log_warning(
                        f"Could not convert smiles to 3D SDF with obabel: {cmpd.smiles}"
                    )
                    continue

            cmpd.sdf_path = out_file

        return predock_cmpds
