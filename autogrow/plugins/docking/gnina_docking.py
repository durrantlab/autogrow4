"""
This module implements a docking plugin for GNINA docking software.

Classes:
    GNINADocking: A docking class for GNINA docking software.
"""

import __future__
import os
from autogrow.plugins.docking.vina_docking_base import VinaDockingBase
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from rdkit import Chem
from autogrow.utils.logging import log_warning

class GNINADocking(VinaDockingBase):
    """
    This module implements a docking plugin for GNINA docking software.
    """

    additional_docking_args = ""

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments required by the GNINA docking plugin.

        Returns:
            List[ArgumentVars]: A list of ArgumentVars objects defining the
            command-line arguments for Vina-like docking.
        """
        enabling_arg = [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable docking with GNINA software. This plugin handles file preparation, execution, and result parsing for GNINA.",
            )
        ]
        
        gnina_specific_args = [
            ArgumentVars(
                name="additional_gnina_args",
                type=str,
                default=None,
                help="Path to a text file containing additional parameters of the GNINA docking software.",
            )
        ]
        
        base_args = super().add_arguments()
        
        return enabling_arg + gnina_specific_args + base_args

    def validate(self, params: dict):
        """
        Validate the provided arguments for GNINA docking.

        Args:
            params (dict): A dictionary of plugin parameters to validate.

        Raises:
            ValueError: If the 'additional_gnina_args' was specified but it does notexist.
        """
        super().validate(params)
        if params["additional_gnina_args"]:
            if not os.path.exists(params["additional_gnina_args"]):
                raise ValueError(
                    f"The text file containing additional GNINA arguments does not exist {self.name}"
                )
            with open(params["additional_gnina_args"], "r") as f:
                lines = f.readlines()
                for line in lines:
                    if not line.startswith("#"):
                        line = line.replace("\n", "")
                        self.additional_docking_args = self.additional_docking_args + " " + line

    def get_prepared_receptor(self) -> str:
        """
        Return the .pdb file containing the receptor.
        """
        return self.params["receptor_path"]

    def get_prepared_ligand_and_add_cmd(self, predocked_cmpd, lig_convert_cmds) -> str:
        """
        Return the .sdf file containing the ligand.
        """
        return predocked_cmpd.sdf_path

    def get_output_file(self, lig_filename) -> str:
        """
        Return the path of the output file for the GNINA docking software.

        Args:
            lig_filename: Path of the input ligand.
        """
        last_position = lig_filename.rfind(".sdf")
        lig_filename = lig_filename[:last_position] + "_gnina.sdf"
        return lig_filename

    def select_best_pose(self, lig_output_files):
        """
        Select the best pose from the output file of the GNINA docking software.

        Args:
            lig_output_files: Path of the input ligand.
        """
        for lig_output_file in lig_output_files:
            # Rename output file name to finishing with _gnina_all.sdf
            last_position = lig_output_file.rfind(".sdf")
            lig_new_filename = lig_output_file[:last_position] + "_gnina_all.sdf"
            os.rename(lig_output_file, lig_new_filename)

            # read the best pose from gnina output file and save in a new sdf file
            reader = Chem.SDMolSupplier(lig_new_filename)
            writer = Chem.SDWriter(lig_output_file)
            for compound in reader:
                writer.write(compound)
                break
            writer.close()

    def process_results_before_exit(self, predocked_cmpds, out_files):
        """
        Select the affinity value for each input compound from the corresponding output of the GNINA docking software.

        Args:
            predocked_cmpds (list): List of Compound objects containing the poses.
            out_files (list): List of paths representing each output file of the GNINA docking software corresponding
                to each pose.
        """
        # Iterate over output files to process the results
        for predocked_cmpd, out_file in zip(predocked_cmpds, out_files):
            if out_file is None or not os.path.exists(out_file):
                # Throw out ones that failed to dock
                log_warning(f"Failed to dock {predocked_cmpd.id}")
                continue

            reader = Chem.SDMolSupplier(out_file)
            for compound in reader:
                affinity = float(compound.GetProp("minimizedAffinity"))
                predocked_cmpd.fitness_score = affinity
                predocked_cmpd.docking_score = affinity
                predocked_cmpd.sdf_path = out_file
                break

    def add_additional_args_docking_cmd(self, torun) -> str:
        """
        Add to the GNINA command line the additional arguments read from 'additional_gnina_args'. Return
        the new command line.

        Args:
            torun (str): GNINA command line
        """
        return torun + self.additional_docking_args
