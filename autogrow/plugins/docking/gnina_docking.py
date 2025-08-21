"""
This module implements a docking plugin for GNINA docking software.

Classes:
    GNINADocking: A docking class for GNINA docking software.
"""

import __future__
import os

from autogrow.plugins.docking.vina_like_docking import VinaLikeDocking
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from rdkit import Chem
from autogrow.utils.logging import log_warning


class GNINADocking(VinaLikeDocking):
    """
    This module implements a docking plugin for GNINA docking software.
    """

    fitness_score = ""
    additional_docking_args = ""

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the GNINA docking plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
            group name and a list of ArgumentVars objects defining the
            command-line arguments for Vina-like docking.
        """
        args = super().add_arguments()[1]
        args.append(ArgumentVars(
            name="gnina_fitness_score",
            type=str,
            default="CNN_VS",
            help="The fitness score to be used to guide the AutoGrow optimization process. The possible scores are:"
                 "CNNscore, CNNaffinity, and CNN_VS.",
        ))
        args.append(ArgumentVars(
            name="gnina_no_gpu",
            action="store_true",
            default=False,
            help="If specified, disable GPU acceleration, even if available.",
        ))
        args.append(ArgumentVars(
            name="gnina_additional_args",
            type=str,
            default=None,
            help="Path to a text file containing additional parameters of the GNINA docking software.",
        ))
        return "GNINA Docking Options", args

    def validate(self, params: dict):
        """
        Validate the provided arguments for GNINA docking.

        Args:
            params (dict): A dictionary of plugin parameters to validate.

        Raises:
            ValueError: If the 'additional_gnina_args' was specified but it does notexist.
        """
        super().validate(params)
        if not params["gnina_fitness_score"] in ["CNNscore", "CNNaffinity", "CNN_VS"]:
            raise ValueError(
                f"The {self.name} is not a valid fitness value to be used for the AUtoGrow optimization. "
                f"See --gnina_fitness_score parameter for more information."
            )
        self.fitness_score = params["gnina_fitness_score"]

        if params["gnina_no_gpu"]:
            self.additional_docking_args = self.additional_docking_args + " --no_gpu"

        if "gnina_additional_args" in params:
            if not os.path.exists(params["gnina_additional_args"]):
                raise ValueError(
                    f"The text file containing additional GNINA arguments does not exist {self.name}"
                )
            with open(params["gnina_additional_args"], "r") as f:
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

            # read the best pose from the gnina output file, and save in a new sdf file
            reader = Chem.SDMolSupplier(lig_new_filename)
            writer = Chem.SDWriter(lig_output_file)
            # It is always ranked as the best one
            if self.fitness_score == "CNN_VS":
                for compound in reader:
                    writer.write(compound)
                    break
            # Search the best compound. The greater the fitness value, the better the compound.
            else:
                best_fitness = 0
                best_compound = None
                for current_compound in reader:
                    fitness_val = float(current_compound.GetProp(self.fitness_score))
                    if fitness_val > best_fitness:
                        best_fitness = fitness_val
                        best_compound = current_compound
                writer.write(best_compound)
            # Closing writer
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
                # Get multiplied by -1 to optimize according to the minimum value as it performs with the docking score
                predocked_cmpd.fitness_score = float(compound.GetProp(self.fitness_score)) * -1
                # Getting the docking score
                predocked_cmpd.docking_score = float(compound.GetProp("minimizedAffinity"))
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
