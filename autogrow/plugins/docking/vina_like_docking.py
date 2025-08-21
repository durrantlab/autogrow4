"""
This module implements a docking plugin for Vina-like docking software.

The VinaLikeDocking class in this module provides functionality to perform
molecular docking using Vina, QVina2, Smina, or similar software. It handles
the conversion of input files, execution of docking commands, and processing
of docking results.

Classes:
    VinaLikeDocking: A docking class for Vina-like docking software.
"""

import __future__
import os
from autogrow.plugins.docking.vina_docking_base import VinaDockingBase
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.types import Compound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert, obabel_convert_cmd

class VinaLikeDocking(VinaDockingBase):
    """
    A docking class for Vina-like docking software (Vina, QVina2, Smina, etc.).

    This class provides methods to set up and execute docking runs using
    Vina-like software, including file format conversion, parameter setting,
    and result processing.
    """

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments required by the Vina-like docking plugin.

        Returns:
            List[ArgumentVars]: A list of ArgumentVars objects defining the
            command-line arguments for Vina-like docking.
        """
        enabling_arg = [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable docking with Vina-like software (e.g., Vina, QVina2, Smina). This plugin handles file preparation, execution, and result parsing for this family of docking programs.",
            )
        ]
        
        base_args = super().add_arguments()
        
        return enabling_arg + base_args

    def get_prepared_receptor(self) -> str:
        """
        Convert the .pdb file containing the receptor to the .pdbqt format required by Vina software.
        """
        receptor_pdbqt = self.params["receptor_path"] + "qt"
        if not os.path.exists(receptor_pdbqt):
            recep_conversion_success = obabel_convert(
                self.params["receptor_path"],
                receptor_pdbqt,
                self.params["obabel_path"],
                extra_params="-xrp",
            )

            assert (
                recep_conversion_success
            ), f"Failed to convert receptor to PDBQT: {self.params['receptor_path']}"

        return receptor_pdbqt

    def get_prepared_ligand_and_add_cmd(self, predocked_cmpd, lig_convert_cmds) -> str:
        """
        Convert the .sdf file containing the ligand to the .pdbqt format required by Vina software.

        Args:
            predocked_cmpd (Compound): A Compound object representing the ligand.
            lig_convert_cmds (list): A list where the command to run to convert the ligand to the .pdbqt format
                will be added.
        """
        lig_pdbqt_filename = f"{predocked_cmpd.sdf_path}.pdbqt"

        # Get commands to convert ligand to pdbqt
        cmd = obabel_convert_cmd(
            predocked_cmpd.sdf_path, lig_pdbqt_filename, self.params["obabel_path"],
        )
        lig_convert_cmds.append(cmd)
        return lig_pdbqt_filename

    def get_output_file(self, lig_filename) -> str:
        """
        Return the path of the output file for the GNINA docking software.

        Args:
            lig_filename: Path of the input ligand.
        """
        return lig_filename + ".vina"

    def select_best_pose(self, lig_output_files):
        """
        Do not make anything.
        """
        pass

    def process_results_before_exit(self, predocked_cmpds, out_files):
        """
        Select the affinity value for each input compound from the corresponding output file of the
        VINA docking software.

        Args:
            predocked_cmpds (list): List of Compound objects containing the poses.
            out_files (list): List of paths representing each output file of the VINA docking software corresponding
                to each pose.
        """
        # Iterate over output files to process the results
        for predocked_cmpd, vina_out_file in zip(predocked_cmpds, out_files):
            if vina_out_file is None or not os.path.exists(vina_out_file):
                # Throw out ones that failed to dock
                log_warning(f"Failed to dock {predocked_cmpd.id}")
                continue

            # TODO: Does this work for qvina2, smina, etc.?
            with open(vina_out_file, "r") as f:
                # MODEL 1
                # REMARK VINA RESULT:    -9.279      0.000      0.000
                # REMARK INTER + INTRA:         -13.559

                lines = f.readlines()
                score = next(
                    (
                        float(line.split()[3])
                        for line in lines
                        if "REMARK VINA RESULT:" in line
                    ),
                    99999,
                )
                if score == 99999:
                    log_warning(f"Failed to parse docking score from {vina_out_file}")

            # TODO: Update a history in the future.
            predocked_cmpd.fitness_score = score
            predocked_cmpd.docking_score = score
            predocked_cmpd.sdf_path = f"{vina_out_file}.sdf"


    def add_additional_args_docking_cmd(self, torun) -> str:
        """
        Do not make anything.
        """
        return torun

    def replace_atoms_not_handled_by_forcefield(self, lig_pdbqt_filename):
        """
        Replace atoms not handled by the forcefield to prevent docking errors.

        This method replaces problematic atoms (e.g., B, Si) with 'A' in the
        PDBQT file to allow docking to proceed.

        Args:
            lig_pdbqt_filename (str): The filename of the ligand in PDBQT format.

        Returns:
            bool: True if changes were made and the ligand should be redocked,
            False otherwise.
        """
        # TODO: Never used, but might consider in the future...
        # VINA/QuickVINA and MGL have problems with the forcefields for
        # certain atom types To correct this, Autodock Vina suggests replacing
        # the
        atoms_to_replace = [
            "B \n",
            "B\n",
            "Si \n",
            "Si\n",
        ]  # add the \n at the end so we replace the end portion of the line
        printout_of_file = ""
        printout_info = ""
        retry = False
        line_count = 0
        with open(lig_pdbqt_filename, "r") as f:
            for line in f:
                line_count = line_count + 1
                if "HETATM" in line:
                    for x in atoms_to_replace:
                        if x in line:
                            line = line.replace(x, "A \n")
                            retry = True

                            printout_info = f"{printout_info}Changing '{str(x.strip())}' to 'A ' in line: {line_count} of {lig_pdbqt_filename}"
                printout_of_file = printout_of_file + line

        if retry is True:
            print(printout_info)
            with open(lig_pdbqt_filename, "w") as f:
                f.write(printout_of_file)
        else:
            printout_info = "\nCheck the docking message for 'Parse error on'"
            printout_info = (
                printout_info
                + "\n\t This ligand failed to dock. Please check that all "
                + "atoms are covered by the docking forcefield"
            )
            printout_info = (
                printout_info
                + "\n\t Any atoms not covered by the forcefield should be "
                + "added to atoms_to_replace in the function "
                + "replace_atoms_not_handled_by_forcefield"
            )
            printout_info += f"\n\t Verify for this ligand: {lig_pdbqt_filename}\n"
            print(printout_info)
        return retry

    ##########################################
    # Convert the dock outputs to a usable formatted .smi file
    # This is mandatory for all Docking classes but
    # implementation and approach varies by docking and scoring choice
    ##########################################

    # TODO: Might be something here that could be useful.
    # Convert PDB to acceptable PDBQT file format before converting
    # def convert_pdb_to_pdbqt_acceptable_format(self, filename):
    #     """
    #     Make sure a PDB file is properly formatted for conversion to pdbqt

    #     Inputs:
    #     :param str filename: the file path of the pdb file to be converted
    #     """
    #     # read in the file
    #     output_lines = []
    #     with open(filename, "r") as f:
    #         for line in f:
    #             line = line.replace("\n", "")
    #             if line[:5] == "ATOM " or line[:7] == "HETATM ":
    #                 # fix things like atom names with two letters
    #                 first = line[:11]
    #                 middle = (
    #                     line[11:17].upper().strip()
    #                 )  # Need to remove whitespaces on both ends
    #                 last = line[17:]

    #                 middle_firstpart = ""
    #                 middle_lastpart = middle

    #                 for _ in range(len(middle_lastpart)):
    #                     if middle_lastpart[:1].isupper() is not True:
    #                         break  # you reached the first number

    #                     middle_firstpart = middle_firstpart + middle_lastpart[:1]
    #                     middle_lastpart = middle_lastpart[1:]
    #                 # now if there are more than two letters in
    #                 # middle_firstpart, keep just two
    #                 if len(middle_firstpart) > 2:
    #                     middle_lastpart = middle_firstpart[2:] + middle_lastpart
    #                     middle_firstpart = middle_firstpart[:2]

    #                 if middle_firstpart not in ["BR", "ZN", "FE", "MN", "CL", "MG"]:
    #                     # so just keep the first letter for the element part
    #                     # of the atom name
    #                     middle_lastpart = middle_firstpart[1:] + middle_lastpart
    #                     middle_firstpart = middle_firstpart[:1]

    #                 middle = middle_firstpart.rjust(3) + middle_lastpart.ljust(3)

    #                 line = first + middle + last

    #                 # make sure all parts of the molecule belong to the same
    #                 # chain and resid
    #                 line = f"{line[:17]}LIG X 999{line[26:]}"

    #             output_lines.append(line)
    #     with open(filename, "w") as f:

    #         for line in output_lines:
    #             f.write(line + "\n")
