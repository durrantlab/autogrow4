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

from autogrow.plugins.docking import DockingBase
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.types import Compound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert, obabel_convert_cmd


class VinaLikeDocking(DockingBase):
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
        return [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable docking with Vina-like software (e.g., Vina, QVina2, Smina). This plugin handles file preparation, execution, and result parsing for this family of docking programs.",
            ),
            ArgumentVars(
                name="docking_executable",
                type=str,
                default=None,
                help="path to the docking_executable (vina, qvina2, smina, etc.)",
            ),
            ArgumentVars(
                name="center_x",
                type=float,
                default=None,
                help="x-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
            ),
            ArgumentVars(
                name="center_y",
                type=float,
                default=None,
                help="y-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
            ),
            ArgumentVars(
                name="center_z",
                type=float,
                default=None,
                help="z-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
            ),
            ArgumentVars(
                name="size_x",
                type=float,
                default=None,
                help="dimension of box to dock into in the x-axis (Angstrom)",
            ),
            ArgumentVars(
                name="size_y",
                type=float,
                default=None,
                help="dimension of box to dock into in the y-axis (Angstrom)",
            ),
            ArgumentVars(
                name="size_z",
                type=float,
                default=None,
                help="dimension of box to dock into in the z-axis (Angstrom)",
            ),
            ArgumentVars(
                name="docking_exhaustiveness",
                default=8,
                type=int,
                help="exhaustiveness of the global search (roughly proportional to time. \
     see docking software for settings.",
            ),
            ArgumentVars(
                name="docking_num_modes",
                default=9,
                type=int,
                help="maximum number of binding modes to generate in docking. \
     See docking software for settings. ",
            ),
            ArgumentVars(
                name="docking_nprocs",
                type=int,
                default=1,
                help="number of processors to use for docking (default is 1, \
     which is best when using multithread_mode=multithreading).",
            ),
        ]

    def validate(self, params: dict):
        """
        Validate the provided arguments for Vina-like docking.

        Args:
            params (dict): A dictionary of plugin parameters to validate.

        Raises:
            ValueError: If the required 'obabel_path' is not defined in params.
        """
        if "obabel_path" not in params:
            raise ValueError(
                f"obabel_path must be defined in the params to use {self.name}"
            )
        if "docking_nprocs" in params and params["docking_nprocs"] == -1:
            # If -1, set to number of processors available
            params["docking_nprocs"] = os.cpu_count()

    def run_docking(self, predocked_cmpds: List[Compound]) -> List[Compound]:
        """
        Run docking using Vina-like software on a list of compounds.

        This method prepares input files, executes docking commands, and
        processes docking results for multiple compounds in parallel.

        Args:
            predocked_cmpds (List[Compound]): A list of
                Compound objects to be docked.

        Returns:
            List[Compound]: A list of Compound objects,
            each containing the docking score and the path to the docked
            (posed) SDF file.
        """
        # Convert receptor (PDB format) to PDBQT format if necessary
        receptor_pdbqt = self.get_prepared_receptor()

        # Convert the ligands to PDBQT format as well if necessary
        lig_convert_cmds = []
        lig_dock_cmds = []
        vina_out_files = []
        vina_out_convert_cmds_1 = []
        vina_out_convert_cmds_2 = []

        for predocked_cmpd in predocked_cmpds:
            if predocked_cmpd.sdf_path is None:
                log_warning(f"Skipping {predocked_cmpd.id} because sdf_path is None")
                vina_out_files.append(None)
                continue

            lig_pdbqt_filename = self.get_prepared_ligand_and_add_cmd(predocked_cmpd, lig_convert_cmds)

            # Get commands to dock ligand
            cmd = self.get_dock_cmd(lig_pdbqt_filename, receptor_pdbqt)
            lig_dock_cmds.append(cmd)

            # Get vina output files
            vina_out_file = self.get_output_file(lig_pdbqt_filename)
            vina_out_files.append(vina_out_file)

            # Also get the docked compound as an SDF file. It is important that
            # all hydrogens be added so prolif filters work. It was challenging
            # to get open babel to add all hydrogens in this case. But I found
            # that first saving to a PDB file and then converting to SDF with -p
            # 7.4 worked well.

            docked_pdb_intermediate = f"{vina_out_file}.pdb"
            cmd1 = obabel_convert_cmd(
                vina_out_file, docked_pdb_intermediate, self.params["obabel_path"]
            )
            vina_out_convert_cmds_1.append(cmd1)

            docked_sdf = f"{vina_out_file}.sdf"
            cmd2 = obabel_convert_cmd(
                docked_pdb_intermediate, docked_sdf, self.params["obabel_path"], "-p 7.4"
            )
            vina_out_convert_cmds_2.append(cmd2)

        assert self.plugin_managers is not None, "Plugin managers is None"
        assert (
            self.plugin_managers.ShellParallelizer is not None
        ), "Shell parallelizer is None"

        # Convert the ligands to PDBQT format if needed
        # TODO: Need to specify nprocs?
        if len(lig_convert_cmds) > 0:
            self.plugin_managers.ShellParallelizer.run(
                cmds=lig_convert_cmds
            )

        # Dock the ligands
        # TODO: Need to specify nprocs?
        self.plugin_managers.ShellParallelizer.run(
            cmds=lig_dock_cmds
        )

        # Postprocessing the docking output. This is useful for docking
        # programs whose output for the input ligand is a file containing
        # all their poses. Thus, it is required to select the best pose
        # to perform the following steps.
        self.select_best_pose(vina_out_files)

        # First step to convert the docked ligands to SDF format
        # TODO: Need to specify nprocs?
        self.plugin_managers.ShellParallelizer.run(
            cmds=vina_out_convert_cmds_1
        )

        # Second step to convert the docked ligands to SDF format
        # TODO: Need to specify nprocs?
        self.plugin_managers.ShellParallelizer.run(
            cmds=vina_out_convert_cmds_2
        )

        self.process_results_before_exit(predocked_cmpds, vina_out_files)
        return predocked_cmpds

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

    #######################################
    # DOCK USING VINA                     #
    #######################################
    def get_dock_cmd(self, lig_pdbqt_filename, receptor_pdbqt_file) -> str:
        """
        Generate the docking command for a given ligand using Vina-like software.

        Args:
            lig_pdbqt_filename (str): The filename of the ligand.
            receptor_pdbqt_file (str): The filename of the receptor.

        Returns:
            str: The complete command string to run the docking for the given ligand.
        """
        params = self.params
        output_file = self.get_output_file(lig_pdbqt_filename)

        torun = (
            f'"{params["docking_executable"]}" '
            f'--center_x {params["center_x"]} --center_y {params["center_y"]} --center_z {params["center_z"]} '
            f'--size_x {params["size_x"]} --size_y {params["size_y"]} --size_z {params["size_z"]} '
            f'--receptor "{receptor_pdbqt_file}" '
            f'--ligand "{lig_pdbqt_filename}" '
            f'--out "{output_file}" --cpu {params["docking_nprocs"]}'
        )

        # Add optional user variables additional variable
        if (
            params["docking_exhaustiveness"] is not None
            and params["docking_exhaustiveness"] != "None"
        ) and type(params["docking_exhaustiveness"]) in [int, float]:
            torun = f"{torun} --exhaustiveness " + str(
                int(params["docking_exhaustiveness"])
            )
        if (
            params["docking_num_modes"] is not None
            and params["docking_num_modes"] != "None"
        ) and type(params["docking_num_modes"]) in [int, float]:
            torun = f"{torun} --num_modes " + str(int(params["docking_num_modes"]))

        torun = self.add_additional_args_docking_cmd(torun)

        # Add output line MUST ALWAYS INCLUDE THIS LINE
        torun = f'{torun} >>"{lig_pdbqt_filename}_docking_output.txt"  2>>"{lig_pdbqt_filename}_docking_output.txt"' \
            if os.name != "nt" else f'{torun} > "{lig_pdbqt_filename}_docking_output.txt" 2>&1'

        return torun

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
