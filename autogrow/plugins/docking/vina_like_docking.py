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

# from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.docking import DockingBase
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.types import Compound, Compound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert, obabel_convert_cmd


class VinaLikeDocking(DockingBase):
    """
    A docking class for Vina-like docking software (Vina, QVina2, Smina, etc.).

    This class provides methods to set up and execute docking runs using
    Vina-like software, including file format conversion, parameter setting,
    and result processing.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the Vina-like docking plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
            group name and a list of ArgumentVars objects defining the
            command-line arguments for Vina-like docking.
        """
        return (
            "Vina-Like Docking Options",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Use docking software from the vina family (vina, qvina2, smina, etc.)",
                ),
                ArgumentVars(
                    name="vina_like_executable",
                    type=str,
                    default=None,
                    help="path to the vina_like_executable (vina, qvina2, smina, etc.)",
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
                    help=" maximum number of binding modes to generate in docking. \
                    See docking software for settings. ",
                ),
            ],
        )

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

        # Convert the ligands to PDBQT format as well
        lig_convert_cmds = []
        lig_dock_cmds = []
        vina_out_files = []
        vina_out_convert_cmds = []

        for predocked_cmpd in predocked_cmpds:
            if predocked_cmpd.sdf_path is None:
                log_warning(f"Skipping {predocked_cmpd.id} because sdf_path is None")
                vina_out_files.append(None)
                continue

            lig_pdbqt_filename = f"{predocked_cmpd.sdf_path}.pdbqt"

            # Get commands to convert ligand to pdbqt
            cmd = obabel_convert_cmd(
                predocked_cmpd.sdf_path, lig_pdbqt_filename, self.params["obabel_path"],
            )
            lig_convert_cmds.append(cmd)

            # Get commands to dock ligand
            cmd = self.get_dock_cmd(lig_pdbqt_filename)
            lig_dock_cmds.append(cmd)

            # Get vina output files
            vina_out_file = f"{lig_pdbqt_filename}.vina"
            vina_out_files.append(vina_out_file)

            # Also get the docked compound as an SDF file
            docked_sdf = f"{vina_out_file}.sdf"
            cmd = obabel_convert_cmd(
                vina_out_file, docked_sdf, self.params["obabel_path"],
            )
            vina_out_convert_cmds.append(cmd)

        # Convert the ligands to PDBQT format
        assert self.plugin_managers is not None, "Plugin managers is None"
        assert (
            self.plugin_managers.ShellParallelizer is not None
        ), "Shell parallelizer is None"
        self.plugin_managers.ShellParallelizer.run(
            cmds=lig_convert_cmds
        )  # TODO: Need to specify nprocs?

        # Dock the ligands
        self.plugin_managers.ShellParallelizer.run(
            cmds=lig_dock_cmds
        )  # TODO: Need to specify nprocs?

        # Convert the docked ligands to SDF format
        self.plugin_managers.ShellParallelizer.run(
            cmds=vina_out_convert_cmds
        )  # TODO: Need to specify nprocs?

        for predocked_cmpd, vina_out_file in zip(predocked_cmpds, vina_out_files):
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
            predocked_cmpd.docking_score = score
            predocked_cmpd.sdf_path = f"{vina_out_file}.sdf"

        return predocked_cmpds

    #######################################
    # DOCK USING VINA                     #
    #######################################
    def get_dock_cmd(self, lig_pdbqt_filename) -> str:
        """
        Generate the docking command for a given ligand using Vina-like software.

        Args:
            lig_pdbqt_filename (str): The filename of the ligand in PDBQT format.

        Returns:
            str: The complete command string to run the docking for the given ligand.
        """
        params = self.params

        receptor_pdbqt_file = f'{params["receptor_path"]}qt'
        torun = (
            f'{params["vina_like_executable"]} '
            f'--center_x {params["center_x"]} --center_y {params["center_y"]} --center_z {params["center_z"]} '
            f'--size_x {params["size_x"]} --size_y {params["size_y"]} --size_z {params["size_z"]} '
            f"--receptor {receptor_pdbqt_file} "
            f"--ligand {lig_pdbqt_filename} "
            f"--out {lig_pdbqt_filename}.vina --cpu 1"
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

        # Add output line MUST ALWAYS INCLUDE THIS LINE
        torun = f"{torun} >>{lig_pdbqt_filename}_docking_output.txt  2>>{lig_pdbqt_filename}_docking_output.txt"

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
