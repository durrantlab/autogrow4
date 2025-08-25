"""
This module implements a base docking plugin for Vina-like docking software.
The VinaDockingBase class in this module provides functionality to perform
molecular docking using Vina, QVina2, Smina, or similar software. It handles
the conversion of input files, execution of docking commands, and processing
of docking results.
Classes:
 VinaDockingBase: An abstract base docking class for Vina-like docking software.
"""
import __future__
import os
from abc import abstractmethod, ABC
from autogrow.plugins.docking import DockingBase
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.types import Compound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert, obabel_convert_cmd


class VinaDockingBase(DockingBase, ABC):
    """
    An abstract base class for Vina-like docking software.
    """

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments for Vina-like docking.
        Returns:
            List[ArgumentVars]: A list of ArgumentVars objects defining the
            command-line arguments for Vina-like docking.
        """
        # This is a base class so it does not have its own enabling argument
        return [
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
        receptor_prepared = self.get_prepared_receptor()
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
            lig_prepared_filename = self.get_prepared_ligand_and_add_cmd(predocked_cmpd, lig_convert_cmds)
            # Get commands to dock ligand
            cmd = self.get_dock_cmd(lig_prepared_filename, receptor_prepared)
            lig_dock_cmds.append(cmd)
            # Get vina output files
            vina_out_file = self.get_output_file(lig_prepared_filename)
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

    def get_dock_cmd(self, lig_filename: str, receptor_file: str) -> str:
        """
        Generate the docking command for a given ligand using Vina-like software.
        Args:
            lig_filename (str): The filename of the ligand.
            receptor_file (str): The filename of the receptor.
        Returns:
            str: The complete command string to run the docking for the given ligand.
        """
        params = self.params
        output_file = self.get_output_file(lig_filename)
        torun = (
            f'"{params["docking_executable"]}" '
            f'--center_x {params["center_x"]} --center_y {params["center_y"]} --center_z {params["center_z"]} '
            f'--size_x {params["size_x"]} --size_y {params["size_y"]} --size_z {params["size_z"]} '
            f'--receptor "{receptor_file}" '
            f'--ligand "{lig_filename}" '
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
        torun = f'{torun} >>"{lig_filename}_docking_output.txt"  2>>"{lig_filename}_docking_output.txt"' \
            if os.name != "nt" else f'{torun} > "{lig_filename}_docking_output.txt" 2>&1'
        return torun

    @abstractmethod
    def get_prepared_receptor(self) -> str:
        pass

    @abstractmethod
    def get_prepared_ligand_and_add_cmd(self, predocked_cmpd, lig_convert_cmds) -> str:
        pass

    @abstractmethod
    def get_output_file(self, lig_filename) -> str:
        pass

    @abstractmethod
    def select_best_pose(self, lig_output_files):
        pass

    @abstractmethod
    def process_results_before_exit(self, predocked_cmpds, out_files):
        pass
    
    @abstractmethod
    def add_additional_args_docking_cmd(self, torun) -> str:
        pass