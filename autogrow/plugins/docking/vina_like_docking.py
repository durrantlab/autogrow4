import __future__
import glob
import os
import sys

from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter
from autogrow.plugins.docking import DockingBase
from typing import Any, Dict, List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.types import PostDockedCompoundInfo, PreDockedCompoundInfo, ScoreType
import autogrow.docking.delete_failed_mol as Delete


class VinaLikeDocking(DockingBase):
    # TODO: I feel like file conversion should not be external to this class...
    file_conversion_class_object: Optional[ParentPDBQTConverter] = None

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
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
                    help="exhaustiveness of the global search (roughly proportional to time. \
                    see docking software for settings.",
                ),
                ArgumentVars(
                    name="docking_num_modes",
                    default=9,
                    help=" maximum number of binding modes to generate in docking. \
                    See docking software for settings. ",
                ),
                ArgumentVars(
                    name="docking_timeout_limit",
                    type=float,
                    default=120,
                    help="The maximum amount of time allowed to dock a single ligand into a \
                    pocket in seconds. Many factors influence the time required to dock, such as: \
                    processor speed, the docking software, rotatable bonds, exhaustiveness docking,\
                    and number of docking modes... \
                    The default docking_timeout_limit is 120 seconds, which is excess for most \
                    docking events using QuickVina2Docking under default settings. If run with \
                    more exhaustive settings or with highly flexible ligands, consider increasing \
                    docking_timeout_limit to accommodate. Default docking_timeout_limit is 120 seconds",
                ),
            ],
        )

        # parser.add_argument(
        #     "--center_x",
        #     "-x",
        #     type=float,
        #     help="x-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
        # )
        # parser.add_argument(
        #     "--center_y",
        #     "-y",
        #     type=float,
        #     help="y-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
        # )
        # parser.add_argument(
        #     "--center_z",
        #     "-z",
        #     type=float,
        #     help="z-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
        # )

        # parser.add_argument(
        #     "--size_x",
        #     type=float,
        #     help="dimension of box to dock into in the x-axis (Angstrom)",
        # )
        # parser.add_argument(
        #     "--size_y",
        #     type=float,
        #     help="dimension of box to dock into in the y-axis (Angstrom)",
        # )
        # parser.add_argument(
        #     "--size_z",
        #     type=float,
        #     help="dimension of box to dock into in the z-axis (Angstrom)",
        # )

        # parser.add_argument(
        #     "--vina_like_executable",
        #     metavar="vina_like_executable",
        #     default=None,
        #     help="path to the vina_like_executable (vina, qvina2, smina, etc.)",
        # )
        # parser.add_argument(
        #     "--docking_exhaustiveness",
        #     metavar="docking_exhaustiveness",
        #     default=None,
        #     help="exhaustiveness of the global search (roughly proportional to time. \
        #     see docking software for settings. Unless specified Autogrow uses the \
        #     docking softwares default setting. For AutoDock Vina 1.1.2 that is 8",
        # )
        # parser.add_argument(
        #     "--docking_num_modes",
        #     metavar="docking_num_modes",
        #     default=None,
        #     help=" maximum number of binding modes to generate in docking. \
        #     See docking software for settings. Unless specified Autogrow uses the \
        #     docking softwares default setting. For AutoDock Vina 1.1.2 that is 9",
        # )
        # parser.add_argument(
        #     "--docking_timeout_limit",
        #     type=float,
        #     default=120,
        #     help="The maximum amount of time allowed to dock a single ligand into a \
        #     pocket in seconds. Many factors influence the time required to dock, such as: \
        #     processor speed, the docking software, rotatable bonds, exhaustiveness docking,\
        #     and number of docking modes... \
        #     The default docking_timeout_limit is 120 seconds, which is excess for most \
        #     docking events using QuickVina2Docking under default settings. If run with \
        #     more exhaustive settings or with highly flexible ligands, consider increasing \
        #     docking_timeout_limit to accommodate. Default docking_timeout_limit is 120 seconds",
        # )

    def validate(self, params: dict):
        """Validate the provided arguments."""
        # TODO: Implement this function
        pass

    def run_docking(
        self, lig_pdbqt_filename, file_conversion_class_object: ParentPDBQTConverter
    ) -> Optional[str]:
        """
        this function runs the docking. Returns None if it worked and the name
        if it failed to dock.

        Inputs:
        :param str pdbqt_filename: the pdbqt file of a ligand to dock and
            score

        Returns:
        :returns: str smile_name: name of smiles if it failed to dock returns
            None if it docked properly
        """

        self.file_conversion_class_object = file_conversion_class_object

        # log("Docking compounds using AutoDock Vina...")
        self.dock_ligand(lig_pdbqt_filename)

        # check that it docked
        pdb_filename = lig_pdbqt_filename.replace("qt", "")

        did_it_dock, smile_name = self.check_docked(pdb_filename)

        if did_it_dock is False:
            # Docking failed

            if smile_name is None:
                print("Missing pdb and pdbqt files for : ", lig_pdbqt_filename)

            return smile_name

        return None

    #######################################
    # DOCK USING VINA                     #
    #######################################
    def dock_ligand(self, lig_pdbqt_filename):
        """
        Dock the ligand pdbqt files in a given directory using AutoDock Vina

        Inputs:
        :param str lig_pdbqt_filename: the ligand pdbqt filename
        """
        params = self.params

        timeout_option = params["timeout_vs_gtimeout"]
        docking_timeout_limit = params["docking_timeout_limit"]
        # do the docking of the ligand Run with a timeout_option limit.
        # Default setting is 5 minutes. This is excessive as most things run
        # within 30seconds This will prevent stalling out. timeout or gtimeout
        receptor_pdbqt_file = f'{params["filename_of_receptor"]}qt'
        torun = (
            f'{timeout_option} {docking_timeout_limit} {params["vina_like_executable"]} '
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

        print(f"\tDocking: {lig_pdbqt_filename}")
        results = self.execute_docking_vina(torun)

        if results is None or results == 256:
            made_changes = self.replace_atoms_not_handled_by_forcefield(
                lig_pdbqt_filename
            )
            if made_changes is True:
                results = self.execute_docking_vina(torun)
                if results == 256 or results is None:
                    print(
                        f"\nLigand failed to dock after corrections: {lig_pdbqt_filename}\n"
                    )
            else:
                print(f"\tFinished Docking: {lig_pdbqt_filename}")
        else:
            print(f"\tFinished Docking: {lig_pdbqt_filename}")

    def replace_atoms_not_handled_by_forcefield(self, lig_pdbqt_filename):
        """
        Replaces atoms not handled by the forcefield to prevent errors. Atoms
        include B and Si.

        Inputs:
        :param str lig_pdbqt_filename: the ligand pdbqt filename

        Returns:
        :returns: bool retry: If True it will be ligand will be redocked, if
            False its dones and wont be docked again.
        """

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

    def execute_docking_vina(self, command):
        """
        Run a single docking execution command

        Inputs:
        :param str command: string of command to run.

        Returns:
        :returns: int result: the exit output for the command. If its None of
            256 it failed.
        """

        try:
            result = os.system(command)
        except Exception:
            result = None
            print(f"Failed to execute: {command}")
        return result

    def check_docked(self, pdb_file):
        """
        given a pdb_file name, test if a pdbqt.vina was created. If it failed
        to dock delete the file pdb and pdbqt file for that ligand -then
        return false

        if it docked properly return True

        Inputs:
        :param str pdb_file: pdb file path

        Returns:
        :returns: bool bool: false if not vina was unsuccessful
        :returns: str smile_name: name of the pdb file
        """

        if not os.path.exists(pdb_file):
            # PDB file doesn't exist
            return False, None
        assert (
            self.file_conversion_class_object is not None
        ), "file_conversion_class_object must be passed to VinaDocking"
        smile_name = self.file_conversion_class_object.get_smile_name_from_pdb(pdb_file)
        if not os.path.exists(f"{pdb_file}qt.vina"):
            # so this pdbqt.vina file didn't exist
            if self.params["debug_mode"] is False:
                print(
                    "Docking unsuccessful: Deleting "
                    + os.path.basename(pdb_file)
                    + "..."
                )

                # REMOVE Failed molecules. Delete ones that were not docked
                # successfully
                Delete.delete_all_associated_files(pdb_file)
                # # delete pdbqt_file
                pdbqt_file = f"{pdb_file}qt"
                Delete.delete_all_associated_files(pdbqt_file)

                return False, smile_name

            # Failed to dock but in debug mode
            print(f"Docking unsuccessful: {os.path.basename(pdb_file)}...")
            return False, smile_name

        # Successfully docked
        return True, smile_name

    ##########################################
    # Convert the dock outputs to a usable formatted .smi file
    # This is mandatory for all Docking classes but
    # implementation and approach varies by docking and scoring choice
    ##########################################
