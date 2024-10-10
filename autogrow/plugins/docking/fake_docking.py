import os
import random
from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter
from autogrow.plugins.docking import DockingBase
from autogrow.types import PreDockedCompound


class FakeDocking(DockingBase):
    """
    RUN FAKE DOCKING

    Inputs:
    :param class ParentDocking: Parent docking class to inherit from
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
        return (
            "Fake Docking Options",
            [
                ArgumentVars(
                    name="FakeDocking",
                    action="store_true",
                    default=False,
                    help="Use fake docking (for testing)",
                ),
                #     ArgumentVars(
                #         name="vina_like_executable",
                #         default=None,
                #         help="path to the vina_like_executable (vina, qvina2, smina, etc.)",
                #     ),
                #     ArgumentVars(
                #         name="center_x",
                #         type=float,
                #         default=None,
                #         help="x-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
                #     ),
                #     ArgumentVars(
                #         name="center_y",
                #         type=float,
                #         default=None,
                #         help="y-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
                #     ),
                #     ArgumentVars(
                #         name="center_z",
                #         type=float,
                #         default=None,
                #         help="z-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
                #     ),
                #     ArgumentVars(
                #         name="size_x",
                #         type=float,
                #         default=None,
                #         help="dimension of box to dock into in the x-axis (Angstrom)",
                #     ),
                #     ArgumentVars(
                #         name="size_y",
                #         type=float,
                #         default=None,
                #         help="dimension of box to dock into in the y-axis (Angstrom)",
                #     ),
                #     ArgumentVars(
                #         name="size_z",
                #         type=float,
                #         default=None,
                #         help="dimension of box to dock into in the z-axis (Angstrom)",
                #     ),
                #     ArgumentVars(
                #         name="docking_exhaustiveness",
                #         default=8,
                #         help="exhaustiveness of the global search (roughly proportional to time. \
                #         see docking software for settings.",
                #     ),
                #     ArgumentVars(
                #         name="docking_num_modes",
                #         default=9,
                #         help=" maximum number of binding modes to generate in docking. \
                #         See docking software for settings. ",
                #     ),
                #     ArgumentVars(
                #         name="docking_timeout_limit",
                #         type=float,
                #         default=120,
                #         help="The maximum amount of time allowed to dock a single ligand into a \
                #         pocket in seconds. Many factors influence the time required to dock, such as: \
                #         processor speed, the docking software, rotatable bonds, exhaustiveness docking,\
                #         and number of docking modes... \
                #         The default docking_timeout_limit is 120 seconds, which is excess for most \
                #         docking events using QuickVina2Docking under default settings. If run with \
                #         more exhaustive settings or with highly flexible ligands, consider increasing \
                #         docking_timeout_limit to accommodate. Default docking_timeout_limit is 120 seconds",
                #     ),
            ],
        )

    def validate(self, params: dict):
        """Validate the provided arguments."""
        # TODO: Implement this function
        pass

    def run_docking(self, predocked_cmpd: PreDockedCompound) -> Optional[float]:
        """
        run_docking is needs to be implemented in each class.

        Inputs:
        :param PreDockedCompound predocked_cmpd: A PreDockedCompound object.

        Returns:
        :returns: float score: The score of the docking.
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

        # Successfully docked
        return True, smile_name

    def dock_ligand(self, lig_pdbqt_filename):
        """
        Generate a random docking score for the ligand and create a mock output file.

        Inputs:
        :param str lig_pdbqt_filename: the ligand pdbqt filename
        """
        # Generate a random docking score (lower is better, typically between -12 and 0)
        random_score = round(random.uniform(-12, 0), 2)

        # Create a mock Vina output file
        output_filename = f"{lig_pdbqt_filename}.vina"
        with open(output_filename, "w") as f:
            f.write(f"REMARK VINA RESULT:    {random_score}    0.000    0.000\n")
            f.write("MODEL 1\n")
            f.write("ENDMDL\n")

        print(f"\tGenerated random docking score for: {lig_pdbqt_filename}")
        print(f"\tRandom docking score: {random_score}")
