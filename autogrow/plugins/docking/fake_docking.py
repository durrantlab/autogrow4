import os
import random
from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.docking import DockingBase
from autogrow.types import PostDockedCompound, PreDockedCompound


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
        pass

    def run_docking(
        self, predocked_cmpd: PreDockedCompound
    ) -> Optional[PostDockedCompound]:
        """
        run_docking is needs to be implemented in each class.

        Inputs:
        :param PreDockedCompound predocked_cmpd: A PreDockedCompound object.

        Returns:
        :returns: PostDockedCompound: A PostDockedCompound object, containing
            the score and a docked (posed) SDF file.
        """

        return predocked_cmpd.to_post_docked_compound(
            random.uniform(-12, -8), predocked_cmpd.sdf_3d_path or ""
        )
