"""
This module implements a fake docking plugin for testing purposes.

The FakeDocking class in this module simulates a docking process by assigning
random scores to compounds without performing actual docking calculations. 
This is useful for testing the AutoGrow pipeline without the computational 
overhead of real docking simulations.

Classes:
    FakeDocking: A fake docking class that inherits from DockingBase.
"""

import os
import random
from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.docking import DockingBase
from autogrow.types import PostDockedCompound, PreDockedCompound


class FakeDocking(DockingBase):
    """
    A fake docking class for testing purposes.

    This class implements a simple docking simulation that assigns random
    scores to compounds without actually performing docking calculations.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the fake docking plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
            group name and a list of ArgumentVars objects defining the
            command-line arguments.
        """
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
                #     )
            ],
        )

    def validate(self, params: dict):
        """
        Validate the provided arguments for the fake docking plugin.

        Args:
            params (dict): A dictionary of plugin parameters to validate.
        """
        pass

    def run_docking(
        self, predocked_cmpds: List[PreDockedCompound]
    ) -> List[PostDockedCompound]:
        """
        Perform fake docking on a list of compounds.

        This method assigns random docking scores to the input compounds
        without performing actual docking calculations.

        Args:
            predocked_cmpds (List[PreDockedCompound]): A list of
                PreDockedCompound objects to be "docked".

        Returns:
            List[PostDockedCompound]: A list of PostDockedCompound objects,
            each containing a random score between -12 and -8, and the original
            3D SDF file path.
        """
        return [
            predocked_cmpd.to_post_docked_compound(
                random.uniform(-12, -8), predocked_cmpd.sdf_path or ""
            )
            for predocked_cmpd in predocked_cmpds
        ]
