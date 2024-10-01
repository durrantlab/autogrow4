"""
The child classes from ParentExample
"""
import __future__

import os
import sys
from typing import Any, Dict, Optional

from autogrow.docking.docking_class.docking_class_children.vina_docking import (
    VinaDocking,
)


class QuickVina2Docking(VinaDocking):
    """
    RUN QuickVina2 Docking

    Inputs:
    :param class ParentDocking: Parent docking class to inherit from
    """

    def __init__(
        self,
        params: Optional[Dict[str, Any]] = None,
        receptor_file: Optional[str] = None,
        file_conversion_class_object: Optional[Any] = None,
        test_boot: bool = True,
    ) -> None:
        """
        get the specifications for Vina/QuickVina2 from vrs load them into
        the self variables we will need and convert the receptor to the proper
        file format (ie pdb-> pdbqt)

        Inputs:
        :param dict params: Dictionary of User variables
        :param str receptor_file: the path for the receptor pdb
        :param obj file_conversion_class_object: object which is used to
            convert files from pdb to pdbqt
        :param bool test_boot: used to initialize class without objects for
            testing purpose
        """

        if not test_boot:

            assert params is not None, "params must be passed to QuickVina2Docking"

            self.params = params
            self.debug_mode = params["debug_mode"]
            self.file_conversion_class_object = file_conversion_class_object

            # VINA SPECIFIC VARS
            receptor_file = params["filename_of_receptor"]
            # mgl_python = params["mgl_python"]
            # receptor_template = params["prepare_receptor4.py"]
            # number_of_processors = params["number_of_processors"]
            # docking_executable = params["docking_executable"]

            ###########################

            self.receptor_pdbqt_file = f"{receptor_file}qt"

            self.params["docking_executable"] = self.get_docking_executable_file(
                self.params
            )

    def get_docking_executable_file(self, params: Dict[str, Any]) -> str:
        """
        This retrieves the docking executable files Path.

        Inputs:
        :param dict params: Dictionary of User variables

        Returns:
        :returns: str docking_executable: String for the docking executable
            file path
        """

        # This must already be true if we are here params["dock_choice"] ==
        # "QuickVina2Docking"

        if params["docking_executable"] is None:
            # get default docking_executable for QuickVina2
            script_dir = str(os.path.dirname(os.path.realpath(__file__)))
            docking_executable_directory = (
                (script_dir.split(f"{os.sep}docking_class")[0] + os.sep)
                + "docking_executables"
            ) + os.sep

            if sys.platform in ["linux", "linux2"]:
                # Use linux version of Autodock Vina
                docking_executable = (
                    docking_executable_directory
                    + "q_vina_2"
                    + os.sep
                    + "q_vina_2_1_linux"
                    + os.sep
                    + "qvina2.1"
                )

            elif sys.platform == "darwin":
                # Use OS X version of Autodock Vina
                docking_executable = (
                    docking_executable_directory
                    + "q_vina_2"
                    + os.sep
                    + "q_vina_2_1_mac"
                    + os.sep
                    + "qvina2.1"
                )

            elif sys.platform == "win32":
                # Windows...
                raise Exception("Windows is currently not supported")
            else:
                raise Exception("This OS is currently not supported")

        else:
            # if user specifies a different QuickVina executable
            docking_executable = params["docking_executable"]

        if os.path.exists(docking_executable) is False:
            printout = f"Docking executable could not be found at: {docking_executable}"
            print(printout)
            raise Exception(printout)

        return docking_executable
