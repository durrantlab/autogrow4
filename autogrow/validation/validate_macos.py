import platform
import os
import sys
from typing import Any, Dict


def run_macos_notarization(params: Dict[str, Any]) -> None:
    """
    Some MacOS require docking software to be notarized. This will require
    an internet signal. This function runs notarization on vina and qvina2
    docking. This is important for MacOS newer than 10.15 and newer than

    For MacOS newer than 10.15, this will require an internet connection.

    Inputs:
    :param dict params: dict of user variables which will govern how the programs runs
    """
    if sys.platform.lower() != "darwin" or params["dock_choice"] not in [
        "VinaDocking",
        "QuickVina2Docking",
    ]:
        return

    current_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
    vina_exe = current_dir + os.sep.join(
        [
            "..",
            "docking",
            "docking_executables",
            "vina",
            "autodock_vina_1_1_2_mac",
            "bin",
            "vina",
        ]
    )
    qvina2_exe = current_dir + os.sep.join(
        [
            "..",
            "docking",
            "docking_executables",
            "q_vina_2",
            "q_vina_2_1_mac",
            "qvina2.1",
        ]
    )

    # Check executables exist. raise exception if not
    if os.path.exists(vina_exe) is False or os.path.exists(qvina2_exe) is False:
        printout = "Docking executables could not be found."
        raise Exception(printout)

    both_docking_exe_work = _validate_docking_executables(params, vina_exe, qvina2_exe)

    if both_docking_exe_work is False:
        # Ensure permissions are unrestricted
        try:
            _run_cmd_on_docking_execs("chmod -R a+rwx ", vina_exe, qvina2_exe)
        except Exception as e:
            printout = "Permissions could not be adjusted on docking executable files."
            print(printout)
            raise Exception(printout) from e

        # Check Platform information
        mac_version = platform.mac_ver()[0].split(".")
        if int(mac_version[0]) < 10:
            _throw_error(
                "We do not provide support for MacOS less than 10.7.\n",
                "Please run using docker version of AutoGrow.",
            )
        if int(mac_version[0]) == 10:
            if int(mac_version[1]) < 7:
                _throw_error(
                    "We do not support for MacOS less than 10.7.\n",
                    "Please run using docker version of AutoGrow.",
                )

            # TODO: Revisit this later. Might just want to require paths rather
            # than packaging the binaries.

            # if int(mac_version[1]) > 15:
            #     # 10.15 is Catalina which requires notarizing docking software

            #     printout = (
            #         "We have not tested MacOS higher than 10.15.\n"
            #         + "Please run using docker version of AutoGrow."
            #     )
            #     print(printout)
            #     raise Exception(printout)

            try:
                _run_cmd_on_docking_execs(
                    "xattr -w com.apple.quarantine ", vina_exe, qvina2_exe
                )
            except Exception:
                _throw_error(
                    "Please install xattr. Can be installed using the command:",
                    "\n\tpip install xattr",
                )


def _run_cmd_on_docking_execs(cmd: str, vina_exe: str, qvina2_exe: str) -> None:
    """
    Executes a given command on both Vina and QuickVina2 executables.

    This function runs the specified command on both the Vina and QuickVina2 
    executables using os.system(). It's typically used for operations like 
    changing permissions or setting attributes.

    Inputs:
    :param str cmd: The command to be executed on the executables
    :param str vina_exe: The path to the Vina executable
    :param str qvina2_exe: The path to the QuickVina2 executable

    Returns:
    :returns: None

    Note:
    This function uses os.system() which can be a security risk if the 
    input is not properly sanitized. Ensure that 'cmd' is safe to execute.
    """
    command = f"{cmd}{vina_exe}"
    os.system(command)
    command = f"{cmd}{qvina2_exe}"
    os.system(command)


def _throw_error(msg1: str, msg2: str) -> None:
    """
    Prints an error message and raises an exception.

    This function concatenates two message strings, prints the resulting message,
    and then raises an Exception with the same message. It's used to handle
    error conditions in a consistent manner.

    Inputs:
    :param str msg1: The first part of the error message
    :param str msg2: The second part of the error message

    Returns:
    :returns: None

    Raises:
    :raises Exception: Always raises an exception with the concatenated error message

    Note:
    Although this function has a return type of None, it will never actually
    return due to the raised exception. The None return type is used for
    compatibility with type checking tools.
    """
    printout = msg1 + msg2
    print(printout)
    raise Exception(printout)


def _validate_docking_executables(
    params: Dict[str, Any], vina_exe: str, qvina2_exe: str
) -> bool:
    """
    This will test if docking executables are compatible with OS.
    This is only required for MacOS.

    Test will output the version of Vina and QVina2.1 executables to txt file
    in the root_output_folder (docking_exe_MACOS_test.txt)
    If both executables are compatible with this MacOS there should be the following
    2 lines in the txt file:
        AutoDock Vina 1.1.2 (May 11, 2011)
        QuickVina 2.1 (24 Dec, 2017)

    Returns True if both work and returns False.

    Inputs:
    :param dict params: dict of user variables which will govern how the programs runs
    :param str vina_exe: path to vina executable
    :param str qvina2_exe: path to quick vina 2 executable
    Returns:
    :returns: bool bool: returns True if both docking executables work; False if either fails
    """
    test_vina_outfile = (
        params["root_output_folder"] + os.sep + "docking_exe_MACOS_test.txt"
    )

    try:
        command = f"{vina_exe} --version > {test_vina_outfile} 2>> {test_vina_outfile}"
        os.system(command)
        command = (
            f"{qvina2_exe} --version >> {test_vina_outfile} 2>> {test_vina_outfile}"
        )
        os.system(command)
    except Exception:
        printout = "Docking executables could not be found."
        # is not compatible on this OS. \nPlease use docker "
        return False

    with open(test_vina_outfile, "r") as test_file:
        lines = test_file.readlines()
    if "AutoDock Vina 1.1.2" not in lines[0]:
        printout = (
            "Vina docking is not compatible on this OS. \nPlease use docker or "
            + "try provide a Vina executable compatible with the OS.\n"
        )
        print(printout)
        if params["dock_choice"] == "VinaDocking":
            return False

    if "QuickVina 2.1" not in lines[1]:
        printout = (
            "QuickVina 2.1 docking is not compatible on this OS. \nPlease use docker"
            + " or try provide a Vina executable compatible with the OS.\n"
        )
        print(printout)
        if params["dock_choice"] == "QuickVina2Docking":
            return False
    return True
