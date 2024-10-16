#!/usr/bin/env python
"""
This script will handle running AutoGrow4 in a docker container.
It handles generating the docker image and container, handling the user variables,
executing AutoGrow4, and copying the files from the container to the desired directory.

This script requires a JSON file that contains all the parameters that would be
required to run AutoGrow4 on a host system (ie paths on your computer).
Necessary files, such as the receptor pdb, will be copied into the docker.

To run AutoGrow from within docker. Launches docker
  image. Accepts the exact same parameters as AutoGrow4, with the following
  exceptions:
    1) User variables must be supplied in JSON format.
        - Please see documentation within the tutorial manual and an example can be found:
          -  ./examples/sample_autogrow_docker_json.json

    Required variables within the JSON file:
    - `-root_output_folder`: folder path on host system that results will be copied to.
    - `-source_compound_file`: Path on host system to the tab-delineate .smi
        file that will seed generation 1.
    - `-receptor_path`: Path on host system of the receptor to be tested.
    - `-center_x`, `-center_y`, `-center_z`: x,y,z coordinates of center of pocket to be tested.
    - `-size_x`, `-size_y`, `-size_z`: dimensions of the pocket in x,y,z coordinates.
    Variable that will be ignored:
    - `-openbabel_bin_directory` should not be specified.
    - `-mgltools_directory` should not be specified.

The resulting AutoGrow4 output to the desired root_output_folder.

An example JSON is provided in: ./sample_autogrow_docker_json.json


To run AutoGrow4 in a docker, please run the `autogrow_in_docker.py` script:
    Example on Linux/MacOS:
        #  cd to this directory in a bash terminal
        1) cd autogrow4/docker/
        # Run autogrow_in_docker.py with sudo and supply a json file using the
        # normal pathing of your system.
        # Please note that the docker downloads its own copy of obabel and MGLTools
        # so you do not need to provide those paths.
        2) `sudo python autogrow_in_docker.py -j ./examples/sample_autogrow_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

    Example on Windows OS:
        1) open a docker enabled and bash enabled terminal with administrative privileges
        #  cd to this directory in a bash terminal
        3) cd autogrow4/docker/
        4)  `python autogrow_in_docker.py -j ./examples/sample_autogrow_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

"""
import __future__

import os
import shutil
import json
import argparse
import sys
from typing import Any, Dict, List, Optional


def change_permissions(file_path: str) -> None:
    """
    This will open the permissions for a given file.

    Inputs:
    :param str file_path: Path to a file to open permissions to.
    """
    os.chmod(file_path, 0o777)


def change_permissions_recursively(file_or_folder_path: str) -> None:
    """
    This will open the permissions for a given file/folder.

    Skip permissions change if Windows.

    Inputs:
    :param str file_or_folder_path: Path to a file/folder to open permissions to.
    """
    if os.name in {"nt", "ce"}:
        return
    if sys.platform.lower() in ["linux", "linux2"]:
        # chmod -R recursively open the permissions
        os.system(f"chmod -R a+rwx {file_or_folder_path}")
    elif sys.platform.lower() == "darwin":
        # chmod -R recursively open the permissions
        os.system(f"chmod -R a+rwx {file_or_folder_path}")
    elif os.path.isdir(file_or_folder_path):
        # chmod may not be a valid command on other OS systems. So let's do this
        #  the manual way.
        _apply_permissions_recursively(file_or_folder_path)
    else:
        change_permissions(file_or_folder_path)


def _apply_permissions_recursively(file_or_folder_path: str) -> None:
    directory_path_list = []
    file_list = []
    for top_dir, dir_list, list_of_files in os.walk(file_or_folder_path, topdown=False):
        directory_path_list.extend(iter([os.path.join(top_dir, d) for d in dir_list]))
        file_list.extend(iter([os.path.join(top_dir, fil) for fil in list_of_files]))
    # Convert mods on all files within a directory
    file_list = list(set(file_list))
    directory_path_list = list(set(directory_path_list))
    for file_path in file_list:
        change_permissions(file_path)
    for dir_path in directory_path_list:
        change_permissions(dir_path)


def adjust_dockerfile() -> None:
    """
    This will open Dockerfile and check if the entrypoint has been switched
    to the windows version (run_autogrow_in_container.bash) if not it will
    modify the Dockerfile to use the windows version of the script.

    This only should run on Windows OS.

    Change:
        ENTRYPOINT ["bash", "/autogrow/run_autogrow_in_container.bash"]
    To:
        # ENTRYPOINT ["bash", "/autogrow/run_autogrow_in_container.bash"]
        ENTRYPOINT ["bash", "/autogrow/run_autogrow_in_container_windows.bash"]
    """
    printout = ""
    normal_entry = "/autogrow/run_autogrow_in_container.bash"
    windows_entry = "/autogrow/run_autogrow_in_container_windows.bash"
    replacement_line = (
        'ENTRYPOINT ["bash", "/autogrow/run_autogrow_in_container_windows.bash"]\n'
    )
    print("Modifying the Dockerfile to run for Windows. Changing Entrypoint.")
    with open(os.path.abspath("Dockerfile"), "r") as f:
        for line in f:
            if "ENTRYPOINT" in line and normal_entry in line:
                if "#" not in line:
                    line = f"# {line}" + "\n"
                printout = printout + line
                printout = printout + replacement_line
                continue
            if "ENTRYPOINT" in line and windows_entry in line:
                continue
            printout = printout + line
    with open(os.path.abspath("Dockerfile"), "w") as f:
        f.write(printout)


def make_docker() -> None:
    """
    This will create the docker to run AutoGrow4.
    This is also where all of the files are copied into the image.

    If docker image can not be created it will raise an exception.
    """
    if os.name in {"nt", "ce"}:
        # so it's running under windows. multiprocessing disabled
        adjust_dockerfile()

    print("Creating new docker image for AutoGrow4")
    output_and_log_dir = os.path.abspath("output_and_log_dir") + os.sep
    log_file = f"{output_and_log_dir}log.txt"
    printout = (
        "\nAttempting to create the docker container. If 1st time running "
        + "this script it may take a few minutes. Output details are piped to: "
        + f"{log_file}\n"
    )

    print(printout)
    try:
        os.system(f"docker build -t autogrow4 . > {log_file}")
    except Exception as e:
        printout = (
            "\nCan not create a docker file. Please make sure to run the "
            + "script with sudo/administrative privileges.\nIf Linux/MacOS:\n"
            + "\t'sudo python autogrow_in_docker.py -j PATH/JSON.json'\n"
            + "If Windows:\n\tRun from bash and "
            + "docker enabled terminal with administrative privileges.\n"
            + "Please also make sure docker is installed on the system."
        )
        print(printout)
        raise Exception(printout) from e

    # Remove the temporary autogrow4 directory
    shutil.rmtree("autogrow4")


def check_for_required_params(json_vars: Dict[str, Any]) -> Dict[str, Any]:
    """
    Confirm all the required inputs were provided.

    Required Variables go here.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name:
                           value}.

    Returns:
    :returns: dict json_vars: The updated json_vars with input-file paths
                              changed to the output dir/inputs/ subdirectory.
    """

    keys_from_input = list(json_vars.keys())

    list_of_required_inputs = [
        "receptor_path",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "root_output_folder",
        "source_compound_file",
    ]

    missing_variables = [
        variable
        for variable in list_of_required_inputs
        if variable not in keys_from_input
    ]
    if missing_variables:
        _handle_missing_required_params(missing_variables)
    # Make sure the dimmensions are in floats. If in int convert to float.
    for x in ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]:
        if type(json_vars[x]) in [float, int]:
            continue
        printout = f"\n{x} must be a float value.\n"
        print(printout)
        raise Exception(printout)

    #######################################
    # Check that all required files exist #
    #######################################

    # convert paths to abspath, in case necessary
    json_vars["receptor_path"] = os.path.abspath(json_vars["receptor_path"])
    json_vars["root_output_folder"] = os.path.abspath(json_vars["root_output_folder"])
    json_vars["source_compound_file"] = os.path.abspath(
        json_vars["source_compound_file"]
    )
    # Check receptor_path exists
    if os.path.isfile(json_vars["receptor_path"]) is False:
        raise NotImplementedError(
            "Receptor file can not be found. File must be a .PDB file."
        )
    if ".pdb" not in json_vars["receptor_path"]:
        raise NotImplementedError("receptor_path must be a .PDB file.")

    # Check root_output_folder exists
    if os.path.exists(json_vars["root_output_folder"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(json_vars["root_output_folder"])
            os.makedirs(json_vars["root_output_folder"] + os.sep + "inputs")
            change_permissions_recursively(json_vars["root_output_folder"])
        except Exception as e:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            ) from e

        if os.path.exists(json_vars["root_output_folder"]) is False:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

    if os.path.isdir(json_vars["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )

    # Check source_compound_file exists
    if os.path.isfile(json_vars["source_compound_file"]) is False:
        raise NotImplementedError(
            "source_compound_file must be a tab delineated .smi file. \
            source_compound_file can not be found: \
            {}.".format(
                json_vars["source_compound_file"]
            )
        )
    if ".smi" not in json_vars["source_compound_file"]:
        raise NotImplementedError(
            "source_compound_file must be a \
            tab delineated .smi file."
        )

    # You need to copy the input files to the output directory, so it can be
    # easily access from the docker container. (JDD addition.)
    shutil.copy2(
        json_vars["source_compound_file"],
        json_vars["root_output_folder"]
        + os.sep
        + "inputs"
        + os.sep
        + os.path.basename(json_vars["source_compound_file"]),
    )
    json_vars["source_compound_file"] = "/Outputfolder/inputs/" + os.path.basename(
        json_vars["source_compound_file"]
    )
    shutil.copy2(
        json_vars["receptor_path"],
        json_vars["root_output_folder"]
        + os.sep
        + "inputs"
        + os.sep
        + os.path.basename(json_vars["receptor_path"]),
    )
    json_vars["receptor_path"] = "/Outputfolder/inputs/" + os.path.basename(
        json_vars["receptor_path"]
    )

    return json_vars


def _handle_missing_required_params(missing_variables: List[str]) -> None:
    printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunAutogrow -h"
    printout += "\nThe following required variables are missing: "
    for variable in missing_variables:
        printout = printout + "\n\t" + variable
    print("")
    print(printout)
    print("")
    raise NotImplementedError("\n" + printout + "\n")


def find_previous_runs(folder_name_path: str) -> Optional[int]:
    """
    This will check if there are any previous runs in the output directory.
        - If there are it will return the interger of the number label of the last Run folder path.
            - ie if there are folders Run_0, Run_1, Run_2 the function will return int(2)
        - If there are no previous Run folders it returns None.

    Inputs:
    :param str folder_name_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files

    Returns:
    :returns: int last_run_number: the int of the last run number or None if no previous runs.
    """

    path_exists = True
    i = 0
    while path_exists:
        folder_path = f"{folder_name_path}{i}{os.sep}"
        if os.path.exists(folder_path):
            i = i + 1
        else:
            path_exists = False

    if i == 0:
        # There are no previous runs in this directory
        last_run_number = None
        return None

    # A previous run exists. The number of the last run.
    return i - 1


def get_run_number(root_folder_path: str, start_a_new_run: bool) -> str:
    """
    Determine run number for the new directory.
        If start_a_new_run is True    Start a fresh new run.
            -If no previous runs exist in the root_folder_path
            -If there are previous runs in the root_folder_path
                incremental increasing the name by 1 from the last
                run in the same output directory.
        If start_a_new_run is False    Find the last run folder and return that path
            -If no previous runs exist in the root_folder_path then  "Run_0"

    Inputs:
    :param str root_folder_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files
    :param bol start_a_new_run: True or False to determine if we continue from
        the last run or start a new run
        - This is set as a params["start_a_new_run"]
        - The default is params["start_a_new_run"] = True
    Returns:
    :returns: str run_num: the string of the run number "Run_*"
    """

    folder_name_path = f"{root_folder_path}Run_"
    print(folder_name_path)

    last_run_number = find_previous_runs(folder_name_path)

    if last_run_number is None:
        # There are no previous simulation runs in this directory
        print("There are no previous runs in this directory.")
        print("Starting a new run named Run_0.")

        # make a folder for the new generation
        run_number = 0

    elif not start_a_new_run:
        # Continue from the last simulation run
        run_number = last_run_number
    else:  # start_a_new_run is True
        # Start a new fresh simulation
        # Make a directory for the new run by increasing run number by +1
        # from last_run_number
        run_number = last_run_number + 1

    folder_name_path = f"{root_folder_path}Run_{str(run_number)}{os.sep}"
    print("The Run number is: ", run_number)
    print("The Run folder path is: ", folder_name_path)
    print("")
    return f"Run_{run_number}"


def get_output_folder(json_vars: Dict[str, Any]) -> tuple:
    """
    Find the folder for where to place output runs on host system.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.
    Returns:
    :returns: str root_output_folder: the string of the directory for
        puting output folders
    :returns: str run_num: the string of the run number "Run_*"
    """
    start_a_new_run = json_vars.get("start_a_new_run", False)
    root_output_folder = os.path.abspath(json_vars["root_output_folder"]) + os.sep
    change_permissions_recursively(root_output_folder)
    run_num = get_run_number(root_output_folder, start_a_new_run)

    return root_output_folder, run_num


def move_files_to_temp_dir(json_vars: Dict[str, Any]) -> None:
    """
    This will move all files needed to a temp_user_files directory and will created a modified
    json_vars dict called docker_json_vars which will be used for pathing within
    the docker.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.
    """

    docker_json_vars = {}
    temp_dir_path = _setup_temp_directory("temp_user_files")
    output_and_log_dir = _setup_temp_directory("output_and_log_dir")
    print("copying files into temp directory: temp_user_files")
    # get files from json_vars
    for var_name, var_item in json_vars.items():
        if str(type(var_item)) not in ["<type 'unicode'>", "<type 'unicode'>"]:
            continue
        var_item = str(var_item)
        # This could be a different variable that is not a path
        # ie) dock_choice: QuickVina2 would be a string that is not a path
        if os.path.exists(var_item) is False:
            continue
        if "babel" in var_name.lower():
            print("obabel from within the docker will be used")
            continue
        if var_name == "root_output_folder":
            continue
        basename = os.path.basename(var_item)
        temp_path = temp_dir_path + basename
        if os.path.isdir(var_item):
            shutil.copytree(var_item, temp_path)
            docker_json_vars[var_name] = f"/UserFiles/{basename}/"
            continue

        if os.path.isfile(var_item):
            shutil.copyfile(var_item, temp_path)
            docker_json_vars[var_name] = f"/UserFiles/{basename}"

    for var_name in json_vars:
        if var_name not in docker_json_vars.keys():
            docker_json_vars[var_name] = json_vars[var_name]

    # Add docker babel and MGL paths
    docker_json_vars["obabel_path"] = "/usr/bin/obabel"

    # Set output folder
    docker_json_vars["root_output_folder"] = "/Outputfolder/"

    with open(f"{temp_dir_path}docker_json_vars.json", "w") as file_item:
        json.dump(docker_json_vars, file_item, indent=4)

    # update permissions so files can be manipulated without sudo/admin
    change_permissions_recursively(temp_dir_path)
    change_permissions_recursively(output_and_log_dir)

    # Copy over AutoGrow4 files into a temp directory
    temp_autogrow4_path = os.path.abspath("autogrow4") + os.sep

    script_dir = str(os.path.dirname(os.path.realpath(__file__)))
    autogrow4_top_dir = str(os.path.dirname(script_dir))
    if os.path.exists(temp_autogrow4_path):
        shutil.rmtree(temp_autogrow4_path)
    os.mkdir(temp_autogrow4_path)
    autogrow4_top_dir += os.sep
    change_permissions_recursively(temp_autogrow4_path)

    # Copy all files in autogrow4 directory into a temp except the Docker folder
    for fol_to_copy in [
        "autogrow",
        "source_compounds",
        "accessory_scripts",
        "tutorial",
    ]:
        shutil.copytree(
            autogrow4_top_dir + fol_to_copy, temp_autogrow4_path + fol_to_copy
        )
    shutil.copyfile(
        f"{autogrow4_top_dir}run_autogrow.py", f"{temp_autogrow4_path}run_autogrow.py",
    )
    # Open permissions
    change_permissions_recursively(temp_autogrow4_path)


def _setup_temp_directory(arg0: str) -> str:
    # make or remove and make the temp_user_files dir
    result = os.path.abspath(arg0) + os.sep
    if os.path.exists(result):
        shutil.rmtree(result)
    os.mkdir(result)
    change_permissions_recursively(result)

    return result


def handle_json_info(params: Dict[str, Any]) -> tuple:
    """
    This will open the json file.
        1) check that JSON file has basic info
            -receptor, size/center...
        2) copy files to a temp directory
            -receptor, .smi files ...
        3) make a JSON file with modified information for within docker
    Inputs:
    :param dict params: Dictionary of User specified variables

    Returns:
    :param dict json_vars: Dictionary of User specified variables
    :returns: str root_output_folder: the string of the directory for
        puting output folders
    :returns: str run_num: the string of the run number "Run_*"
    """
    print("Handling files")

    json_file = params["json_file"]
    if os.path.exists(json_file) is False:
        printout = f"\njson_file is required. Can not find json_file: {json_file}.\n"
        print(printout)
        raise Exception(printout)
    json_vars = json.load(open(json_file))
    json_vars = check_for_required_params(json_vars)
    move_files_to_temp_dir(json_vars)

    # get output folder
    outfolder_path, run_num = get_output_folder(json_vars)

    return json_vars, outfolder_path, run_num


def run_autogrow_docker_main(params: Dict[str, Any]) -> None:
    """
    This function runs the processing to:
        1) check that JSON file has basic info
            -receptor, size/center...
        2) copy files to a temp directory
            -receptor, .smi files ...
        3) make a JSON file with modified information for within docker
        4) Build docker image and link files to output folder
            -This includes an adjustment to the Dockerfile if
            running it on a Windows OS
        5) execute run_autogrow.py from within the docker container
        6) export the files back to the final end dir

    Inputs:
    :param dict params: Dictionary of User specified variables
    """
    printout = (
        "\n\nThis script builds a docker for AutoGrow4 and runs AutoGrow4 "
        + "within the docker. The setup may take a few minutes the first time being run "
        + "and AutoGrow may take a long time depending on the settings.\n\n"
    )
    print(printout)

    # Check that we are in the correct directory if not raise exception
    script_dir = str(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if os.path.abspath(os.getcwd()) != os.path.abspath(script_dir):
        printout = (
            f"\nMust execute this script from this directory: {script_dir}\n"
            + f"Before running please 'cd {script_dir}'\n"
        )
        print(printout)
        raise Exception(printout)

    # Run parts 1-3
    # 1) check that JSON file has basic info
    #     -receptor, size/center...
    # 2) copy files to a temp directory
    #     -receptor, .smi files ...
    # 3) make a JSON file with modified information for within docker
    json_vars, outfolder_path, run_num = handle_json_info(params)

    # Run build docker image
    make_docker()

    # Run part 5) run AutoGrow in the container
    print("\nRunning AutoGrow4 in Docker")

    command = (
        f"docker run --rm -it -v {outfolder_path}:/Outputfolder/"
        + f" autogrow4  --name autogrow4  --{run_num}"
    )
    # Execute AutoGrow4
    print(command)
    os.system(command)
    change_permissions_recursively(outfolder_path)
    print(f"AutoGrow Results placed in: {outfolder_path}")


PARSER = argparse.ArgumentParser()

# Allows the run commands to be submitted via a .json file.
PARSER.add_argument(
    "--json_file",
    "-j",
    metavar="param.json_file",
    required=True,
    help="Name of a json file containing all parameters. \
    Overrides other arguments. This takes all the parameters described in \
    run_autogrow.py. MGLTools and openbabel paths can be ignored as they are \
    already installed in the docker image.",
)
PARSER.add_argument(
    "--override_sudo_admin_privileges",
    metavar="param.override_sudo_admin_privileges",
    default=False,
    help="Docker normally requires `sudo` (linux/macos) or `Administrator` \
    privileges (windows/cygwin). If an system does not have such privileges, \
    or does not require such privileges, setting this to True, will skip the \
    check for privileges. This variable is provided via commandline, and \
    IS NOT RECOMMENDED for most OS. ",
)
ARGS_DICT = vars(PARSER.parse_args())

print("")
print("BE SURE TO RUN THIS SCRIPT WITH SUDO (LINUX/MACOS) OR ADMINISTRATOR")
print("(WINDOWS) PRIVILEGES!")
print("")

# Check that this is running with appropriate privileges.
# i.e., sudo (linux/macos) or Administrator privileges (Windows/cygwin)
if ARGS_DICT["override_sudo_admin_privileges"] == False:
    if sys.platform.lower() in ["darwin", "linux", "linux2"]:
        if os.getuid() != 0:
            printout = (
                "\n\nMust run this script with `sudo` privileges.\n\t"
                + "Please retry running with `sudo` privileges.\n\n"
            )
            print(printout)
            raise Exception(printout)
    elif sys.platform.lower() in ["win32", "cygwin"]:
        import ctypes

        if ctypes.windll.shell32.IsUserAnAdmin() != 1:  # type: ignore
            printout = (
                "\n\nMust run this script from a terminal with `Administrator` privileges.\n\t"
                + "Please retry running with `Administrator` privileges.\n\n"
            )
            print(printout)
            raise Exception(printout)
    else:
        print("")
        print("BE SURE TO RUN THIS SCRIPT WITH SUDO (LINUX/MACOS) OR ADMINISTRATOR")
        print("(WINDOWS) PRIVILEGES!")
        print("")
else:

    print("\n##############################################################")
    print("WARNING: Skipping check for privileges.")
    print("\tBE SURE TO RUN THIS SCRIPT WITH APPROPRIATE PRIVILEGES:")
    print("\tSUDO (LINUX/MACOS) OR ADMINISTRATOR (WINDOWS) PRIVILEGES!")
    print("\tFailure to do so may result in Docker failures.")
    print("##############################################################\n")


run_autogrow_docker_main(ARGS_DICT)
