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
          -  ./sample_autogrow_docker_json.json

    Required variables within the JSON file:
    - `-root_output_folder`: folder path on host system that results will be copied to.
    - `-source_compound_file`: Path on host system to the tab-delineate .smi
        file that will seed generation 1.
    - `-filename_of_receptor`: Path on host system of the receptor to be tested.
    - `-center_x`, `-center_y`, `-center_z`: x,y,z coordinates of center of pocket to be tested.
    - `-size_x`, `-size_y`, `-size_z`: dimensions of the pocket in x,y,z coordinates.
    Variable that will be ignored:
    - `-openbabel_bin_directory` should not be specified.
    - `-mgltools_directory` should not be specified.

The resulting AutoGrow4 output will be zipped and transfered to the desired
root_output_folder.

An example JSON is provided in: ./sample_autogrow_docker_json.json


To run AutoGrow4 in a docker, please run the `autogrow_in_docker.py` script:
    Example on Linux/MacOS:
        #  cd to this directory in a bash terminal
        1) cd autogrow4/Docker/
        # Run autogrow_in_docker.py with sudo and supply a json file using the
        # normal pathing of your system.
        # Please note that the docker downloads its own copy of obabel and MGLTools
        # so you do not need to provide those paths.
        2) `sudo python autogrow_in_docker.py -j ./sample_autogrow_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

    Example on Windows OS:
        1) open a docker enabled and bash enabled terminal with administrative priveledges
        #  cd to this directory in a bash terminal
        3) cd autogrow4/Docker/
        4)  `python autogrow_in_docker.py -j ./sample_autogrow_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

"""
import os
import shutil
import json
import argparse

def make_docker():
    """
    This will create the docker to run AutoGrow4.
    This is also where all of the files are copied into the image.

    If docker image can not be created it will raise an exception.
    """
    print("Creating new docker image for AutoGrow4")
    script_dir = str(os.path.dirname(os.path.realpath(__file__))) + os.sep
    output_and_log_dir = os.path.abspath("output_and_log_dir") + os.sep
    log_file = "{}log.txt".format(output_and_log_dir)
    printout = "\nAttempting to create the docker container. If 1st time running " + \
        "this script it may take a few minutes. Output details are piped to: " + \
        "{}\n".format(log_file)

    build_bash_script = script_dir + "build.bash"
    try:
        os.system("cd {}".format(script_dir))
        os.system("bash {} {} > {}".format(build_bash_script, script_dir, log_file))
    except:
        printout = "\nCan not create a docker file. Please make sure to run the " + \
            "script with sudo/administrative priveledges.\nIf Linux/MacOS:" + \
            "\t'sudo bash build.bash'\nnIf Windows:\tRun from bash and " + \
            "docker enabled terminal with administrative priveledges.\n" + \
            "Please also make sure docker is installed on the system."
        print(printout)
        raise Exception(printout)
#
def check_for_required_inputs(json_vars):
    """
    Confirm all the required inputs were provided.

    Required Variables go here.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.
    """
    keys_from_input = list(json_vars.keys())

    list_of_required_inputs = [
        "filename_of_receptor",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "root_output_folder",
        "source_compound_file",
    ]

    missing_variables = []
    for variable in list_of_required_inputs:
        if variable in keys_from_input:
            continue
        missing_variables.append(variable)

    if len(missing_variables) != 0:
        printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunAutogrow -h"
        printout = printout + "\nThe following required variables are missing: "
        for variable in missing_variables:
            printout = printout + "\n\t" + variable
        print("")
        print(printout)
        print("")
        raise NotImplementedError("\n" + printout + "\n")

    # Make sure the dimmensions are in floats. If in int convert to float.
    for x in ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]:
        if type(json_vars[x]) in [float, int]:
            continue
        printout = "\n{} must be a float value.\n".format(x)
        print(printout)
        raise Exception(printout)


    #######################################
    # Check that all required files exist #
    #######################################

    # convert paths to abspath, in case necessary
    json_vars["filename_of_receptor"] = os.path.abspath(
        json_vars["filename_of_receptor"]
    )
    json_vars["root_output_folder"] = os.path.abspath(
        json_vars["root_output_folder"]
    )
    json_vars["source_compound_file"] = os.path.abspath(
        json_vars["source_compound_file"]
    )
    # Check filename_of_receptor exists
    if os.path.isfile(json_vars["filename_of_receptor"]) is False:
        raise NotImplementedError(
            "Receptor file can not be found. File must be a .PDB file."
        )
    if ".pdb" not in json_vars["filename_of_receptor"]:
        raise NotImplementedError("filename_of_receptor must be a .PDB file.")

    # Check root_output_folder exists
    if os.path.exists(json_vars["root_output_folder"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(json_vars["root_output_folder"])
            os.system("chmod -R a+rwx {}".format(json_vars["root_output_folder"]))
        except:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

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
            {}.".format(json_vars["source_compound_file"])
        )
    if ".smi" not in json_vars["source_compound_file"]:
        raise NotImplementedError(
            "source_compound_file must be a \
            tab delineated .smi file."
        )
#
def find_previous_runs(folder_name_path):
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
    while path_exists is True:
        folder_path = "{}{}{}".format(folder_name_path, i, os.sep)
        if os.path.exists(folder_path):
            i = i + 1
        else:
            path_exists = False

    if i == 0:
        # There are no previous runs in this directory
        last_run_number = None
        return None

    # A previous run exists. The number of the last run.
    last_run_number = i - 1
    return last_run_number
#
def set_run_directory(root_folder_path, start_a_new_run):
    """
    Determine and make the folder for the run directory.
        If start_a_new_run is True    Start a frest new run.
            -If no previous runs exist in the root_folder_path then make a new
                folder named root_folder_path + "Run_0"
            -If there are previous runs in the root_folder_path then make a
                new folder incremental increasing the name by 1 from the last
                run in the same output directory.
        If start_a_new_run is False    Find the last run folder and return that path
            -If no previous runs exist in the root_folder_path then make a new
            folder named root_folder_path + "Run_0"

    Inputs:
    :param str root_folder_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files
    :param bol start_a_new_run: True or False to determine if we continue from
        the last run or start a new run
        - This is set as a vars["start_a_new_run"]
        - The default is vars["start_a_new_run"] = True
    Returns:
    :returns: str folder_path: the string of the newly created directory for
        puting output folders
    """

    folder_name_path = root_folder_path + "Run_"
    print(folder_name_path)

    last_run_number = find_previous_runs(folder_name_path)

    if last_run_number is None:
        # There are no previous simulation runs in this directory
        print("There are no previous runs in this directory.")
        print("Starting a new run named Run_0.")

        # make a folder for the new generation
        run_number = 0
        folder_path = "{}{}{}".format(folder_name_path, run_number, os.sep)
        os.makedirs(folder_path)

    else:
        if start_a_new_run is False:
            # Continue from the last simulation run
            run_number = last_run_number
            folder_path = "{}{}{}".format(folder_name_path, last_run_number, os.sep)
        else:  # start_a_new_run is True
            # Start a new fresh simulation
            # Make a directory for the new run by increasing run number by +1
            # from last_run_number
            run_number = last_run_number + 1
            folder_path = "{}{}{}".format(folder_name_path, run_number, os.sep)
            os.mkdir(folder_path)
            os.system("chmod -R a+rwx {}".format(folder_path))

    if run_number == last_run_number:
        temp_dir_path = os.path.abspath("temp_user_files") + os.sep
        temp_path = temp_dir_path + "old_runs" + os.sep + "Run_0"
        if os.path.exists(temp_path) is True:
            shutil.rmtree(temp_path)

        shutil.copytree(folder_path, temp_path)

    print("The Run number is: ", run_number)
    print("The Run folder path is: ", folder_path)
    print("")
    return folder_path
#
def get_output_folder(json_vars):
    """
    Find the folder for where to place output runs on host system.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.
    Returns:
    :returns: str folder_path: the string of the newly created directory for
        puting output folders
    """
    if "start_a_new_run" in json_vars.keys():
        start_a_new_run = json_vars["start_a_new_run"]
    else:
        start_a_new_run = False

    root_output_folder = os.path.abspath(json_vars["root_output_folder"]) + os.sep
    folder_path = set_run_directory(root_output_folder, start_a_new_run)

    return folder_path
#
def move_files_to_temp_dir(json_vars):
    """
    This will move all files needed to a temp_user_files directory and will created a modified
    json_vars dict called docker_json_vars which will be used for pathing within
    the docker.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.
    """
    docker_json_vars = {}
    # make or remove and make the temp_user_files dir
    temp_dir_path = os.path.abspath("temp_user_files") + os.sep
    if os.path.exists(temp_dir_path):
        shutil.rmtree(temp_dir_path)
    os.mkdir(temp_dir_path)
    os.system("chmod -R a+rwx {}".format(temp_dir_path))

    # make or remove and make an output_and_log_dir
    output_and_log_dir = os.path.abspath("output_and_log_dir") + os.sep
    if os.path.exists(output_and_log_dir):
        shutil.rmtree(output_and_log_dir)
    os.mkdir(output_and_log_dir)
    os.system("chmod -R a+rwx {}".format(output_and_log_dir))

    print("copying files into temp directory: temp_user_files")
    # get files from json_vars
    for var_name in json_vars.keys():
        var_item = json_vars[var_name]
        if str(type(var_item)) not in ["<type 'unicode'>", "<type 'unicode'>"]:
            continue
        var_item = str(var_item)
        # This could be a different variable that is not a path
        # ie) dock_choice: QuickVina2 would be a string that is not a path
        if os.path.exists(var_item) is False:
            continue
        if "mgl" in var_name.lower():
            print("MGLTools from within the docker will be used")
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
            docker_json_vars[var_name] = "/UserFiles/" + basename + "/"
            continue

        if os.path.isfile(var_item):
            shutil.copyfile(var_item, temp_path)
            docker_json_vars[var_name] = "/UserFiles/" + basename

    for var_name in json_vars.keys():
        if var_name not in docker_json_vars.keys():
            docker_json_vars[var_name] = json_vars[var_name]

    # Add docker babel and MGL paths
    docker_json_vars["mgltools_directory"] = "/mgltools_x86_64Linux2_1.5.6"
    docker_json_vars["obabel_path"] = "/usr/bin/obabel"

    # Set output folder
    docker_json_vars["root_output_folder"] = "/Outputfolder/"

    with open(temp_dir_path + "docker_json_vars.json", "w") as file_item:
        json.dump(docker_json_vars, file_item, indent=4)

    # update permissions so files can be manipulated without sudo/admin
    os.system("chmod -R a+rwx {}".format(temp_dir_path))
    os.system("chmod -R a+rwx {}".format(output_and_log_dir))
#
def handle_json_info(vars):
    """
    This will open the json file.
        1) check that JSON file has basic info
            -receptor, size/center...
        2) copy files to a temp directory
            -receptor, .smi files ...
        3) make a JSON file with modified information for within docker
    Inputs:
    :param dict vars: Dictionary of User specified variables

    Returns:
    :param dict json_vars: Dictionary of User specified variables
    :returns: str folder_path: the string of the newly created directory for
        puting output folders
    """
    print("Handling files")

    json_file = vars["json_file"]
    if os.path.exists(json_file) is False:
        printout = "\njson_file is required. Can not find json_file: {}.\n".format(json_file)
        print(printout)
        raise Exception(printout)
    json_vars = json.load(open(json_file))
    check_for_required_inputs(json_vars)
    move_files_to_temp_dir(json_vars)

    # get output folder
    outfolder_path = get_output_folder(json_vars)

    return json_vars, outfolder_path
#
def run_autogrow_docker_main(vars):
    """
    This function runs the processing to:
        1) check that JSON file has basic info
            -receptor, size/center...
        2) copy files to a temp directory
            -receptor, .smi files ...
        3) make a JSON file with modified information for within docker
        4) run build.bash
            which transfers the necessary files to docker container
        5) execute RunAutogrow.py from within the docker containiner
        6) export the files back to the final end dir

    Inputs:
    :param dict vars: Dictionary of User specified variables
    """
    printout = "\n\nThis script builds a docker for AutoGrow4 and runs AutoGrow4 " + \
        "within the docker. The setup may take a few minutes the first time being run " + \
        "and AutoGrow may take a long time depending on the settings.\n\n"
    print(printout)

    # Check that we are in the correct directory if not raise exception
    script_dir = str(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if os.path.abspath(os.getcwd()) != os.path.abspath(script_dir):
        printout = "\nMust execute this script from this directory: {}\n".format(script_dir)
        printout = printout + "Before running please 'cd {}'\n".format(script_dir)
        print(printout)
        raise Exception(printout)

    # Run parts 1-3
    # 1) check that JSON file has basic info
    #     -receptor, size/center...
    # 2) copy files to a temp directory
    #     -receptor, .smi files ...
    # 3) make a JSON file with modified information for within docker
    json_vars, outfolder_path = handle_json_info(vars)

    # Run part 4) run build.bash
    make_docker()

    # Run part 5) run AutoGrow in the container
    # docker cp foo.txt mycontainer:/foo.txt
    print("\nRunning AutoGrow4 in Docker")
    tmp_path = os.path.abspath("temp_user_files")
    script_dir = str(os.path.dirname(os.path.realpath(__file__))) + os.sep
    execute_outside_docker = script_dir + "execute_autogrow_from_outside_docker.sh"
    command = "bash {} {} {}".format(execute_outside_docker, tmp_path, outfolder_path)
    os.system(command)
    print("AutoGrow Results placed in: {}".format(outfolder_path))
#

PARSER = argparse.ArgumentParser()

# Allows the run commands to be submitted via a .json file.
PARSER.add_argument(
    "--json_file",
    "-j",
    metavar="param.json_file",
    required=True,
    help="Name of a json file containing all parameters. \
    Overrides other arguments. This takes all the parameters described in \
    RunAutogrow.py. MGLTools and openbabel paths can be ignored as they are \
    already installed in the docker image.",
)

ARGS_DICT = vars(PARSER.parse_args())
run_autogrow_docker_main(ARGS_DICT)
