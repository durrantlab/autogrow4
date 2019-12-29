#!/usr/bin/env python

import os
import sys
import shutil
import json
import argparse
import glob

def test_docker_exists():
    """
    This will test whether a docker image name autogrow exists
    and it will name the container autogrow. making sure the container
    is named autogrow will ensure we are able to copy the files into
    the container later on.

    If docker image has not been created it will raise an exception.
    """
    try:
        os.system("docker run autogrow --name autogrow -d -test")
        os.system("docker commit ba0d61f03294 autogrow")
    except:
        printout = "\nCan not find docker file. Please run 'sudo bash build.bash'\n"
        print(printout)
        raise Exception(printout)
#
############################################
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
        else:
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
def move_files_to_temp_dir(json_vars):
    """
    This will move all files needed to a tmp directory and will created a modified
    json_vars dict called docker_json_vars which will be used for pathing within
    the docker.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.

    Returns:
    :returns: dict docker_json_vars: A modified version of the json dictionary
        that is to be used within the docker container.
    """
    docker_json_vars = {}
    # make or remove and make the tmp dir
    temp_dir_path = os.path.abspath(".tmp") + os.sep
    if os.path.exists(temp_dir_path):
        shutil.rmtree(temp_dir_path)
    os.mkdir(temp_dir_path)

    # get files from json_vars
    for var_name in json_vars.keys():
        var_item = json_vars[var_name]
        if type(var_item) != str:
            continue
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
            docker_json_vars[var_name] = "/UserFiles/" + basename + os.sep
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

    return docker_json_vars
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
    :param dict argv: Dictionary of User specified variables

    Returns:
    :param dict argv: Dictionary of User specified variables
    :param dict json_vars: Dictionary of User specified variables
    :returns: dict docker_json_vars: A modified version of the json dictionary
        that is to be used within the docker container.
    """
    json_file = vars["json_file"]
    if os.path.exists(json_file) is False:
        printout = "\njson_file is required. Can not find json_file: {}.\n".format(json_file)
        print(printout)
        raise Exception(printout)

    json_vars = json.load(open(json_file))
    check_for_required_inputs(json_vars)
    docker_json_vars = move_files_to_temp_dir(json_vars)

    # mover files from tmp dir to docker
    # docker cp foo.txt mycontainer:/foo.txt
    temp_dir = os.path.abspath(".tmp") + os.sep + "*"
    for file_stuff in glob.glob(temp_dir):
        #@@@@@JAKE NEED TO FIND A WAY TO REPLACE THE CONTAINER ID
        os.system("sudo docker cp  {} autogrow:/UserFiles/".format(file_stuff))
        os.system("docker commit ba0d61f03294 autogrow")
    
    return json_vars, docker_json_vars

def run_autogrow_docker_main(vars):
    """
    This function runs the processing to:
        1) make sure the docker exists and is named autogrow
        2) check that JSON file has basic info
            -receptor, size/center...
        3) copy files to a temp directory
            -receptor, .smi files ...
        4) make a JSON file with modified information for within docker
        5) transfer files to docker container
        6) execute RunAutogrow.py from within the docker containiner
        7) copy files back to the final end dir

    Inputs:
    :param dict argv: Dictionary of User specified variables
    """
    # run 1st part to test the docker
    test_docker_exists()

    # Run parts 2-5
    # 2) check that JSON file has basic info
    #     -receptor, size/center...
    # 3) copy files to a temp directory
    #     -receptor, .smi files ...
    # 4) make a JSON file with modified information for within docker
    # 5) transfer files to docker container
    json_vars, docker_json_vars = handle_json_info(vars)

    # Run part 6) run AutoGrow in the container
    # docker cp foo.txt mycontainer:/foo.txt
    command = "sudo docker run --rm -it -v {}".format(os.path.abspath(".tmp"))
    command = command + "--name autogrow "
    command = command + "/UserFiles/docker_json_vars.json"
    os.system(command)
    print("done")
    sys.exit(0)
    # Move output to cwd
    # shutil.move("./.tmp/autogrow_output.zip", "./autogrow_output.zip")

    # # Delete tmp files
    # shutil.rmtree("./.tmp/")

    # -d `# detached` \
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

# sudo docker run autogrow -it \
#     -v /home/jacob/Documents/autogrow4/Docker/.tmp:/autogrow_work_dir/ \
#     -filename_of_receptor ../tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
#     -output_dir /autogrow_work_dir/autogrow_output/