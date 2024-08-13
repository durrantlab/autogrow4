import os

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

    folder_name_path = f"{root_folder_path}Run_"
    print(folder_name_path)

    last_run_number = _find_previous_runs(folder_name_path)

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
            os.makedirs(folder_path)

    print("The Run number is: ", run_number)
    print("The Run folder path is: ", folder_path)
    print("")
    return folder_path

def _find_previous_runs(folder_name_path):
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
    return i - 1  # last_run_number


