import os

def config_paths(params: dict):
    _make_paths_abs(params)
    _make_dirs_end_in_sep(params)
    _create_output_folder(params["root_output_folder"])

def _make_paths_abs(params: dict):
    # convert paths to abspath, in case necessary
    for pname in [
        "filename_of_receptor",
        "root_output_folder",
        "source_compound_file",
        "nn1_script",
        "nn2_script",
        "mgltools_directory",
    ]:
        if pname not in list(params.keys()):
            continue
        params[pname] = os.path.abspath(params[pname])

def _make_dirs_end_in_sep(params: dict):
    dir_params = [
        "root_output_folder",
        "mgltools_directory",
    ]

    for dir_param in dir_params:
        if dir_param not in params:
            continue
        if params[dir_param][-1] != os.sep:
            params[dir_param] = params[dir_param] + os.sep


def _create_output_folder(out_folder: str):
    # Check root_output_folder exists
    if os.path.exists(out_folder) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(out_folder)
        except:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )
