import os

def config_params_under_windows(params: dict, printout: str) -> str:
    if os.name not in {"nt", "ce"}:
        # It's not windows
        return ""

    # so it's running under windows. multiprocessing disabled
    params["number_of_processors"] = 1
    printout += "\nWARNING: Multiprocessing is disabled on windows machines.\n"

    # convert path names with spaces if this is windows
    params_that_are_paths = [
        "filename_of_receptor",
        "root_output_folder",
        "nn1_script",
        "nn2_script",
        "mgltools_directory",
        "prepare_ligand4.py",
        "prepare_receptor4.py",
        "mgl_python",
    ]

    for pname in params_that_are_paths:
        if pname in params and " " in params[pname]:
            params[pname] = f'"{params[pname]}"'

    return printout