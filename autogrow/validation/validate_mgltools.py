import os

def validate_mgltools(params: dict, printout: str):
    if params["conversion_choice"] != "MGLToolsConversion":
        # mgl tools not used, no need to validate
        return

    # If MGLTools is being used handle its paths
    if "mgltools_directory" not in params.keys():
        printout = (
            "\nmgltools_directory was not provided but conversion_choice"
            + " is set to MGLToolsConversion. Please"
            + " provide the path to the mgltools_directory\n"
        )
        raise NotImplementedError(printout)
    if not os.path.exists(params["mgltools_directory"]):
        raise NotImplementedError("mgltools_directory does not exist")
    if not os.path.isdir(params["mgltools_directory"]):
        raise NotImplementedError(
            "mgltools_directory is not a directory. Check your input parameters."
        )

    # Make sure scripts and executables exist
    # If MGLTools is being used handle its paths
    if not os.path.exists(params["prepare_ligand4.py"]) and not os.path.exists(
        params["prepare_ligand4.py"].replace('"', "")
    ):
        printout += (
            "\nERROR: Could not find prepare_ligand4.py at "
            + params["prepare_ligand4.py"]
            + "\n"
        )
        raise NotImplementedError(printout)
    if not os.path.exists(params["prepare_receptor4.py"]) and not os.path.exists(
        params["prepare_receptor4.py"].replace('"', "")
    ):
        printout += (
            "\nERROR: Could not find prepare_receptor4.py at "
            + params["prepare_receptor4.py"]
            + "\n"
        )
        raise NotImplementedError(printout)
    if not os.path.exists(params["mgl_python"]) and not os.path.exists(
        params["mgl_python"].replace('"', "")
    ):
        printout += (
            "\nERROR: Could not find pythonsh at " + params["mgl_python"] + "\n"
        )
        raise NotImplementedError(printout)
