import os
from typing import Any, Dict


def config_mgltools(params: Dict[str, Any]) -> None:
    """
    Configures MGLTools paths in the user parameters.

    If the conversion choice is set to "MGLToolsConversion", this function ensures that the
    necessary MGLTools-related scripts are correctly specified in the parameters. It sets
    default paths for 'prepare_ligand4.py', 'prepare_receptor4.py', and 'mgl_python' based
    on the 'mgltools_directory' if they are not already provided.

    Inputs:
    :param params: Dictionary of user variables.

    Returns:
    :returns: None. The function modifies the 'params' dictionary in place.
    """
    if params["conversion_choice"] != "MGLToolsConversion":
        # Not set to use mgltools
        return

    # find other mgltools-related scripts
    if params.get("prepare_ligand4.py", "") == "":
        params["prepare_ligand4.py"] = (
            params["mgltools_directory"]
            + "MGLToolsPckgs"
            + os.sep
            + "AutoDockTools"
            + os.sep
            + "Utilities24"
            + os.sep
            + "prepare_ligand4.py"
        )
    if params.get("prepare_receptor4.py", "") == "":
        params["prepare_receptor4.py"] = (
            params["mgltools_directory"]
            + "MGLToolsPckgs"
            + os.sep
            + "AutoDockTools"
            + os.sep
            + "Utilities24"
            + os.sep
            + "prepare_receptor4.py"
        )
    if params.get("mgl_python", "") == "":
        params["mgl_python"] = (
            params["mgltools_directory"] + "bin" + os.sep + "pythonsh"
        )
