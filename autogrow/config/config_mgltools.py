import os

def config_mgltools(params: dict):
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
