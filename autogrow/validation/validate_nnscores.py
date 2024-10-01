import os
from typing import Any, Dict


def validate_nnscores(params: Dict[str, Any], printout: str) -> str:
    """
    Validates the NNScore scripts and their compatibility with other choices.

    This function checks the existence of NN1 and NN2 scripts, and ensures 
    that the correct docking and conversion choices are made when using NNScores.

    Args:
    :param params: A dictionary containing configuration parameters
    :param printout: A string to accumulate error messages

    Returns:
    :return: The updated printout string

    Raises:
    :raises NotImplementedError: If NN1 or NN2 scripts are not found
    """
    if not os.path.exists(params["nn1_script"]) and not os.path.exists(
        params["nn1_script"].replace('"', "")
    ):
        printout += (
            "\nERROR: Could not find "
            + os.path.basename(params["nn1_script"])
            + " at "
            + params["nn1_script"]
            + "\n"
        )
        raise NotImplementedError(printout)
    if not os.path.exists(params["nn2_script"]) and not os.path.exists(
        params["nn2_script"].replace('"', "")
    ):
        printout += (
            "\nERROR: Could not find "
            + os.path.basename(params["nn2_script"])
            + " at "
            + params["nn2_script"]
            + "\n"
        )
        raise NotImplementedError(printout)

    # CHECK THAT NN1/NN2 are using only traditional Vina Docking
    if params["scoring_choice"] in ["NN1", "NN2"]:
        if params["dock_choice"] != "VinaDocking":
            # printout =
            _nnscore_mgltools_conversion_warning(
                "\nPlease switch dock_choice option to VinaDocking"
            )
        # IF ALTERNATIVE CONVERSION OF PDB2PDBQT CHECK THAT NN1/NN2 are using only MGLTOOLS
        if params["conversion_choice"] != "MGLToolsConversion":
            # printout =
            _nnscore_mgltools_conversion_warning(
                "Please switch conversion_choice option to MGLToolsConversion"
            )

    return printout


def _nnscore_mgltools_conversion_warning(msg: str) -> None:
    """
    Generates a warning message for NNScore and MGLTools conversion issues.

    This function creates a detailed warning message about the requirements
    for using NN1/NN2 scoring methods, and raises an exception with this message.

    Args:
    :param msg: Additional message to be included in the warning

    Raises:
    :raises Exception: Always raises an exception with the generated warning message
    """
    result = (
        "\n\nNeural Networks 1 and 2 (NN1/NN2) are trained on data "
        + "using PDBQT files converted by MGLTools \n"
    )
    result += "and docked using Autodock Vina 1.1.2.\n"
    result += (
        "\nUsing conversion or docking software besides" + " these will not work. \n"
    )
    result += msg + " or deselect NN1/NN2 as the scoring_choice.\n"
    raise Exception(result)
