import os
import sys
from autogrow.config.config_utils import copy_new_custom_py_file, get_path_to_custom_script, make_complete_children_dict
from autogrow.validation.validate_custom_classes import validate_custom_param_type

def setup_custom_dock_and_conversion_scoring_options(params: dict) -> dict:
    """
    This function handles selecting the user defined Custom options
    for Custom docking Conversion, and scoring scripts.

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: dict params: Dictionary of User variables with the added options
    """
    if (
        params["conversion_choice"] != "Custom"
        and params["dock_choice"] != "Custom"
        and params["scoring_choice"] != "Custom"
    ):
        # It's not a custom approach.
        return params
    
    master_need_restart = False
    master_printout = ""
    if params["conversion_choice"] == "Custom":
        params, need_restart, printout = _setup_custom_conversion_script(params)
        if need_restart is True:
            master_need_restart = True
            master_printout += printout
    if params["dock_choice"] == "Custom":
        params, need_restart, printout = _setup_custom_docking_script(params)
        if need_restart is True:
            master_need_restart = True
            master_printout = master_printout + printout
    if params["scoring_choice"] == "Custom":
        params, need_restart, printout = _setup_custom_scoring_script(params)
        if need_restart is True:
            master_need_restart = True
            master_printout = master_printout + printout

    if master_need_restart:
        print(master_printout)
        # Technically Exit intentionally but maybe should be a raise Exception
        sys.exit(0)

    return params


def _setup_custom_scoring_script(params):
    """
    This will handle Custom scoring_scripts

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: dict params: Dictionary of User variables modified with the
        params["dock_choice"] set to the new custom dock_choice
    :returns: bool need_restart: If True AutoGrow will need to be restarted
        after all other files are incorporated
    :returns: str printout: "" or a message to be print prior to
        being restarted if needed
    """
    need_restart = False
    printout = ""
    if params["custom_scoring_script"] is None:
        return params, need_restart, printout

    validate_custom_param_type(
        "custom_scoring_script",
        params,
        list,
        str,
        "list of [name_Scoring_script1, Path/to/name_Scoring_script1.py]",
    )

    full_children_dict = make_complete_children_dict("ParentScoring")
    custom_class = params["custom_scoring_script"]
    if custom_class[0] not in full_children_dict.keys():
        new_file = get_path_to_custom_script(
            custom_class,
            "custom_scoring_script",
            "list of [name_scoring_script1, Path/to/name_scoring_script1.py]",
            ["docking", "scoring", "scoring_classes", "scoring_functions"],
        )

        # Add copy the script to the scoring_choices folder
        copy_new_custom_py_file(custom_class[1], new_file)

        need_restart = True

    params["scoring_choice"] = custom_class[0]
    return params, need_restart, printout



def _setup_custom_conversion_script(params):
    """
    This will handle Custom Conversion_scripts

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: dict params: Dictionary of User variables modified with
        the params["conversion_choice"] set to the new custom conversion_choice
    :returns: bool need_restart: If True AutoGrow will need to be restarted
        after all other files are incorporated
    :returns: str printout: "" or a message to be print prior to being
        restarted if needed
    """
    need_restart = False
    printout = ""

    if params["custom_conversion_script"] is None:
        return params, need_restart, printout

    validate_custom_param_type(
        "custom_conversion_script",
        params,
        list,
        str,
        "list of [name_Conversion_script1, Path/to/name_Conversion_script1.py]",
    )

    full_children_dict = make_complete_children_dict("parent_pdbqt_converter")
    custom_class = params["custom_conversion_script"]
    if custom_class[0] not in full_children_dict.keys():

        new_file = get_path_to_custom_script(
            custom_class,
            "custom_conversion_script",
            "list of [name_Conversion_script1, Path/to/name_Conversion_script1.py]",
            ["docking", "docking_class", "docking_file_conversion"],
        )

        # Add copy the script to the docking_file_conversion folder
        copy_new_custom_py_file(custom_class[1], new_file)

        need_restart = True

    params["conversion_choice"] = custom_class[0]

    return params, need_restart, printout


def _setup_custom_docking_script(params):
    """
    This will handle Custom Docking_scripts

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: dict params: Dictionary of User variables modified with the
        params["dock_choice"] set to the new custom dock_choice
    :returns: bool need_restart: If True AutoGrow will need to be estarted
         after all other files are incorporated
    :returns: str printout: "" or a message to be print prior to being
        restarted if needed
    """
    need_restart = False
    printout = ""

    if params["custom_docking_script"] is None:
        return params, need_restart, printout

    validate_custom_param_type(
        "custom_docking_script",
        params,
        list,
        str,
        "list of [name_Docking_script1, Path/to/name_Docking_script1.py]",
    )

    full_children_dict = make_complete_children_dict("ParentDocking")
    custom_class = params["custom_docking_script"]
    if custom_class[0] not in full_children_dict.keys():
        new_file = get_path_to_custom_script(
            custom_class,
            "custom_docking_script",
            "list of [name_Docking_script1, Path/to/name_Docking_script1.py]",
            ["docking", "docking_class", "docking_class_children"],
        )

        # Add copy the script to the children folder
        copy_new_custom_py_file(custom_class[1], new_file)

        need_restart = True

    params["dock_choice"] = custom_class[0]
    return params, need_restart, printout

