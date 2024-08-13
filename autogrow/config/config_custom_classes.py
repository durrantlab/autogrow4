import sys

from autogrow.config.config_utils import copy_new_custom_py_file, get_path_to_custom_script, make_complete_children_dict
from autogrow.validation.validate_custom_classes import validate_custom_param_type

############################################
########   Custom Option Settings   ########
############################################


def handle_custom_params_if_argparsed(input_params):
    """
    There are several Custom options such as filters, docking software which
    are lists because a single filter can use multiple options at once. This
    function is used to properly import and parse those user variables if using
    the commandline argparse. If JSON is used, these will already be lists.

    This function will handle those if there are used and return
    the modified input_params dict

    Inputs:
    :param dict input_params: The parameters. A dictionary of
        {parameter name: value}.
    Returns:
    :returns: dict input_params: The parameters. A dictionary of
        {parameter name: value}.
    """

    for pname in [
        "alternative_filter",  # Custom Filters # TODO: Why not called custom_filter for consistency?
        "custom_conversion_script",
        "custom_docking_script",
        "custom_scoring_script",
    ]:
        if pname not in input_params.keys():
            input_params[pname] = None
        if input_params[pname] not in [None, [], "", "[]"]:
            _parse_custom_filter_input(input_params, pname)

    return input_params


def _parse_custom_filter_input(input_params, arg1):
    orginal = input_params[arg1][0]
    orginal = orginal.replace("[[", "[").replace("]]", "]")
    new_alternative_filter = []
    for custom_filter in orginal.split("]"):
        custom_filter = custom_filter.replace("[", "").replace("]", "")
        custom_filter = [x for x in custom_filter.split(",") if x != ""]
        if len(custom_filter) == 2:
            new_alternative_filter.append(custom_filter)
    input_params[arg1] = new_alternative_filter



