import contextlib
from typing import Any, Dict

from autogrow.config.config_filters import setup_filters
from autogrow.config.config_multiprocessing import config_multiprocessing
from autogrow.config.config_paths import config_paths
from autogrow.config.config_run_directory import set_run_directory
from autogrow.config.defaults import define_defaults



def setup_params(orig_params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Set up parameters, correct types, and set defaults.

    Inputs:
    :param orig_params: Dictionary of original parameters.

    Returns:
    :returns: corrected_params: Dictionary of corrected parameters.
    """
    _cast_some_params(orig_params)
    _set_missing_first_generation_params(orig_params)
    config_paths(orig_params)
    config_multiprocessing(orig_params)

    # Start with getting the default values
    default_params = define_defaults()

    orig_params = setup_filters(orig_params)

    # Check if the user wants to continue a run or start a new run. Make new run
    # directory if necessary. return the Run folder path The run folder path
    # will be where we place our generations and output files
    orig_params["output_directory"] = set_run_directory(
        orig_params["root_output_folder"]
    )
    corrected_params = _correct_param_to_default_types(orig_params, default_params)

    # Now add defaults to corrected_params
    for key in list(default_params.keys()):
        if key not in list(corrected_params.keys()):
            corrected_params[key] = default_params[key]

    return corrected_params


def _correct_param_to_default_types(
    params: Dict[str, Any], default_params_for_ref: Dict[str, Any]
) -> Dict[str, Any]:
    """
    This checks that all the user variables loaded in use that same or comparable
    datatypes as the defaults in params. This prevents type issues later in the
    simulation.

    Given the many uservars and the possibility for intentional differences,
    especially as the program is developed, this function tries to be
    NOT OPINIONATED, only correcting for several obvious and easy to correct issues
    of type discrepancies occur between params[key] and default_params_for_ref[key]
        ie
            1) params[key] = "true" and default_params_for_ref[key] = False
                this script will not change params[key] to False... it will
                convert "true" to True
                ---> params[key]=True
            2) params[key] = "1.01" and default_params_for_ref[key] = 2.1
                this script will change params[key] from "1.01" to float(1.01)

    Inputs:
    :param dict params: Dictionary of user specified variables, to
        correct.
    :param dict default_params_for_ref: Dictionary of program defaults, used to 
        identify the proper types.
    
    Returns:
    :returns: dict params: Dictionary of corrected user specified variables
    """
    for key in list(params.keys()):
        if key not in list(default_params_for_ref.keys()):
            # Examples may be things like receptor_path or dimensions of
            # the docking box. Just skip these
            continue

        if type(params[key]) == type(default_params_for_ref[key]):
            # The types are the same, so you can go on to the text one.
            continue

        # Several variable default is None which means checks are processed
        # elsewhere...
        if default_params_for_ref[key] is None:
            # check argv[key] is "none" or "None"
            if type(params[key]) != str:
                continue

            if params[key].lower() == "none":
                params[key] = None
        elif type(default_params_for_ref[key]) in [int, float]:
            if type(params[key]) in [int, float]:
                # this is fine
                continue

            if type(params[key]) == str:
                try:
                    temp_item = float(params[key])
                    if type(temp_item) == float:
                        params[key] = temp_item
                    else:
                        _wrong_type_error(key, params, default_params_for_ref)
                except Exception:
                    _wrong_type_error(key, params, default_params_for_ref)
            else:
                import pdb; pdb.set_trace()
                _wrong_type_error(key, params, default_params_for_ref)
        elif type(default_params_for_ref[key]) == bool:
            if params[key] is None:
                # Do not try to handle this. May make sense.
                continue
            if type(params[key]) == str:
                if params[key].lower() in ["true", "1"]:
                    params[key] = True
                elif params[key].lower() in ["false", "0"]:
                    params[key] = False
                elif params[key].lower() in ["none"]:
                    params[key] = None
                else:
                    _wrong_type_error(key, params, default_params_for_ref)
            else:
                _wrong_type_error(key, params, default_params_for_ref)
    return params


def _cast_some_params(input_params: Dict[str, Any]) -> None:
    """
    Some parameters must be cast to different types.

    Required Variables go here.

    Inputs:
    :param dict input_params: The parameters. A dictionary of {parameter name: value}.
    """

    # Make sure the dimmensions are in floats. If in int convert to float.
    for x in ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]:
        if type(input_params[x]) == float:
            continue
        if type(input_params[x]) == int:
            input_params[x] = float(input_params[x])
        else:
            printout = f"\n{x} must be a float value.\n"
            print(printout)
            raise Exception(printout)

    # Check Docking exhaustiveness and modes...
    _convert_param_to_int_if_needed("docking_exhaustiveness", input_params)
    _convert_param_to_int_if_needed("docking_num_modes", input_params)


def _set_missing_first_generation_params(params: Dict[str, Any]) -> None:
    # Check parameters specific to the first generation. If not defined, use the
    # default of 10. If defined, use the same number for the first generation.

    for pname in [
        "top_mols_to_seed_next_generation",
        "number_of_crossovers",
        "number_of_mutants",
        "number_elitism_advance_from_previous_gen",
    ]:
        if f"{pname}_first_generation" not in list(params.keys()):
            # The first-generation version is not defined.
            if pname not in list(params.keys()):
                # Subsequent-generation versions are also not defined. Use
                # defined default of 10.
                params[pname] = 10
                params[f"{pname}_first_generation"] = 10
            else:
                # Subsequent-generation versions are defined. Use the same
                # number for the first generation.
                params[f"{pname}_first_generation"] = params[pname]


def _wrong_type_error(key: str, argv: Dict[str, Any], params: Dict[str, Any]) -> None:
    """
    This function is used to raise an error when the user has inputted the wrong
    type for a variable. It is used in the check_value_types function.

    Inputs:
    :param str key: the key of the variable that is the wrong type
    :param dict argv: the dictionary of user specified variables
    :param dict params: the dictionary of program defaults
    """
    printout = (
        "This parameter is the wrong type. \n \t Check :"
        + f" {key} type={type(argv[key])}\n"
    )
    printout += f"\t Should be type={type(params[key])}\n\t"
    printout += "Please check Autogrow documentation using -h"
    raise IOError(printout)


def _convert_param_to_int_if_needed(
    param_name: str, input_params: Dict[str, Any]
) -> None:
    """
    Converts the parameter to an integer if needed.

    If the parameter is a string representation of an integer or "None", it converts it accordingly.

    Inputs:
    :param param_name: Name of the parameter to convert.
    :param input_params: Dictionary of input parameters.
    """
    if param_name not in list(input_params.keys()):
        return

    # Make sure "None" is same as None
    if input_params[param_name] == "None":
        input_params[param_name] = None

    if input_params[param_name] is None:
        return

    with contextlib.suppress(Exception):
        input_params[param_name] = int(input_params[param_name])

    if type(input_params[param_name]) in [
        int,
        float,
    ]:
        # Conversion successful
        return

    raise Exception(
        f"{param_name} needs to be an integer. \
        If you do not know what to use, leave this blank and the \
        default for the docking software will be used."
    )
