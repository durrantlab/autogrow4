"""
Configuration system for managing AutoGrow settings and parameters.

This package handles all aspects of configuration for AutoGrow, including:

- Parameter parsing and validation
- Path configuration and management
- Multiprocessing setup
- Default value definitions
- JSON configuration utilities
"""

import contextlib
from typing import Any, Dict
import os

from autogrow.config.config_multiprocessing import config_multiprocessing
from autogrow.config.config_paths import config_paths
from autogrow.utils.logging import log_info


def setup_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Set up parameters, correct types, and set defaults.

    This function processes the original parameters, casts certain parameters to
    the correct types, sets up filters, configures paths and multiprocessing,
    and applies default values where necessary.

    Args:
        orig_params (Dict[str, Any]): Dictionary of original parameters.

    Returns:
        Dict[str, Any]: Dictionary of corrected and completed parameters.
    """
    _cast_some_params(params)
    _set_missing_first_generation_params(params)
    config_paths(params)
    config_multiprocessing(params)

    # Check if the user wants to continue a run or start a new run. Make new run
    # directory if necessary. return the Run folder path The run folder path
    # will be where we place our generations and output files
    if not os.path.exists(params["output_directory"]):
        os.makedirs(params["output_directory"])
        log_info(f"Making the output folder path: {params['output_directory']}")

    return params


def _cast_some_params(input_params: Dict[str, Any]) -> None:
    """
    Cast specific parameters to their required types.

    This function ensures that certain parameters, particularly those related to
    dimensions and docking settings, are cast to the correct data types.

    Args:
        input_params (Dict[str, Any]): The parameters dictionary to be modified.

    Raises:
        Exception: If a required parameter cannot be cast to the correct type.
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
    """
    Set missing parameters specific to the first generation.

    This function checks for parameters specific to the first generation and
    sets them if they're not defined. If the first-generation version is not
    defined, it uses either the default value of 10 or the value defined for
    subsequent generations.

    Args:
        params (Dict[str, Any]): The parameters dictionary to be modified.
    """
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


def _convert_param_to_int_if_needed(
    param_name: str, input_params: Dict[str, Any]
) -> None:
    """
    Convert a parameter to an integer if needed.

    This function attempts to convert the specified parameter to an integer. If
    the parameter is "None" or None, it leaves it as None. If conversion fails,
    it raises an exception.

    Args:
        param_name (str): Name of the parameter to convert.
        input_params (Dict[str, Any]): Dictionary of input parameters.

    Raises:
        Exception: If the parameter cannot be converted to an integer and is not
            None.
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
