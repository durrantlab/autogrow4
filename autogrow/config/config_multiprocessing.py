"""
Multiprocessing Configuration Module

This module handles the configuration of multiprocessing settings for AutoGrow.
"""
from typing import Any, Dict


def config_multiprocessing(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Configure multiprocessing settings and establish a Parallelizer object.

    This function handles the multiprocessing configuration based on user
    parameters. It creates a Parallelizer object and adds it to the params
    dictionary. If the multithread mode is set to "serial", it ensures that
    only one processor is used.

    Args:
        params (Dict[str, Any]): Dictionary of user variables governing how the
            program runs.

    Returns:
        Dict[str, Any]: Updated dictionary of user variables with the
            Parallelizer object added.

    Notes:
        If multithread_mode is set to "serial", the procs_per_node is
        forcibly set to 1, regardless of the user's input.
    """
    # Handle Serial overriding procs_per_node
    # serial fixes it to 1 processor
    if params["multithread_mode"].lower() == "serial":
        params["multithread_mode"] = "serial"
        if params["procs_per_node"] != 1:
            print(
                "Because --multithread_mode was set to serial, "
                + "this will be run on a single processor."
            )
        params["procs_per_node"] = 1

    # Avoid EOF error
    from autogrow.utils.parallelizer import Parallelizer

    params["parallelizer"] = Parallelizer(
        params["multithread_mode"], params["procs_per_node"]
    )

    return params
