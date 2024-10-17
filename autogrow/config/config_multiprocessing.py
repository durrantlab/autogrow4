from typing import Any, Dict


def config_multiprocessing(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    This function handles the multiprocessing functions. It establishes a Paralellizer object
    and adds it to the params dictionary.
    Inputs:
    :param dict params: dict of user variables which will govern how the programs runs
    Returns:
    :returns: dict params: dict of user variables which will govern how the programs runs
    """

    # Handle Serial overriding number_of_processors
    # serial fixes it to 1 processor
    if params["multithread_mode"].lower() == "serial":
        params["multithread_mode"] = "serial"
        if params["number_of_processors"] != 1:
            print(
                "Because --multithread_mode was set to serial, "
                + "this will be run on a single processor."
            )
        params["number_of_processors"] = 1

    # Avoid EOF error
    from autogrow.utils.parallelizer import Parallelizer

    params["parallelizer"] = Parallelizer(
        params["multithread_mode"], params["number_of_processors"]
    )

    return params
