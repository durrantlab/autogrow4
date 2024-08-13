def config_multiprocessing(params):
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

    # Handle mpi errors if mpi4py isn't installed
    if params["multithread_mode"].lower() == "mpi":
        params["multithread_mode"] = "mpi"
        try:
            import mpi4py
        except:
            printout = (
                "mpi4py not installed but --multithread_mode is set to"
                + " mpi. \n Either install mpi4py or switch "
            )
            printout += "multithread_mode to multithreading or serial"
            raise ImportError(printout)

        try:
            import func_timeout
            from func_timeout import func_timeout, FunctionTimedOut
        except:
            printout = (
                "func_timeout not installed but --multithread_mode is "
                + "set to mpi. \n Either install func_timeout "
            )
            printout += "or switch multithread_mode to"
            printout += " multithreading or serial"
            raise ImportError(printout)

    # Avoid EOF error
    from autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import (
        Parallelizer,
    )

    # launch mpi workers
    if params["multithread_mode"] == "mpi":
        params["parallelizer"] = Parallelizer(
            params["multithread_mode"], params["number_of_processors"]
        )

        if params["parallelizer"] is None:
            printout = "EOF ERRORS FAILED TO CREATE A PARALLIZER OBJECT"
            print(printout)
            raise Exception(printout)

    else:
        params["parallelizer"] = Parallelizer(
            params["multithread_mode"], params["number_of_processors"], True
        )

    return params
