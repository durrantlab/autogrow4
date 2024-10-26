"""AutoGrow's main execution module.

This module handles the core AutoGrow workflow including population generation,
docking operations, and results visualization. It manages the genetic
algorithm's generational progression and coordinates all major operations.
"""
import __future__

import datetime
import logging
import multiprocessing
import os
import sys
from typing import Any, Dict, Optional

from autogrow import program_info
from autogrow.config.argparser import get_user_params
import autogrow.docking.execute_docking as DockingClass
from autogrow.operators.populate_generation import populate_generation
from autogrow.plugins.plugin_managers import setup_plugin_managers
from autogrow.utils.logging import LogLevel, create_logger, log_info, log_warning


def main(params: Optional[Dict[str, Any]] = None) -> None:
    """Executes the main AutoGrow workflow.

    Orchestrates the genetic algorithm process including population generation,
    molecular docking, and optional plotting. Handles initialization, logging of
    parameters, and manages the generation-by-generation progression of the
    algorithm.

    Args:
        params (Dict[str, Any], optional): Configuration parameters controlling
            the algorithm's behavior. If None, parameters will be parsed from
            command line arguments.

    Raises:
        Exception: If the simulation has already completed the specified number
            of generations.
        ValueError: If population generation fails due to insufficient diversity
            or other constraints.
    """
    start_time = str(datetime.datetime.now())

    multiprocessing.freeze_support()

    create_logger(logging.DEBUG)

    if params is None:
        params = get_user_params()

    # Setup all plugin managers
    setup_plugin_managers(params)

    printout = f"\n(RE)STARTING AUTOGROW 4.0: {str(datetime.datetime.now())}\n"

    printout += program_info()

    printout += "\nUse the -h tag to get detailed help regarding program usage.\n"
    print(printout)
    sys.stdout.flush()

    log_info("Parameters")
    with LogLevel():
        for key in list(params.keys()):
            log_info(f"{key}: {str(params[key])}")

    # Unpack necessary variables
    # output_directory is the root output folder for the run
    num_gens_to_make = params["num_generations"]

    start_gen_num = 1
    if start_gen_num > num_gens_to_make:
        log_warning(
            "This simulation has already been completed to the user defined number \
                of generations. Please check your user variables."
        )
        raise Exception(
            "This simulation has already been completed to the user defined number \
                of generations. Please check your user variables."
        )

    # This is the main loop which will control and execute all commands This
    # is broken into 3 main sections:
    # 1)  operations which populating the new generation with ligands which
    #     both pass the userdefined filter and convert from 1D smiles to 3D
    #     PDB
    # 2)  Docking which handles converting from PDBs to Docking specific
    #     formats and running the actual Docking simulations
    # 3)  Ranking the generation based on the Docking scores
    for gen_num in range(start_gen_num, num_gens_to_make + 1):
        sys.stdout.flush()

        # Get directory for smi to go
        cur_gen_dir = f"{params['output_directory']}generation_{gen_num}{os.sep}"

        log_info(f"Creating Generation {gen_num}")

        with LogLevel():
            smi_new_gen_path, new_gen_predock_cmpds = populate_generation(
                params, gen_num, cur_gen_dir
            )
            sys.stdout.flush()

            if new_gen_predock_cmpds is None:
                raise ValueError(
                    "Population failed to make enough mutants or crossovers... \
                                    Errors could include not enough diversity, too few seeds to the generation, \
                                    the seed mols are unable to cross-over due to lack of similarity,\
                                    or all of the seed lack functional groups for performing reactions."
                )

            # Run file conversions of PDB to docking specific file type and
            # Begin Docking unweighted_ranked_smile_file is the file name
            # where the unweighted ranked but score .smi file resides
            unweighted_ranked_smile_file = DockingClass.run_docking_common(
                params, gen_num, cur_gen_dir, smi_new_gen_path, new_gen_predock_cmpds,
            )

        sys.stdout.flush()

    # if params["generate_plot"] is True:
    #     matplotlib_is_callable = False
    #     try:
    #         import matplotlib  # type: ignore

    #         matplotlib_is_callable = True
    #     except Exception:
    #         matplotlib_is_callable = False
    #     if not matplotlib_is_callable:
    #         print("Can not make figure as matplotlib is not installed")
    #     else:
    #         print("Plotting")
    #         import autogrow.plotting.generate_line_plot as plot

    #         plot.generate_figures(params)

    sys.stdout.flush()

    printout = f"\nAutoGrow4 run started at:   {start_time}\nAutoGrow4 "
    printout = f"{printout}run completed at: {str(datetime.datetime.now())}\n"
    print(printout)

    print("AUTOGROW FINISHED")
