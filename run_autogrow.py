# !/usr/bin/env python

"""This is the executable file for Autogrow 4.0.3. This script should come
first. It should obtain and verify all the parameters work. This than should
pass these parameters variables to the main execution function titled
AutogrowMainExecute.py found in MainFunctions

If you use AutoGrow 4.0.3 in your research, please cite the following reference:
Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm
for de novo drug design and lead optimization. J Cheminform 12, 25 (2020).
[doi: 10.1186/s13321-020-00429-4]
"""


import __future__
import logging
import multiprocessing
import datetime
import sys
from autogrow.config.argparser import get_argparse_vars
from autogrow import program_info
import autogrow.autogrow_main_execute as AutogrowMainExecute
from autogrow.config import load_commandline_parameters
from autogrow.plugins.crossover import CrossoverBase, CrossoverPluginManager
from autogrow.plugins.docking import DockingBase, DockingPluginManager
from autogrow.plugins.mutation import MutationBase, MutationPluginManager
from autogrow.plugins.plugin_manager_base import get_all_plugin_managers
from autogrow.plugins.selectors import SelectorBase, SelectorPluginManager
from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase, SmiTo3DSdfPluginManager
from autogrow.plugins.smiles_filters import SmilesFilterBase, SmilesFilterPluginManager
from autogrow.utils.logging import LogLevel, create_logger, log_info


################
# Run AutoGrow #
################


def _load_plugin_managers() -> None:
    # Note that this just loads the plugins (and sets up args). It doesn't
    # actually create the plugin objects yet.

    # Set up Plugins
    SmilesFilterPluginManager(SmilesFilterBase)
    SelectorPluginManager(SelectorBase)
    DockingPluginManager(DockingBase)
    MutationPluginManager(MutationBase)
    CrossoverPluginManager(CrossoverBase)
    SmiTo3DSdfPluginManager(SmiTo3DSdfBase)


def _setup_plugin_managers(params) -> None:
    # This sets up the plugin managers, after the params have been loaded.

    plugin_managers = get_all_plugin_managers()

    for name, plugin_managers in plugin_managers.items():
        plugin_managers.setup_plugin_manager(params)


def main():
    create_logger(logging.DEBUG)

    _load_plugin_managers()

    args_dict = get_argparse_vars()

    # Setup all plugin managers
    _setup_plugin_managers(args_dict)

    _run_autogrow_with_params(args_dict)


# TODO: Rename this here and in `main`
def _run_autogrow_with_params(args_dict):
    start_time = str(datetime.datetime.now())
    # load the commandline parameters

    printout = f"\n(RE)STARTING AUTOGROW 4.0: {str(datetime.datetime.now())}\n"

    printout += program_info()

    printout += "\nUse the -h tag to get detailed help regarding program usage.\n"
    print(printout)
    sys.stdout.flush()

    # output the paramters used
    params, printout = load_commandline_parameters(args_dict)

    # print out the UserVars for the record
    # print("\n=====================================================")
    # print("==============   Parameters as list:  ===============")
    log_info("Parameters")
    with LogLevel():
        for key in list(params.keys()):
            log_info(f"{key}: {str(params[key])}")

    # print("\n=====================================================")
    # print("===========   Parameters as dictionary:  ============")
    # print(params)
    # print("=====================================================")
    # print("=====================================================\n\n")

    AutogrowMainExecute.main_execute(params)

    # Print completion message

    printout = f"\nAutoGrow4 run started at:   {start_time}\nAutoGrow4 "
    printout = f"{printout}run completed at: {str(datetime.datetime.now())}\n"
    print(printout)

    print("AUTOGROW FINISHED")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
