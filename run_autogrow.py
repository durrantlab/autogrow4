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
import multiprocessing
import datetime
import sys
from autogrow.config.argparser import get_argparse_vars
from autogrow import program_info
import autogrow.autogrow_main_execute as AutogrowMainExecute
from autogrow.config import load_commandline_parameters


################
# Run AutoGrow #
################


def main():

    args_dict = get_argparse_vars()

    _run_autogrow_with_params(args_dict)


# TODO Rename this here and in `main`
def _run_autogrow_with_params(args_dict):
    start_time = str(datetime.datetime.now())
    # load the commandline parameters

    printout = f"(RE)STARTING AUTOGROW 4.0: {str(datetime.datetime.now())}"

    printout += program_info()

    printout += "\nUse the -h tag to get detailed help regarding program usage.\n"
    print(printout)
    sys.stdout.flush()

    # output the paramters used
    printout += "\nPARAMETERS" + "\n"
    printout += " ========== " + "\n"

    params, printout = load_commandline_parameters(args_dict)

    # print out the UserVars for the record
    print("\n=====================================================")
    print("==============   Parameters as list:  ===============")
    for key in list(params.keys()):
        print(key, params[key])
    print("\n=====================================================")
    print("===========   Parameters as dictionary:  ============")
    print(params)
    print("=====================================================")
    print("=====================================================\n\n")

    AutogrowMainExecute.main_execute(params)

    # Print completion message

    printout = f"\nAutoGrow4 run started at:   {start_time}\nAutoGrow4 "
    printout = f"{printout}run completed at: {str(datetime.datetime.now())}\n"
    print(printout)

    print("AUTOGROW FINISHED")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()