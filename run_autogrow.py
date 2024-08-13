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

# NOTE: imports of files are burried below to prevent EOF issues in MPI mode

################
# Run AutoGrow #
################

def main():
    from autogrow.config.argparser import get_argparse_vars
    args_dict = get_argparse_vars()

    if args_dict["cache_prerun"] is False:
        start_time = str(datetime.datetime.now())
        # load the commandline parameters

        printout = f"(RE)STARTING AUTOGROW 4.0: {str(datetime.datetime.now())}"

        from autogrow import program_info
        printout += program_info()

        printout += "\nUse the -h tag to get detailed help regarding program usage.\n"
        print(printout)
        sys.stdout.flush()

        # output the paramters used
        printout += "\nPARAMETERS" + "\n"
        printout += " ========== " + "\n"

        from autogrow.config import load_commandline_parameters
        user_vars, printout = load_commandline_parameters(args_dict)

        # print out the UserVars for the record
        print("\n=====================================================")
        print("==============   Parameters as list:  ===============")
        for key in list(user_vars.keys()):
            print(key, user_vars[key])
        print("\n=====================================================")
        print("===========   Parameters as dictionary:  ============")
        print(user_vars)
        print("=====================================================")
        print("=====================================================\n\n")

        # Run AUTOGROW. Import move here to prevent EOF in MPI mode. importing
        # files before the Parallelizer class is established in MPI mode can
        # produce errors
        import autogrow.autogrow_main_execute as AutogrowMainExecute

        AutogrowMainExecute.main_execute(user_vars)

        # Print completion message

        printout = f"\nAutoGrow4 run started at:   {start_time}\nAutoGrow4 "
        printout = f"{printout}run completed at: {str(datetime.datetime.now())}\n"
        print(printout)

        print("AUTOGROW FINISHED")

        # # kill mpi workers
        user_vars["parallelizer"].end(user_vars["multithread_mode"])

    else:  # cache prerun. This is necessary to prevent race conditions in mpi mode.
        import autogrow.user_vars
        import autogrow.autogrow_main_execute as AutogrowMainExecute
        import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
