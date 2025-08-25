"""
AutoGrow Configuration and Command-Line Argument Parsing Module.

This module provides functionality for setting up and parsing command-line
arguments for the AutoGrow program, an automated drug optimization and
generation tool. It defines the structure of the argument parser, including
various argument groups for different aspects of the program such as general
settings, input/output, receptor information, genetic algorithm parameters,
conversion settings, and scoring options.

The module also supports dynamic addition of argument groups through a plugin
system, allowing for extensibility of the command-line interface.

Key components:
- get_user_params: The main function for parsing and processing user parameters

The module uses argparse for command-line parsing and supports loading
parameters from a JSON file. It also performs parameter validation and
type correction.

Note:
    When adding new parameters or modifying existing ones, make sure to update
    the relevant _add_*_params functions and consider any impacts on the
    plugin system.
"""
import argparse
import copy
import json
from typing import Any, Dict

from autogrow.config import setup_params
from autogrow.config.custom_argparser import CustomArgumentParser, CustomArgumentGroup
from autogrow.config.json_config_utils import (
    convert_json_params_from_unicode,
    save_vars_as_json,
)
from autogrow.config.argument_vars import plugin_arg_groups_to_add
from autogrow.plugins.registry_base import PluginManagerRegistry
from autogrow.validation import validate_all

# parser = argparse.ArgumentParser(
parser = CustomArgumentParser(
    description="AutoGrow: An automated drug optimization and generation tool."
)


def get_user_params() -> Dict[str, Any]:
    """
    Parse command-line arguments and return a dictionary of user parameters.

    This function sets up argument groups, adds arguments to each group, parses
    the command-line arguments, and processes them. It also handles loading
    parameters from a JSON file if specified.

    Returns:
        Dict[str, Any]: A dictionary containing all the parsed and processed
        user parameters.
    """
    global parser

    # Create plugin manager instances but don't register arguments yet
    plugin_managers = PluginManagerRegistry.get_all_plugin_managers()
    
    # Register arguments once
    for plugin_manager in plugin_managers:
        plugin_manager.register_plugin_arguments()

    # General Settings
    general = parser.add_argument_group(
        "General Settings (basic configuration for the program)"
    )
    _add_general_params(general)

    # Input/Output Settings
    io = parser.add_argument_group(
        "Input/Output Settings (directories and files for input and output)"
    )
    _add_io_params(io)

    # Receptor Information
    receptor = parser.add_argument_group(
        "Receptor Information (details about the receptor for docking)"
    )
    _add_receptor_params(receptor)

    # Genetic Algorithm Options
    ga = parser.add_argument_group(
        "Genetic Algorithm Options (settings for the genetic algorithm)"
    )
    _add_ga_params(ga)

    ga_first_gen = parser.add_argument_group(
        "Genetic Algorithm Options Applied to the First Generation (settings for the first generation)"
    )
    _add_ga_first_gen_params(ga_first_gen)
    # Conversion Settings
    conversion = parser.add_argument_group(
        "Conversion Settings (options for file conversion)"
    )
    _add_conversion_params(conversion)

    # Scoring Settings TODO: LATER: Restore this (after refactoring) later.
    # scoring = parser.add_argument_group(
    #     "Scoring Settings (options for scoring docked compounds)"
    # )
    # _add_scoring_params(scoring)

    # Miscellaneous
    # misc = parser.add_argument_group("Miscellaneous (other settings)")
    # _add_misc_params(misc)

    # Now add in plugin arg groups
    for title, arg_vars in plugin_arg_groups_to_add:
        group = parser.add_argument_group(title)
        for arg_var in arg_vars:
            if arg_var.name[:2] != "--":
                arg_var.name = f"--{arg_var.name}"

            if arg_var.action is None:
                group.add_argument(
                    arg_var.name,
                    type=arg_var.type,
                    default=arg_var.default,
                    help=arg_var.help,
                )
            else:
                group.add_argument(
                    arg_var.name,
                    action=arg_var.action,
                    default=arg_var.default,
                    help=arg_var.help,
                )

    args_dict = vars(parser.parse_args())

    if "json" in args_dict and args_dict["json"] is not None:
        # If there's a --json parameter, load all parameters from that file,
        # only using them if they are not specified at the command line.
        new_args_dict = json.load(open(args_dict["json"]))
        new_args_dict = convert_json_params_from_unicode(new_args_dict)
        print(
            "WARNING: Loaded parameters from JSON file will be only used if they are not specified at the command line."
        )
        for k, v in new_args_dict.items():
            if k is not args_dict or args_dict[k] is None:
                args_dict[k] = v

    # copying args_dict so we can delete out of while iterating through the
    # original args_dict
    new_args_dict = copy.deepcopy(args_dict)

    # Remove any None values from the dictionary
    new_args_dict = {k: v for k, v in new_args_dict.items() if v is not None}

    new_args_dict = setup_params(new_args_dict)

    validate_all(new_args_dict)

    # Save variables in vars dict to a .json file for later usage and reference
    # It saves the file to the output_directory + "vars.json"
    # -If AutoGrow has been run multiple times for the same directory it
    # will save the new vars file as append a number to the file name
    # starting with 2. The util scripts will only look at the original "vars.json"
    #     ie) output_directory + "vars_2.json"
    save_vars_as_json(new_args_dict)

    # output the paramters used
    # new_args_dict, printout = load_commandline_parameters(new_args_dict)

    return new_args_dict


def _add_general_params(parser: CustomArgumentGroup):
    """
    Add general parameters to the argument group.

    This function adds arguments for JSON file input, number of processors,
    and multithread mode to the general settings group.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # Allows the run commands to be submitted via a .json file.
    parser.add_argument(
        "--json",
        "-j",
        metavar="param.json",
        help="Path to a JSON file containing parameters. Overrides other arguments.",
    )
    # processors and multithread mode
    parser.add_argument(
        "--procs_per_node",
        "-p",
        type=int,
        metavar="N",
        default=-1,
        help="Number of processors to use for parallel calculations. Set to -1 for all available CPUs.",
    )
    parser.add_argument(
        "--multithread_mode",
        default="multithreading",
        choices=["multithreading", "serial"],
        help="Multithreading mode. 'serial' forces single-processor execution, overriding --procs_per_node.",
    )
    # for postprocessing
    parser.add_argument(
        "--process_input_compounds",
        action="store_true",
        default=False,
        help="Include an analysis of the input compounds in the output, including docking and generating plots for generation 0.",
    )


def _add_io_params(parser: CustomArgumentGroup):
    """
    Add input/output parameters to the argument group.

    This function adds arguments for root output folder and source compound
    file to the input/output settings group.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # Input/Output directories
    parser.add_argument(
        "--output_directory",
        "-o",
        type=str,
        help="Path to the output directory where all results will be saved.",
    )
    parser.add_argument(
        "--source_compound_file",
        "-s",
        type=str,
        help="Path to the source compounds file to seed the first generation. Can be a tab-delineated .smi file or an .sdf file. These molecules will seed the first generation.",
    )


def _add_receptor_params(parser: CustomArgumentGroup):
    """
    Add receptor-related parameters to the argument group.

    This function adds an argument for the receptor file path to the receptor
    information group.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # receptor information
    parser.add_argument(
        "--receptor_path",
        "-r",
        metavar="receptor.pdb",
        help="Path to the receptor file in PDB format.",
    )


def _add_ga_first_gen_params(parser: CustomArgumentGroup):
    """
    Add genetic algorithm parameters for the first generation.

    This function adds arguments specific to the first generation of the
    genetic algorithm, such as number of molecules to seed the next generation,
    diversity molecules, crossovers, mutants, and elitism.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # Seeding next gen and diversity
    parser.add_argument(
        "--top_mols_to_seed_next_generation_first_generation",
        type=int,
        help="Number of molecules from generation 0 to seed generation 1. Overrides --top_mols_to_seed_next_generation for the first generation only. Should be less than the sum of crossovers and mutants for the first generation.",
    )
    parser.add_argument(
        "--diversity_mols_to_seed_first_generation",
        type=int,
        default=10,
        help="Number of diverse molecules to select from generation 0 to seed generation 1. Should be less than the sum of crossovers and mutants for the first generation.",
    )
    parser.add_argument(
        "--number_of_crossovers_first_generation",
        type=int,
        help="Number of new molecules to create via crossover in generation 1. Overrides --number_of_crossovers for the first generation only.",
    )
    parser.add_argument(
        "--number_of_mutants_first_generation",
        type=int,
        help="Number of new molecules to create via mutation in generation 1. Overrides --number_of_mutants for the first generation only.",
    )
    parser.add_argument(
        "--number_elitism_advance_from_previous_gen_first_generation",
        type=int,
        help="Number of elite molecules to carry over from generation 0 to 1, based purely on fitness score. Overrides --number_elitism_advance_from_previous_gen for the first generation only.",
    )
    # parser.add_argument(
    #     "--dock_source_compounds_first",
    #     choices=[True, False, "True", "False", "true", "false"],
    #     default=False,
    #     help="If True source ligands will be docked prior to seeding generation 1. \
    #     If True and the source_compound file already has docking/fitness metric score \
    #     in -2 column of .smi file, it will not redock but reuse the scores from \
    #     the source_compound_file.\
    #     If True and no fitness metric score in -2 column of .smi file, it will \
    #     dock each ligand from the source_compound_file and displayed as generation 0.\
    #     If False, generation 1 will be randomly seeded by the source compounds with \
    #     no preference and there will be no generation 0. \
    #     If performing multiple simulations using same source compounds and protein, \
    #     we recommend running once this and using the generation 0 ranked file as the \
    #     source_compound_file for future simulations. \
    #     Default is True.",
    # )

def _add_ga_params(parser: CustomArgumentGroup):
    """
    Add general genetic algorithm parameters.

    This function adds arguments for the number of generations, elitism,
    redocking elite compounds, and diversity seed depreciation.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # Populations settings
    parser.add_argument(
        "--num_generations",
        type=int,
        default=10,
        help="Total number of generations to run.",
    )
    parser.add_argument(
        "--number_elitism_advance_from_previous_gen",
        type=int,
        default=10,
        help="Number of top-ranked molecules (elites) to carry over to the next generation without modification, based purely on fitness score. Does not select for diversity.",
    )

    parser.add_argument(
        "--diversity_seed_depreciation_per_gen",
        type=int,
        default=2,
        help="Amount by which to decrease the number of diverse seeds per generation.",
    )
    parser.add_argument(
        "--top_mols_to_seed_next_generation",
        type=int,
        default=10,
        help="Number of molecules to select from each generation to seed the next. Should be less than the sum of crossovers and mutants. Can be overridden for the first generation with --top_mols_to_seed_next_generation_first_generation.",
    )
    parser.add_argument(
        "--number_of_crossovers",
        type=int,
        default=10,
        help="Number of new molecules to create via crossover in each generation. Can be overridden for the first generation with --number_of_crossovers_first_generation.",
    )
    parser.add_argument(
        "--number_of_mutants",
        type=int,
        default=10,
        help="Number of new molecules to create via mutation in each generation. Can be overridden for the first generation with --number_of_mutants_first_generation.",
    )


def _add_conversion_params(parser: CustomArgumentGroup):
    """
    Add file conversion parameters.

    This function adds an argument for the path to the OpenBabel executable.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # Path to Open Babel file conversion for docking inputs
    parser.add_argument(
        "--obabel_path",
        help="Path to the Open Babel executable (e.g., '/usr/bin/obabel'). Required for most file conversions. Tip: Run `which obabel` on Linux.",
    )


def _add_scoring_params(parser: CustomArgumentGroup):
    """
    Add scoring-related parameters.

    This function adds arguments for the scoring choice and ligand efficiency
    rescoring option. TODO: Need to implement this in the future.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # scoring
    # parser.add_argument(
    #     "--scoring_choice",
    #     metavar="scoring_choice",
    #     choices=["VINA"],
    #     default="VINA",
    #     help="The scoring_choice to use to assess the ligands docking fitness. \
    #     Default is using Vina/QuickVina2 ligand affinity.",
    # )
    parser.add_argument(
        "--rescore_lig_efficiency",
        action="store_true",
        default=False,
        help="This will divide the final scoring_choice output by the number of \
        non-Hydrogen atoms in the ligand. This adjusted ligand efficiency score will \
        override the scoring_choice value. This is compatible with all scoring_choice options.",
    )


# TODO: Use these parameter later when you implement gypsum option.
# def _add_gypsum_params(parser: CustomArgumentGroup):
#     # gypsum # max variance is the number of conformers made per ligand
#     parser.add_argument(
#         "--max_variants_per_compound",
#         type=int,
#         default=3,
#         help="number of conformers made per ligand. \
#         See Gypsum-DL publication for details",
#     )
#     parser.add_argument(
#         "--gypsum_thoroughness",
#         "-t",
#         type=str,
#         help="How widely Gypsum-DL will search for \
#         low-energy conformers. Larger values increase \
#         run times but can produce better results. \
#         See Gypsum-DL publication for details",
#     )
#     parser.add_argument(
#         "--min_ph",
#         metavar="MIN",
#         type=float,
#         default=6.4,
#         help="Minimum pH to consider.See Gypsum-DL \
#         and Dimorphite-D publication for details.",
#     )
#     parser.add_argument(
#         "--max_ph",
#         metavar="MAX",
#         type=float,
#         default=8.4,
#         help="Maximum pH to consider.See Gypsum-DL \
#         and Dimorphite-D publication for details.",
#     )
#     parser.add_argument(
#         "--pka_precision",
#         metavar="D",
#         type=float,
#         default=1.0,
#         help="Size of pH substructure ranges. See Dimorphite-DL \
#         publication for details.",
#     )
#     parser.add_argument(
#         "--gypsum_timeout_limit",
#         type=float,
#         default=15,
#         help="Maximum time gypsum is allowed to run for a given ligand in seconds. \
#         On average Gypsum-DL takes on several seconds to run for a given ligand, but \
#         factors such as mol size, rotatable bonds, processor speed, and gypsum \
#         settings (ie gypsum_thoroughness or max_variants_per_compound) will change \
#         how long it takes to run. If increasing gypsum settings it is best to increase \
#         the gypsum_timeout_limit. Default gypsum_timeout_limit is 15 seconds",
#     )


def _add_misc_params(parser: CustomArgumentGroup):
    """
    Add miscellaneous parameters.

    This function adds an argument for generating a plot at the end of the run.

    Args:
        parser (CustomArgumentGroup): The argument group to add the
            parameters to.
    """
    # # Make a line plot of the simulation at the end of the run.
    # parser.add_argument(
    #     "--generate_plot",
    #     choices=[True, False, "True", "False", "true", "false"],
    #     default=True,
    #     help="Make a line plot of the simulation at the end of the run.",
    # )
