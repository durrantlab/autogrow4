"""
AutoGrow Configuration and Command-Line Argument Parsing Module

This module provides functionality for setting up and parsing command-line
arguments for the AutoGrow program, an automated drug optimization and
generation tool. It defines the structure of the argument parser, including
various argument groups for different aspects of the program such as general
settings, input/output, receptor information, genetic algorithm parameters,
conversion settings, and scoring options.

The module also supports dynamic addition of argument groups through a plugin
system, allowing for extensibility of the command-line interface.

Key components:
- ArgumentVars: A dataclass for defining argument properties
- register_argparse_group: A function for plugins to register new argument groups
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
from dataclasses import dataclass
import json
from typing import Any, Dict, List, Optional, Tuple, Type

from autogrow.config import setup_params
from autogrow.config.json_config_utils import (
    convert_json_params_from_unicode,
    save_vars_as_json,
)
from autogrow.utils.logging import log_warning
from autogrow.validation import validate_all

parser = argparse.ArgumentParser(
    description="AutoGrow: An automated drug optimization and generation tool."
)

plugin_arg_groups_to_add = []


@dataclass
class ArgumentVars:
    """
    A data class to represent argument variables for command-line parsing.

    This class is used to define the properties of command-line arguments,
    which can be used to dynamically add arguments to the ArgumentParser.

    Attributes:
        name (str): The name of the argument.
        default (Any): The default value of the argument.
        help (str): The help text describing the argument.
        type (Optional[Type]): The expected type of the argument value.
            Defaults to None.
        action (Optional[str]): The action to be taken when the argument is
            encountered. Defaults to None.
    """

    name: str
    default: Any
    help: str
    type: Optional[Type] = None
    action: Optional[str] = None


def register_argparse_group(title: str, arg_vars: List[ArgumentVars]):
    """
    Register an argument group with its associated arguments.

    This function is used by plugins to add their own argument groups and
    arguments to the main parser.

    Args:
        title (str): The title of the argument group.
        arg_vars (List[ArgumentVars]): A list of ArgumentVars objects
            representing the arguments to be added to the group.
    """
    global plugin_arg_groups_to_add
    plugin_arg_groups_to_add.append((title, arg_vars))


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
    global plugin_arg_groups_to_add

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
        "General Genetic Algorithm Options (settings for the genetic algorithm)"
    )
    _add_ga_params(ga)

    ga_first_gen = parser.add_argument_group(
        "Genetic Algorithm Options Applied to the First Generation (settings for the first generation)"
    )
    _add_ga_first_gen_params(ga_first_gen)

    ga_subsequent_gen = parser.add_argument_group(
        "Genetic Algorithm Options Applied to Subsequent Generations (settings for all generations after the first)"
    )
    _add_ga_subsequent_gen_params(ga_subsequent_gen)

    # Conversion Settings
    conversion = parser.add_argument_group(
        "Conversion Settings (options for file conversion)"
    )
    _add_conversion_params(conversion)

    # Scoring Settings TODO: Restore this (after refactoring) later.
    # scoring = parser.add_argument_group(
    #     "Scoring Settings (options for scoring docked compounds)"
    # )
    # _add_scoring_params(scoring)

    # Miscellaneous
    misc = parser.add_argument_group("Miscellaneous (other settings)")
    _add_misc_params(misc)

    # Now add in plugin arg groups
    titles = []
    titles_to_arg_vars = {}
    for title, arg_vars in plugin_arg_groups_to_add:
        if title not in titles:
            titles.append(title)

        if title not in titles_to_arg_vars:
            titles_to_arg_vars[title] = []

        titles_to_arg_vars[title].extend(arg_vars)

    for title in titles:
        arg_vars = titles_to_arg_vars[title]
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
        # ignoring the command line.
        new_args_dict = json.load(open(args_dict["json"]))
        new_args_dict = convert_json_params_from_unicode(new_args_dict)
        print(
            "WARNING: Loaded parameters from JSON file. Ignoring other parameters specified at the command line."
        )
    else:
        # No --json, so process using command line

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


def _add_general_params(parser: argparse._ArgumentGroup):
    """
    Add general parameters to the argument group.

    This function adds arguments for JSON file input, number of processors,
    and multithread mode to the general settings group.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # Allows the run commands to be submitted via a .json file.
    parser.add_argument(
        "--json",
        "-j",
        metavar="param.json",
        help="Name of a json file containing all parameters. \
        Overrides other arguments.",
    )

    # processors and multithread mode
    parser.add_argument(
        "--number_of_processors",
        "-p",
        type=int,
        metavar="N",
        default=1,
        help="Number of processors to use for parallel calculations. Set to -1 for all available CPUs.",
    )
    parser.add_argument(
        "--multithread_mode",
        default="multithreading",
        choices=["multithreading", "serial"],
        help="Determine what style \
        multithreading: multithreading or serial. serial will override \
        number_of_processors and force it to be on a single processor.",
    )


def _add_io_params(parser: argparse._ArgumentGroup):
    """
    Add input/output parameters to the argument group.

    This function adds arguments for root output folder and source compound
    file to the input/output settings group.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # Input/Output directories
    parser.add_argument(
        "--root_output_folder",
        "-o",
        type=str,
        help="The Path to the folder which all output files will be placed.",
    )
    parser.add_argument(
        "--source_compound_file",
        "-s",
        type=str,
        help="PATH to the file containing the source compounds. It must be \
        tab-delineated .smi file. These ligands will seed the first generation.",
    )


def _add_receptor_params(parser: argparse._ArgumentGroup):
    """
    Add receptor-related parameters to the argument group.

    This function adds an argument for the receptor file path to the receptor
    information group.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # receptor information
    parser.add_argument(
        "--receptor_path",
        "-r",
        metavar="receptor.pdb",
        help="The path to the receptor file. Should be .pdb file.",
    )


def _add_ga_first_gen_params(parser: argparse._ArgumentGroup):
    """
    Add genetic algorithm parameters for the first generation.

    This function adds arguments specific to the first generation of the
    genetic algorithm, such as number of molecules to seed the next generation,
    diversity molecules, crossovers, mutants, and elitism.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # Seeding next gen and diversity
    parser.add_argument(
        "--top_mols_to_seed_next_generation_first_generation",
        type=int,
        help="Number of mols that seed next generation, for the first generation.\
        Should be less than number_of_crossovers_first_generation + number_of_mutations_first_generation\
        If not defined it will default to top_mols_to_seed_next_generation",
    )
    parser.add_argument(
        "--diversity_mols_to_seed_first_generation",
        type=int,
        default=10,
        help="Should be less than number_of_crossovers_first_generation \
        + number_of_mutations_first_generation",
    )
    parser.add_argument(
        "--number_of_crossovers_first_generation",
        type=int,
        help="The number of ligands which will be created via crossovers in the \
        first generation. If not defined it will default to number_of_crossovers",
    )
    parser.add_argument(
        "--number_of_mutants_first_generation",
        type=int,
        help="The number of ligands which will be created via mutation in \
        the first generation. If not defined it will default to number_of_mutants",
    )
    parser.add_argument(
        "--number_elitism_advance_from_previous_gen_first_generation",
        type=int,
        help="The number of ligands chosen for elitism for the first generation \
        These will advance from the previous generation directly into the next \
        generation.  This is purely advancing based on Docking/Rescore fitness. \
        This does not select for diversity. If not defined it will default to \
        number_elitism_advance_from_previous_gen",
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


def _add_ga_subsequent_gen_params(parser: argparse._ArgumentGroup):
    """
    Add genetic algorithm parameters for subsequent generations.

    This function adds arguments for the number of molecules to seed the next
    generation, number of crossovers, and number of mutants for generations
    after the first.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    parser.add_argument(
        "--top_mols_to_seed_next_generation",
        type=int,
        default=10,
        help="Number of mols that seed next generation, for all generations after the first.\
        Should be less than number_of_crossovers_first_generation \
        + number_of_mutations_first_generation",
    )
    parser.add_argument(
        "--number_of_crossovers",
        type=int,
        default=10,
        help="The number of ligands which will be created via crossover in each \
        generation besides the first",
    )
    parser.add_argument(
        "--number_of_mutants",
        type=int,
        default=10,
        help="The number of ligands which will be created via mutation in each \
        generation besides the first.",
    )


def _add_ga_params(parser: argparse._ArgumentGroup):
    """
    Add general genetic algorithm parameters.

    This function adds arguments for the number of generations, elitism,
    redocking elite compounds, and diversity seed depreciation.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # Populations settings
    parser.add_argument(
        "--num_generations",
        type=int,
        default=10,
        help="The number of generations to be created.",
    )
    parser.add_argument(
        "--number_elitism_advance_from_previous_gen",
        type=int,
        default=10,
        help="The number of ligands chosen for elitism. These will advance from \
        the previous generation directly into the next generation. \
        This is purely advancing based on Docking/Rescore \
        fitness. This does not select for diversity.",
    )
    parser.add_argument(
        "--redock_elite_from_previous_gen",
        choices=[True, False, "True", "False", "true", "false"],
        default=False,
        help="If True than ligands chosen via Elitism (ie advanced from last generation) \
        will be passed through Gypsum and docked again. This provides a better exploration of conformer space \
        but also requires more computation time. If False, advancing ligands are simply carried forward by \
        copying the PDBQT files.",
    )

    parser.add_argument(
        "--diversity_seed_depreciation_per_gen",
        type=int,
        default=2,
        help="Each gen diversity_mols_to_seed_first_generation will decrease this amount",
    )


def _add_conversion_params(parser: argparse._ArgumentGroup):
    """
    Add file conversion parameters.

    This function adds an argument for the path to the OpenBabel executable.
    TODO: Describe nuance here.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # Path to Open Babel file conversion for docking inputs
    parser.add_argument(
        "--obabel_path",
        help="The path to the open babel executable. \
        Path may look like PATH/envs/py37/bin/obabel; \
        may be found on Linux by running: which obabel",
    )


def _add_scoring_params(parser: argparse._ArgumentGroup):
    """
    Add scoring-related parameters.

    This function adds arguments for the scoring choice and ligand efficiency
    rescoring option. TODO: Need to implement this in the future.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # scoring
    parser.add_argument(
        "--scoring_choice",
        metavar="scoring_choice",
        choices=["VINA"],
        default="VINA",
        help="The scoring_choice to use to assess the ligands docking fitness. \
        Default is using Vina/QuickVina2 ligand affinity.",
    )
    parser.add_argument(
        "--rescore_lig_efficiency",
        action="store_true",
        default=False,
        help="This will divide the final scoring_choice output by the number of \
        non-Hydrogen atoms in the ligand. This adjusted ligand efficiency score will \
        override the scoring_choice value. This is compatible with all scoring_choice options.",
    )


# TODO: Use these parameter later when you implement gypsum option.
# def _add_gypsum_params(parser: argparse._ArgumentGroup):
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


def _add_misc_params(parser: argparse._ArgumentGroup):
    """
    Add miscellaneous parameters.

    This function adds an argument for generating a plot at the end of the run.

    Args:
        parser (argparse._ArgumentGroup): The argument group to add the
            parameters to.
    """
    # Make a line plot of the simulation at the end of the run.
    parser.add_argument(
        "--generate_plot",
        choices=[True, False, "True", "False", "true", "false"],
        default=True,
        help="Make a line plot of the simulation at the end of the run.",
    )
