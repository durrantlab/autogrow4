from typing import Any, Dict, Union
import os


def define_defaults() -> Dict[str, Any]:
    """
    Sets the command-line parameters to their default values.

    Returns:
    :returns: dict default_params: a dictionary of all default variables
    """

    default_params: Dict[str, Any] = {}

    #### OPTIONAL FILE-LOCATION VARIABLES ####
    # (RECOMMEND SETTING TO "" SO AUTOGROW CAN AUTOLOCATE THESE FILES)#
    default_params["obabel_path"] = "obabel"

    # Crossover function
    default_params["max_time_mcs_prescreen"] = 1
    default_params["max_time_mcs_thorough"] = 1
    default_params["min_atom_match_mcs"] = 4
    default_params["protanate_step"] = False

    # Mutation Settings
    default_params["rxn_library_path"] = "click_chem_rxns"

    # processors
    default_params["number_of_processors"] = 1
    default_params["multithread_mode"] = "multithreading"

    # Genetic Algorithm Components
    # default_params["selector_choice"] = "Roulette_Selector"
    # default_params["tourn_size"] = 0.1

    # Seeding next gen and diversity
    default_params["top_mols_to_seed_next_generation_first_generation"] = 10
    default_params["top_mols_to_seed_next_generation"] = 10
    default_params["diversity_mols_to_seed_first_generation"] = 10
    default_params["diversity_seed_depreciation_per_gen"] = 2

    # Populations settings
    default_params["num_generations"] = 10
    default_params["number_of_crossovers_first_generation"] = 10
    default_params["number_of_mutants_first_generation"] = 10
    default_params["number_of_crossovers"] = 10
    default_params["number_of_mutants"] = 10
    default_params["number_elitism_advance_from_previous_gen"] = 10
    default_params["number_elitism_advance_from_previous_gen_first_generation"] = 10
    default_params["redock_elite_from_previous_gen"] = False

    # Filters
    default_params["LipinskiStrictFilter"] = False
    default_params["LipinskiLenientFilter"] = False
    default_params["GhoseFilter"] = False
    default_params["GhoseModifiedFilter"] = False
    default_params["MozziconacciFilter"] = False
    default_params["VandeWaterbeemdFilter"] = False
    default_params["PAINSFilter"] = False
    default_params["NIHFilter"] = False
    default_params["BRENKFilter"] = False
    default_params["No_Filters"] = False
    default_params["alternative_filter"] = None

    # docking
    # default_params["dock_choice"] = "QuickVina2Docking"
    # default_params["vina_like_executable"] = None
    # default_params["docking_exhaustiveness"] = None
    # default_params["docking_num_modes"] = None
    # default_params["docking_timeout_limit"] = 120

    # scoring
    default_params["scoring_choice"] = "VINA"
    default_params["rescore_lig_efficiency"] = False

    # gypsum # max variance is the number of conformers made per ligand
    default_params["max_variants_per_compound"] = 3
    default_params["gypsum_thoroughness"] = 3
    default_params["min_ph"] = 6.4
    default_params["max_ph"] = 8.4
    default_params["pka_precision"] = 1.0
    default_params["gypsum_timeout_limit"] = 10

    # Other vars
    default_params["generate_plot"] = True

    return default_params
