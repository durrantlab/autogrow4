from autogrow.utils.shellcmds import determine_bash_timeout_vs_gtimeout
import os


def define_defaults():
    """
    Sets the command-line parameters to their default values.

    Returns:
    :returns: dict vars: a dictionary of all default variables
    """

    default_vars = {}

    # where we are currently (absolute filepath from route)
    # used for relative pathings
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Some variables which can be manually replaced but defaults
    # point to prepackaged locations.
    ## Neural Network executable for scoring binding
    default_vars["nn1_script"] = os.path.join(
        script_dir, "..", "docking", "scoring", "nn_score_exe", "nnscore1", "NNScore.py"
    )
    # Example: vars['nn1_script'] =
    #    "/PATH/autogrow4/autogrow/docking/scoring/nn_score_exe/nnscore1/NNScore.py"

    default_vars["nn2_script"] = os.path.join(
        script_dir,
        "..",
        "docking",
        "scoring",
        "nn_score_exe",
        "nnscore2",
        "NNScore2.py",
    )
    # Example: vars['nn2_script'] =
    #    "/PATH/autogrow4/autogrow/docking/scoring/nnscore2/NNScore2.py"

    #### OPTIONAL FILE-LOCATION VARIABLES ####
    # (RECOMMEND SETTING TO "" SO AUTOGROW CAN AUTOLOCATE THESE FILES)#

    # PARSER.add_argument('--conversion_choice', choices
    #    = ["MGLTools","obabel"], default="MGLTools",
    default_vars["conversion_choice"] = "MGLToolsConversion"
    default_vars["obabel_path"] = "obabel"
    default_vars["custom_conversion_script"] = ""
    # vars['prepare_ligand4.py'] =
    #   "/PATH/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    default_vars["prepare_ligand4.py"] = ""
    # vars['prepare_receptor4.py'] =
    #   "/PATH/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
    default_vars["prepare_receptor4.py"] = ""
    # vars['mgl_python'] = "/PATH/MGLTools-1.5.4/bin/pythonsh"
    default_vars["mgl_python"] = ""

    # Crossover function
    default_vars["start_a_new_run"] = False
    default_vars["max_time_mcs_prescreen"] = 1
    default_vars["max_time_mcs_thorough"] = 1
    default_vars["min_atom_match_mcs"] = 4
    default_vars["protanate_step"] = False

    # Mutation Settings
    default_vars["rxn_library"] = "click_chem_rxns"
    default_vars["rxn_library_file"] = ""
    default_vars["function_group_library"] = ""
    default_vars["complementary_mol_directory"] = ""

    # processors
    default_vars["number_of_processors"] = 1
    default_vars["multithread_mode"] = "multithreading"

    # Genetic Algorithm Components
    default_vars["selector_choice"] = "Roulette_Selector"
    default_vars["tourn_size"] = 0.1

    # Seeding next gen and diversity
    default_vars["top_mols_to_seed_next_generation_first_generation"] = 10
    default_vars["top_mols_to_seed_next_generation"] = 10
    default_vars["diversity_mols_to_seed_first_generation"] = 10
    default_vars["diversity_seed_depreciation_per_gen"] = 2

    # Populations settings
    default_vars["filter_source_compounds"] = True
    default_vars["use_docked_source_compounds"] = True
    default_vars["num_generations"] = 10
    default_vars["number_of_crossovers_first_generation"] = 10
    default_vars["number_of_mutants_first_generation"] = 10
    default_vars["number_of_crossovers"] = 10
    default_vars["number_of_mutants"] = 10
    default_vars["number_elitism_advance_from_previous_gen"] = 10
    default_vars["number_elitism_advance_from_previous_gen_first_generation"] = 10
    default_vars["redock_elite_from_previous_gen"] = False

    # Filters
    default_vars["LipinskiStrictFilter"] = False
    default_vars["LipinskiLenientFilter"] = False
    default_vars["GhoseFilter"] = False
    default_vars["GhoseModifiedFilter"] = False
    default_vars["MozziconacciFilter"] = False
    default_vars["VandeWaterbeemdFilter"] = False
    default_vars["PAINSFilter"] = False
    default_vars["NIHFilter"] = False
    default_vars["BRENKFilter"] = False
    default_vars["No_Filters"] = False
    default_vars["alternative_filter"] = None

    # docking
    default_vars["dock_choice"] = "QuickVina2Docking"
    default_vars["docking_executable"] = None
    default_vars["docking_exhaustiveness"] = None
    default_vars["docking_num_modes"] = None
    default_vars["docking_timeout_limit"] = 120
    default_vars["custom_docking_script"] = ""

    # scoring
    default_vars["scoring_choice"] = "VINA"
    default_vars["rescore_lig_efficiency"] = False
    default_vars["custom_scoring_script"] = ""

    # gypsum # max variance is the number of conformers made per ligand
    default_vars["max_variants_per_compound"] = 3
    default_vars["gypsum_thoroughness"] = 3
    default_vars["min_ph"] = 6.4
    default_vars["max_ph"] = 8.4
    default_vars["pka_precision"] = 1.0
    default_vars["gypsum_timeout_limit"] = 10

    # Other vars
    default_vars["debug_mode"] = False
    default_vars["reduce_files_sizes"] = False
    default_vars["generate_plot"] = True
    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option in ["timeout", "gtimeout"]:
        default_vars["timeout_vs_gtimeout"] = timeout_option
    else:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
             Autogrow or you may need to execute through Bash."
        )

    return default_vars
