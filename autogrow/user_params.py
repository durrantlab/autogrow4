"""user_vars
This should contain the functions for defining input variables.
Both the default variables and the user input variables.
This should also validate them.
"""
import __future__
import os

############################################
###### Variables Handlining Settings #######
############################################
def check_for_required_inputs(input_params):
    """
    Confirm all the required inputs were provided.

    Required Variables go here.

    Inputs:
    :param dict input_params: The parameters. A dictionary of {parameter name: value}.
    """
    # Check numbers which may be defined by first generation
    if "top_mols_to_seed_next_generation_first_generation" not in list(
        input_params.keys()
    ):
        if "top_mols_to_seed_next_generation" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["top_mols_to_seed_next_generation"] = 10
            input_params["top_mols_to_seed_next_generation_first_generation"] = 10
        else:
            input_params[
                "top_mols_to_seed_next_generation_first_generation"
            ] = input_params["top_mols_to_seed_next_generation"]

    if "number_of_crossovers_first_generation" not in list(input_params.keys()):
        if "number_of_crossovers" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_of_crossovers"] = 10
            input_params["number_of_crossovers_first_generation"] = 10
        else:
            input_params["number_of_crossovers_first_generation"] = input_params[
                "number_of_crossovers"
            ]

    if "number_of_mutants_first_generation" not in list(input_params.keys()):
        if "number_of_mutants" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_of_mutants"] = 10
            input_params["number_of_mutants_first_generation"] = 10
        else:
            input_params["number_of_mutants_first_generation"] = input_params[
                "number_of_mutants"
            ]

    if "number_elitism_advance_from_previous_gen_first_generation" not in list(
        input_params.keys()
    ):
        if "number_elitism_advance_from_previous_gen" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_elitism_advance_from_previous_gen"] = 10
            input_params[
                "number_elitism_advance_from_previous_gen_first_generation"
            ] = 10
        else:
            input_params[
                "number_elitism_advance_from_previous_gen_first_generation"
            ] = input_params["number_elitism_advance_from_previous_gen"]

    #######################################
    # Check that all required files exist #
    #######################################

    # convert paths to abspath, in case necessary
    input_params["receptor_path"] = os.path.abspath(input_params["receptor_path"])
    input_params["root_output_folder"] = os.path.abspath(
        input_params["root_output_folder"]
    )
    input_params["source_compound_file"] = os.path.abspath(
        input_params["source_compound_file"]
    )

    # Check receptor_path exists
    if os.path.isfile(input_params["receptor_path"]) is False:
        raise NotImplementedError(
            "Receptor file can not be found. File must be a .PDB file."
        )
    if ".pdb" not in input_params["receptor_path"]:
        raise NotImplementedError("receptor_path must be a .PDB file.")

    # Check root_output_folder exists
    if os.path.exists(input_params["root_output_folder"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(input_params["root_output_folder"])
        except Exception as e:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            ) from e

        if os.path.exists(input_params["root_output_folder"]) is False:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

    if os.path.isdir(input_params["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )

    # Check source_compound_file exists
    if os.path.isfile(input_params["source_compound_file"]) is False:
        raise NotImplementedError(
            "source_compound_file can not be found. \
            File must be a tab delineated .smi file."
        )
    if ".smi" not in input_params["source_compound_file"]:
        raise NotImplementedError(
            "source_compound_file must be a \
            tab delineated .smi file."
        )
