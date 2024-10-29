"""Validation utilities for AutoGrow parameters and dependencies.

This module provides a single entry point for validating all AutoGrow
requirements including parameter validation, dependency checks, and required
input verification.
"""
import os
from autogrow.validation.validate_dependencies import validate_dependencies
from autogrow.validation.validate_params import validate_params


def validate_all(params: dict) -> None:
    """Validates all AutoGrow requirements before execution.

    Performs complete validation of AutoGrow setup by checking:
    1. Parameter values and types via validate_params()
    2. Required system dependencies via validate_dependencies()
    3. Required input files/parameters via _check_for_required_inputs()

    Args:
        params (dict): Dictionary containing all AutoGrow parameters

    Raises:
        Various exceptions from called validation functions if requirements are
        not met

    Note:
        This function runs all validations to catch all potential issues before
        execution begins.
    """
    validate_params(params)
    validate_dependencies()
    _check_for_required_inputs(params)


def _check_for_required_inputs(input_params):
    """Validates and processes required input parameters for AutoGrow.

    Verifies the presence of all required parameters, sets defaults for missing
    values, and validates file paths. Handles parameters for generation seeding,
    crossovers, mutations, and elitism, with special handling for
    first-generation parameters.

    Args:
        input_params (dict): Dictionary containing parameter names and their
        values.
            Required parameters include:
                - receptor_path: Path to PDB receptor file
                - output_directory: Directory for output files
                - source_compound_file: Path to tab-delimited SMI file
            Optional parameters with defaults (value=10):
                - top_mols_to_seed_next_generation
                - number_of_crossovers
                - number_of_mutants
                - number_elitism_advance_from_previous_gen
            First generation parameters are automatically set based on their
            non-first-generation counterparts.

    Raises:
        NotImplementedError: If any of the following conditions are met:
            - Receptor file doesn't exist or isn't a PDB file
            - Output folder can't be created or accessed
            - Source compound file doesn't exist or isn't a SMI file

    Notes:
        - All file paths are converted to absolute paths
        - Output directory is created if it doesn't exist
        - First generation parameters are set to match their general
          counterparts if not explicitly defined
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
    input_params["output_directory"] = (
        os.path.abspath(input_params["output_directory"]) + os.sep
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

    # Check output_directory exists
    if os.path.exists(input_params["output_directory"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(input_params["output_directory"])
        except Exception as e:
            raise NotImplementedError(
                "output_directory could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            ) from e

        if os.path.exists(input_params["output_directory"]) is False:
            raise NotImplementedError(
                "output_directory could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

    if os.path.isdir(input_params["output_directory"]) is False:
        raise NotImplementedError(
            "output_directory is not a directory. \
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
