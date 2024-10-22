"""Validation utilities for AutoGrow parameters and dependencies.

This module provides a single entry point for validating all AutoGrow
requirements including parameter validation, dependency checks, and required
input verification.
"""
from autogrow.user_params import check_for_required_inputs
from autogrow.validation.validate_dependencies import validate_dependencies
from autogrow.validation.validate_params import validate_params


def validate_all(params: dict) -> None:
    """Validates all AutoGrow requirements before execution.

    Performs complete validation of AutoGrow setup by checking:
    1. Parameter values and types via validate_params()
    2. Required system dependencies via validate_dependencies()
    3. Required input files/parameters via check_for_required_inputs()

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
    check_for_required_inputs(params)

    # TODO: STUFF MISSING HERE?
