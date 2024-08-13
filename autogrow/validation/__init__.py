from autogrow.validation.validate_custom_classes import validate_custom_params
from autogrow.validation.validate_dependencies import validate_dependencies
from autogrow.validation.validate_macos import run_macos_notarization
from autogrow.validation.validate_mgltools import validate_mgltools
from autogrow.validation.validate_nnscores import validate_nnscores
from autogrow.validation.validate_params import validate_params


def validate_all(params: dict, printout: str):
    validate_params(params)
    validate_custom_params(params)
    validate_dependencies()
    run_macos_notarization(params)
    validate_mgltools(params, printout)
    printout = validate_nnscores(params, printout)

    # STUFF MISSING HERE?

    return printout
