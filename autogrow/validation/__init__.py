from autogrow.validation.validate_dependencies import validate_dependencies
from autogrow.validation.validate_params import validate_params


def validate_all(params: dict) -> None:
    validate_params(params)
    validate_dependencies()

    # STUFF MISSING HERE?
