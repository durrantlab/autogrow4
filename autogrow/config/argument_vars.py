from dataclasses import dataclass
from typing import Any, List, Optional, Set, Type
import re


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


# Storage for plugin argument groups that need to be added
plugin_arg_groups_to_add = []


def is_snake_case(s: str) -> bool:
    """Checks if a string is in snake_case.

    Args:
        s (str): The string to check.

    Returns:
        bool: True if the string is in snake_case, False otherwise.
    """
    return re.match(r'^[a-z0-9_]+$', s) is not None


def register_argparse_group(title: str, arg_vars: List[ArgumentVars], plugin_name: Optional[str] = None):
    """Register an argument group with its associated arguments.

    This function is used by plugins to add their own argument groups and
    arguments to the main parser. It also validates the naming convention of
    the arguments.

    Args:
        title (str): The title of the argument group.
        arg_vars (List[ArgumentVars]): List of ArgumentVars objects representing
            the arguments to be added to the group.
        plugin_name (Optional[str]): The name of the plugin registering the
            arguments. Used for validation.
    
    Raises:
        ValueError: If an argument name violates the naming convention.
    """
    if plugin_name is not None:
        if not arg_vars:
            pass  # No arguments to validate
        else:
            # 1. First arg should be CamelCase and match plugin name
            enabling_arg = arg_vars[0]
            if enabling_arg.name != plugin_name:
                raise ValueError(
                    f"Plugin '{plugin_name}' has an invalid enabling argument name: "
                    f"'{enabling_arg.name}'. Expected '{plugin_name}'."
                )

            # 2. Other args should be snake_case
            for arg in arg_vars[1:]:
                # remove leading dashes for validation
                arg_name = arg.name.lstrip('-')
                if not is_snake_case(arg_name):
                    raise ValueError(
                        f"Argument name '{arg.name}' in plugin '{plugin_name}' is not snake_case."
                    )
    global plugin_arg_groups_to_add
    plugin_arg_groups_to_add.append((title, arg_vars))