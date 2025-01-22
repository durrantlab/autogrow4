from dataclasses import dataclass
from typing import Any, List, Optional, Set, Type

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
registered_args: Set[str] = set()  # Keep track of registered argument names

def register_argparse_group(title: str, arg_vars: List[ArgumentVars]):
    """Register an argument group with its associated arguments.
    
    This function is used by plugins to add their own argument groups and
    arguments to the main parser.
    
    Args:
        title (str): The title of the argument group.
        arg_vars (List[ArgumentVars]): List of ArgumentVars objects representing
            the arguments to be added to the group.
    """
    global plugin_arg_groups_to_add, registered_args
    
    # Filter out any arguments that have already been registered
    new_arg_vars = []
    for arg_var in arg_vars:
        arg_name = arg_var.name.lstrip('-')  # Remove leading dashes
        if arg_name not in registered_args:
            registered_args.add(arg_name)
            new_arg_vars.append(arg_var)
    
    # Only add the group if it has new arguments
    if new_arg_vars:
        plugin_arg_groups_to_add.append((title, new_arg_vars))