"""
A minimal, from-scratch ArgumentParser that allows adding the same argument
to multiple groups, avoiding the complexities of subclassing argparse.
"""
import sys
import os
import textwrap
from typing import Any, Dict, List, Optional, FrozenSet

class Namespace:
    """A simple class to hold attributes parsed from the command line."""
    def __init__(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self) -> str:
        return f"Namespace({vars(self)})"

class Action:
    """Represents a single command-line argument's properties."""
    def __init__(self, option_strings: List[str], dest: str, **kwargs: Any) -> None:
        self.option_strings = option_strings
        self.dest = dest
        self.action = kwargs.get('action', 'store')
        self.type = kwargs.get('type', None)
        self.help = kwargs.get('help', None)
        self.metavar = kwargs.get('metavar', None)
        self.choices = kwargs.get('choices', None)

        if self.action == 'store_true':
            self.nargs = 0
            self.default = kwargs.get('default', False)
        else:
            self.nargs = 1
            self.default = kwargs.get('default', None)

class CustomArgumentGroup:
    """A group for command-line arguments for display purposes in the help text."""
    def __init__(self, parser: 'CustomArgumentParser', title: str) -> None:
        self._parser = parser
        self.title = title
        self._group_actions: List[Action] = []

    def add_argument(self, *args: str, **kwargs: Any) -> Action:
        """
        Adds an argument to the group by delegating to the main parser.

        Args:
            *args: The name(s) of the argument (e.g., '-f', '--foo').
            **kwargs: Keyword arguments for the argument (e.g., help, type).

        Returns:
            Action: The Action object representing the argument.
        """
        action = self._parser.add_argument(*args, **kwargs)
        if action not in self._group_actions:
            self._group_actions.append(action)
        return action

class CustomArgumentParser:
    """
    A custom ArgumentParser that allows adding the same argument to multiple
    groups, provided the argument definitions are identical.
    """
    def __init__(self, description: Optional[str] = None) -> None:
        self.description = description
        self._actions: List[Action] = []
        self._option_string_actions: Dict[str, Action] = {}
        self._action_groups: List[CustomArgumentGroup] = []
        self._registered_kwargs: Dict[FrozenSet[str], Dict[str, Any]] = {}

        optional_group = self.add_argument_group('optional arguments')
        help_action = self.add_argument(
            '-h', '--help', action='help', help='show this help message and exit'
        )
        if help_action not in optional_group._group_actions:
            optional_group._group_actions.append(help_action)

    def add_argument(self, *args: str, **kwargs: Any) -> Action:
        """
        Adds an argument to the parser and validates any duplicates.

        Args:
            *args: The name(s) of the argument (e.g., '-f', '--foo').
            **kwargs: Keyword arguments for the argument (e.g., help, type).

        Returns:
            Action: The Action object for the argument.
        
        Raises:
            ValueError: If a duplicate argument is found with conflicting definitions.
        """
        option_strings = [arg for arg in args if arg.startswith('-')]
        
        # Check for duplicates by checking if any of the new option strings are already registered.
        for option in option_strings:
            if option in self._option_string_actions:
                existing_action = self._option_string_actions[option]
                existing_kwargs_key = frozenset(existing_action.option_strings)
                existing_kwargs = self._registered_kwargs.get(existing_kwargs_key)

                if existing_kwargs is not None and kwargs != existing_kwargs:
                    raise ValueError(f"Argument {option} redefined with conflicting options")
                
                return existing_action
        
        # Handle special 'help' action
        if kwargs.get('action') == 'help':
            action = Action(option_strings, 'help', **kwargs)
            self._actions.append(action)
            for option in option_strings:
                self._option_string_actions[option] = action
            return action

        # Determine the destination attribute name
        dest = kwargs.pop('dest', None)
        if dest is None:
            for option in option_strings:
                if option.startswith('--'):
                    dest = option.lstrip('-').replace('-', '_')
                    break
            if dest is None and option_strings:
                dest = option_strings[0].lstrip('-').replace('-', '_')
            if dest is None:
                 raise ValueError("Could not determine destination for argument")


        action_obj = Action(option_strings, dest, **kwargs)
        
        self._actions.append(action_obj)
        for option in option_strings:
            self._option_string_actions[option] = action_obj
        
        self._registered_kwargs[frozenset(option_strings)] = kwargs
        
        return action_obj

    def add_argument_group(self, title: str) -> CustomArgumentGroup:
        """
        Creates and returns a new argument group.

        Args:
            title (str): The title for the argument group.

        Returns:
            CustomArgumentGroup: The new argument group object.
        """
        group = CustomArgumentGroup(self, title)
        self._action_groups.append(group)
        return group
    
    def parse_args(self, args: Optional[List[str]] = None) -> Namespace:
        """
        Parses command-line arguments.

        Args:
            args (Optional[List[str]]): A list of arguments to parse. If None,
                sys.argv[1:] is used.

        Returns:
            Namespace: An object containing the parsed arguments as attributes.
        """
        if args is None:
            args = sys.argv[1:]

        if '-h' in args or '--help' in args:
            self._print_help()
            sys.exit(0)
            
        namespace = Namespace()
        
        for action in self._actions:
            if action.action != 'help':
                setattr(namespace, action.dest, action.default)

        i = 0
        while i < len(args):
            arg = args[i]
            value_str = None
            option_string = arg

            if '=' in arg:
                option_string, value_str = arg.split('=', 1)
            
            action = self._option_string_actions.get(option_string)

            if not action:
                script_name = os.path.basename(sys.argv[0] if sys.argv else 'script.py')
                sys.stderr.write(f"usage: {script_name} [-h] ...\n")
                sys.stderr.write(f"{script_name}: error: unrecognized arguments: {arg}\n")
                sys.exit(2)
            if action.action == 'store':
                if value_str is None:  # value is in the next arg
                    i += 1
                    if i >= len(args):
                        raise ValueError(f"Argument {option_string} requires a value")
                    value_str = args[i]
                
                value = action.type(value_str) if action.type else value_str
                if action.choices and value not in action.choices:
                    choices_str = ", ".join(map(str, action.choices))
                    raise ValueError(f"invalid choice: '{value}' (choose from {choices_str})")
                setattr(namespace, action.dest, value)

            elif action.action == 'store_true':
                if value_str is not None:
                    raise ValueError(f"Argument {option_string} does not take a value")
                setattr(namespace, action.dest, True)

            i += 1
            
        return namespace

    def _print_help(self) -> None:
        """Formats and prints the help message to the console."""
        script_name = os.path.basename(sys.argv[0] if sys.argv else 'script.py')

        # Build usage string parts
        usage_parts = []
        # Get all actions, but filter out help since we handle it separately to ensure it comes first.
        actions_for_usage = [a for a in self._actions if a.dest != 'help']
        
        # Start with help
        usage_parts.append("[-h]")

        for action in actions_for_usage:
            # Prefer the long option string for display
            option_string = next((s for s in action.option_strings if s.startswith('--')), action.option_strings[0])

            if action.action == 'store':
                metavar = action.metavar
                if metavar is None:
                    if action.choices:
                        metavar = "{" + ",".join(map(str, action.choices)) + "}"
                    else:
                        metavar = action.dest.upper()
                usage_parts.append(f"[{option_string} {metavar}]")
            elif action.action == 'store_true':
                usage_parts.append(f"[{option_string}]")

        # Manually wrap usage string
        prefix = f"usage: {script_name} "
        terminal_width = 120 # A bit wider to avoid too many wraps
        
        lines = []
        current_line = ""
        
        for part in usage_parts:
            separator = " " if current_line else ""
            if current_line and (len(current_line) + len(separator) + len(part) > (terminal_width - len(prefix))):
                lines.append(current_line)
                current_line = part
            else:
                current_line += separator + part
        
        if current_line:
            lines.append(current_line)

        # Print the usage message.
        if lines:
            print(prefix + lines[0])
            indent = ' ' * len(prefix)
            for line in lines[1:]:
                print(indent + line)

        if self.description:
            print()
            print(self.description)
        print()

        terminal_width = 80
        # The column where the help text should start.
        help_start_col = 30

        for group in self._action_groups:
            if not group._group_actions:
                continue
            print(f"{group.title}:")

            for action in group._group_actions:
                options_str = ", ".join(action.option_strings)
                if action.action == 'store':
                    metavar = action.metavar
                    if metavar is None:
                        if action.choices:
                            metavar = "{" + ",".join(map(str, action.choices)) + "}"
                        else:
                            metavar = action.dest.upper()
                    options_with_metavar = [f"{opt} {metavar}" for opt in action.option_strings]
                    options_str = ", ".join(options_with_metavar)
                else: # for 'store_true' or 'help'
                    options_str = ", ".join(action.option_strings)

                options_display = f"  {options_str}"
                help_text = ' '.join(textwrap.dedent(action.help or "").strip().split())

                if not help_text:
                    print(options_display)
                    continue
                
                help_indent_str = " " * help_start_col
                help_width = terminal_width - help_start_col

                # If the options part is too long, print it on its own line
                # and start the help text on the next.
                if len(options_display) >= help_start_col - 1:
                    print(options_display)
                    wrapped_lines = textwrap.wrap(help_text, width=help_width)
                    for line in wrapped_lines:
                        print(f"{help_indent_str}{line}")
                else:
                    # Help can start on the same line
                    padded_options = options_display.ljust(help_start_col - 1)
                    wrapped_lines = textwrap.wrap(help_text, width=help_width)
                    
                    if wrapped_lines:
                        print(f"{padded_options} {wrapped_lines[0]}")
                        for line in wrapped_lines[1:]:
                            print(f"{help_indent_str}{line}")
                    else:
                        print(padded_options.rstrip())

            print()