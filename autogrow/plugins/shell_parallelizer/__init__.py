"""
Implements base classes and utilities for shell parallelization in AutoGrow.

This module provides the ShellCmdResult dataclass, ShellParallelizerBase
abstract base class, and ShellParallelizerPluginManager for managing shell
parallelization plugins.
"""

from abc import abstractmethod
from dataclasses import dataclass
import glob
import os
import random
from typing import Any, Dict, List, Optional, Tuple, Union, cast
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound, Compound, ScoreType
import autogrow.docking.scoring.execute_scoring_mol as Scoring
import autogrow.docking.ranking.ranking_mol as Ranking
from autogrow.utils.logging import log_debug, log_warning


@dataclass
class ShellCmdResult:
    """
    Represents the result of a shell command execution.

    Attributes:
        cmd (str): The shell command that was executed.
        return_code (int): The return code of the command execution.
        output (str): The combined stdout and stderr output of the command.
    """

    cmd: str
    return_code: int
    output: str


class ShellParallelizerBase(PluginBase):
    """
    An abstract base class for shell parallelizer plugins.

    This class defines the interface for shell parallelizer plugins and provides
    some common utility methods.
    """

    def run(self, **kwargs) -> List[ShellCmdResult]:
        """
        Run the plugin with provided arguments.

        Args:
            **kwargs: Arbitrary keyword arguments. Expected keys:
                - cmds (List[str]): A list of shell commands to run in parallel.
                - nprocs (int, optional): The number of processors to use.

        Returns:
            List[ShellCmdResult]: A list of ShellCmdResult objects for each
                command.
        """
        nprocs = kwargs.get("nprocs", -1)
        return self.run_cmds_in_parallel(cmds=kwargs["cmds"], nprocs=nprocs)

    @abstractmethod
    def run_cmds_in_parallel(
        self, cmds: List[str], nprocs: int = -1
    ) -> List[ShellCmdResult]:
        """
        Run a list of shell commands in parallel.

        This method must be implemented by subclasses.

        Args:
            cmds (List[str]): A list of shell commands to run in parallel.
            nprocs (int, optional): The number of processors to use. Defaults
                to -1 (use all available).

        Returns:
            List[ShellCmdResult]: A list of ShellCmdResult objects for each
                command.
        """
        pass

    def get_nprocs_to_use(self, nprocs: int) -> int:
        """
        Determine the number of processors to use for parallelization.

        Args:
            nprocs (int): The requested number of processors. If -1, uses all
                available.

        Returns:
            int: The number of processors to use.

        Note:
            If the number of CPUs cannot be determined, it defaults to 1 and
            logs a warning.
        """
        if nprocs == -1:
            if os.cpu_count() is None:
                nprocs = 1
                log_warning("Could not determine the number of CPUs. Defaulting to 1.")
            else:
                nprocs = os.cpu_count()  # type: ignore
        return nprocs


class ShellParallelizerPluginManager(PluginManagerBase):
    """
    A plugin manager for shell parallelizer plugins.

    This class manages the selection and execution of shell parallelizer
    plugins.
    """

    def execute(self, **kwargs) -> List[ShellCmdResult]:
        """
        Execute the selected shell parallelizer plugin.

        Args:
            **kwargs: Arbitrary keyword arguments to pass to the selected
                plugin.

        Returns:
            List[ShellCmdResult]: A list of ShellCmdResult objects from the
                executed commands.

        Raises:
            Exception: If no shell parallelizer is specified or if multiple are
                selected.

        Note:
            Only one shell parallelizer can be selected at a time.
        """
        shell_parallelizers = self.get_selected_plugins_from_params()

        if shell_parallelizers is None or len(shell_parallelizers) == 0:
            raise Exception(
                f"You must specify a shell parallelizer program! Choose from {str(self.plugins.keys())}"
            )
        if len(shell_parallelizers) > 1:
            raise Exception(
                f"Only one shell parallelizer can be selected at a time! You selected {shell_parallelizers}"
            )

        # Get the selector plugin to use
        shell_parallelizer = cast(
            ShellParallelizerBase, self.plugins[shell_parallelizers[0]]
        )

        return shell_parallelizer.run(**kwargs)
