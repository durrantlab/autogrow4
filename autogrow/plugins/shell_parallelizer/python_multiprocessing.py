"""
Implements a parallel execution plugin for AutoGrow using Python's
multiprocessing.

This module provides a PythonMultiprocessing class that uses Python's
multiprocessing module to run shell commands in parallel, improving performance
for multi-core systems.
"""

from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.shell_parallelizer import ShellCmdResult, ShellParallelizerBase
import multiprocessing
import subprocess
import os

from autogrow.utils.logging import log_warning


class PythonMultiprocessing(ShellParallelizerBase):
    """
    A plugin that uses Python's multiprocessing to execute shell commands in
    parallel.

    This plugin extends ShellParallelizerBase to provide parallel execution
    capabilities using Python's built-in multiprocessing module.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the Python Multiprocessing
        Plugin.

        This method defines the command-line arguments that can be used to
        configure the Python Multiprocessing Plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("Python Multiprocessing Shell
                  Parallelizer")
                - A list with one ArgumentVars object defining the argument
                  to enable the Python Multiprocessing Plugin

        Note:
            This plugin doesn't require additional parameters beyond its
            activation flag.
        """
        return (
            "Python Multiprocessing Shell Parallelizer",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Use Python's multiprocessing module to run shell commands in parallel.",
                )
            ],
        )

    def validate(self, params: dict):
        """
        Validate the arguments provided for the Python Multiprocessing Plugin.

        This method is a placeholder for argument validation. Currently, the
        Python Multiprocessing Plugin doesn't require any additional validation
        beyond its activation.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.
                Not used in the current implementation.
        """
        pass

    def run_cmd(self, cmd: str) -> ShellCmdResult:
        """
        Run a single shell command and return its output.

        This method executes a shell command using subprocess and captures its
        output and return code.

        Args:
            cmd (str): The shell command to run.

        Returns:
            ShellCmdResult: An object containing the command, return code,
                and output (stdout and stderr combined).
        """
        process = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output, error = process.communicate()
        return_code = process.returncode
        return ShellCmdResult(cmd=cmd, return_code=return_code, output=output + error)

    def run_cmds_in_parallel(
        self, cmds: List[str], nprocs: int = -1
    ) -> List[ShellCmdResult]:
        """
        Run a list of shell commands in parallel using Python's multiprocessing.

        This method uses Python's multiprocessing module to execute multiple
        shell commands concurrently, improving performance on multi-core
        systems.

        Args:
            cmds (List[str]): A list of shell commands to run in parallel.
            nprocs (int, optional): The number of processors to use. Defaults
                to -1, which uses all available processors.

        Returns:
            List[ShellCmdResult]: A list of ShellCmdResult objects, one for
                each command executed.

        Note:
            If the number of CPUs cannot be determined, it defaults to using a
            single processor and logs a warning.
        """
        if nprocs == -1:
            if os.cpu_count() is None:
                nprocs = 1
                log_warning("Could not determine the number of CPUs. Defaulting to 1.")
            else:
                nprocs = os.cpu_count()  # type: ignore

        with multiprocessing.Pool(processes=nprocs) as pool:
            results = pool.map(self.run_cmd, cmds)

        return results
