"""
Implements parallel execution plugin using Python's multiprocessing.

This module provides a PythonMultiprocessing class that uses Python's
multiprocessing module to run shell commands in parallel, improving performance
for multi-core systems.
"""

from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.shell_parallelizer import ShellCmdResult, ShellParallelizerBase
import multiprocessing
import subprocess


class PythonMultiprocessing(ShellParallelizerBase):
    """
    A plugin that uses Python's multiprocessing to execute commands in parallel.

    This plugin extends ShellParallelizerBase to provide parallel execution
    capabilities using Python's built-in multiprocessing module.
    """

    number_procs = -1

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to Python Multiprocessing Plugin.

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
        self.number_procs = int(params["procs_per_node"])
        if self.number_procs == -1:
            self.number_procs = self.get_nprocs_to_use(self.number_procs)

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
        self, cmds: List[str]
    ) -> List[ShellCmdResult]:
        """
        Run a list of shell commands in parallel using Python's multiprocessing.

        This method uses Python's multiprocessing module to execute multiple
        shell commands concurrently, improving performance on multi-core
        systems.

        Args:
            cmds (List[str]): A list of shell commands to run in parallel.

        Returns:
            List[ShellCmdResult]: A list of ShellCmdResult objects, one for
                each command executed.

        Note:
            If the number of CPUs cannot be determined, it defaults to using a
            single processor and logs a warning.
        """

        results = []
        if self.number_procs == 1:
            for cmd in cmds:
                results.append(self.run_cmd(cmd))
        else:
            with multiprocessing.Pool(processes=self.number_procs) as pool:
                results = pool.map(self.run_cmd, cmds)

        return results
