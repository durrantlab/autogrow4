"""
Implements a parallel execution plugin for AutoGrow using GNU parallel.

This module provides a ParallelExecPlugin class that uses GNU parallel to run
shell commands in parallel, improving performance for multi-core systems.
"""

from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.shell_parallelizer import ShellCmdResult, ShellParallelizerBase
import subprocess
import os

from autogrow.utils.logging import log_warning


class ParallelExecPlugin(ShellParallelizerBase):
    """
    A plugin that uses GNU parallel to execute shell commands in parallel.

    This plugin extends ShellParallelizerBase to provide parallel execution
    capabilities using the GNU parallel utility.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the Parallel Exec Plugin.

        This method defines the command-line arguments that can be used to
        configure the Parallel Exec Plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("Parallel Exec Shell
                  Parallelizer")
                - A list of ArgumentVars objects defining the arguments:
                    1. An argument to enable the Parallel Exec Plugin
                    2. The path to the GNU parallel executable

        Note:
            The default path for the GNU parallel executable is set to
            "/usr/bin/parallel".
        """
        return (
            "Parallel Exec Shell Parallelizer",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Use GNU parallel to run shell commands in parallel. See https://www.gnu.org/software/parallel/",
                ),
                ArgumentVars(
                    name="parallel_exec_path",
                    type=str,
                    default="/usr/bin/parallel",
                    help="Path to the GNU parallel executable.",
                ),
            ],
        )

    def validate(self, params: dict):
        """
        Validate the arguments provided for the Parallel Exec Plugin.

        This method checks if the specified GNU parallel executable exists and
        is executable.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.

        Raises:
            ValueError: If the parallel exec executable is not found or not
                executable at the specified path.
        """
        if not os.path.isfile(params["parallel_exec_path"]):
            raise ValueError(
                f"Parallel exec executable not found at {params['parallel_exec_path']}"
            )
        if not os.access(params["parallel_exec_path"], os.X_OK):
            raise ValueError(
                f"Parallel exec executable at {params['parallel_exec_path']} is not executable"
            )

    def run_cmd(self, cmd: str) -> ShellCmdResult:
        """
        Run a single shell command and return its output.

        :param str cmd: The shell command to run.
        :return: A ShellCmdResult containing command, return_code, and output.
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
        Run a single shell command and return its output.

        This method executes a shell command using subprocess and captures its
        output and return code.

        Args:
            cmd (str): The shell command to run.

        Returns:
            ShellCmdResult: An object containing the command, return code,
                and output (stdout and stderr combined).
        """
        nprocs = self.get_nprocs_to_use(nprocs)

        # Create a temporary file to store commands
        with open("commands.txt", "w") as f:
            for cmd in cmds:
                f.write(cmd + "\n")

        # Construct the parallel exec command without the redirection
        parallel_cmd = [
            self.params["parallel_exec_path"],
            "-j",
            str(nprocs),  # Don't include the "<" here
        ]

        # Open the commands.txt file to pass its contents via stdin
        with open("commands.txt", "r") as infile:
            process = subprocess.Popen(
                parallel_cmd,
                stdin=infile,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            output, error = process.communicate()

        # Clean up the temporary file
        os.remove("commands.txt")

        # Parse the output and create ShellCmdResult objects
        results = []
        for cmd, cmd_output in zip(cmds, output.split("\n\n")):
            return_code = 0 if cmd_output else 1  # Assume success if there's output
            results.append(
                ShellCmdResult(cmd=cmd, return_code=return_code, output=cmd_output)
            )

        return results
