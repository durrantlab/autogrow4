from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.shell_parallelizer import ShellCmdResult, ShellParallelizerBase
import subprocess
import os

from autogrow.utils.logging import log_warning


class ParallelExecPlugin(ShellParallelizerBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
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
        """Validate the provided arguments."""
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
        Uses parallel exec to run a list of shell commands in parallel.

        Inputs:
        :param List[str] cmds: A list of shell commands to run in parallel.
        :param int nprocs: The number of processors to use. Default is -1 (use
            all available).

        Returns:
        :returns: List[ShellCmdResult]: A list of ShellCmdResult objects for each command.
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
