from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.shell_parallelizer import ShellCmdResult, ShellParallelizerBase
import multiprocessing
import subprocess
import os

from autogrow.utils.logging import log_warning


class PythonMultiprocessing(ShellParallelizerBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
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
        """Validate the provided arguments."""
        pass

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
        Uses python multiprocessing to run a list of shell commands in
        parallel.

        Inputs:
        :param List[str] cmds: A list of shell commands to run in parallel.
        :param int nprocs: The number of processors to use. Default is -1 (use
            all available).

        Returns:
        :returns: List[ShellCmdResult]: A list of ShellCmdResult objects for each command.
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
