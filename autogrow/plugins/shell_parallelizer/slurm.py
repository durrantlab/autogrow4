"""Implements a Slurm-based parallel execution plugin for AutoGrow using array jobs."""

import os
import pathlib
import sys
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.shell_parallelizer import ShellCmdResult, ShellParallelizerBase
from autogrow.utils.logging import log_info
import random
import string

class Slurm(ShellParallelizerBase):
    """A plugin that uses Slurm array jobs to execute shell commands in parallel.
    
    This plugin submits commands as a single Slurm array job. If the completion
    marker file exists when run, it collects and returns results. Otherwise, it
    submits the job and exits the program.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        return (
            "Slurm Shell Parallelizer",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Use Slurm to run shell commands in parallel",
                ),
                ArgumentVars(
                    name="slurm_template_file",
                    type=str,
                    help="Path to Slurm job template file containing sbatch directives",
                    default=None,
                ),
                ArgumentVars(
                    name="sbatch_path",
                    type=str,
                    help="Path to the sbatch executable",
                    default="sbatch",
                ),
            ],
        )

    def validate(self, params: dict):
        """Validate Slurm template file exists."""
        if "slurm_template_file" not in params:
            raise ValueError(
                "Slurm template file not specified. Use the --slurm_template_file flag."
            )
        if not os.path.exists(params["slurm_template_file"]):
            raise ValueError(
                f"Slurm template file not found: {params['slurm_template_file']}"
            )
        if "sbatch_path" not in params:
            raise ValueError(
                "sbatch executable path not specified. Use the --sbatch_path flag."
            )
        if not os.path.exists(params["sbatch_path"]):
            raise ValueError(f"sbatch executable not found: {params['sbatch_path']}")

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
                to -1 (use all available).

        Returns:
            List[ShellCmdResult]: A list of ShellCmdResult objects, one for
                each command executed.

        Note:
            If the number of CPUs cannot be determined, it defaults to using a
            single processor and logs a warning.
        """

        # Get the current generation directory
        cache_dir = self.params["cur_gen_dir"]

        # Get random prefix (10 random letters)
        prefix = "".join(random.choices(string.ascii_letters, k=10))

        # Define paths for job files
        completion_file = os.path.join(cache_dir, f"{prefix}_slurm_job_complete")
        commands_file = os.path.join(cache_dir, f"{prefix}_slurm_commands.txt")
        array_script = os.path.join(cache_dir, f"{prefix}_slurm_array.sh")

        # If completion file exists, collect and return results
        if os.path.exists(completion_file):
            return self._collect_results(commands_file, cache_dir)

        # Otherwise, submit array job
        self._submit_array_job(cmds, commands_file, array_script, completion_file, self.params["sbatch_path"])

        # Exit program with message
        log_info("\nSlurm array job submitted. Please wait for completion and then restart AutoGrow to continue the process.\n")
        sys.exit(0)

    def _submit_array_job(
        self,
        cmds: List[str],
        commands_file: str,
        array_script: str,
        completion_file: str,
        sbatch_path: str
    ):
        """Create and submit Slurm array job."""
        # Write commands to file
        with open(commands_file, "w") as f:
            for cmd in cmds:
                f.write(f"{cmd}\n")

        # Read template
        template = pathlib.Path(self.params["slurm_template_file"]).read_text()

        # Create array script
        with open(array_script, "w") as f:
            f.write(template + "\n\n")
            # Add array job logic
            f.write(
                """
# Get command for this array task
CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${commands_file}")

# Create output directory for this task
OUT_DIR="${SLURM_SUBMIT_DIR}/task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${OUT_DIR}"

# Run command and capture output
cd "${OUT_DIR}"
eval "${CMD}" > output.txt 2> error.txt

# Last task creates completion file
if [ "${SLURM_ARRAY_TASK_ID}" -eq "${SLURM_ARRAY_TASK_MAX}" ]; then
    touch "${completion_file}"
fi
""".replace(
                    "${commands_file}", commands_file
                ).replace(
                    "${completion_file}", completion_file
                )
            )

        # Submit array job
        num_tasks = len(cmds)
        # os.system(f"{sbatch_path} --array=1-{num_tasks} {array_script}")
        print(f"{sbatch_path} --array=1-{num_tasks} {array_script}")

    def _collect_results(
        self, commands_file: str, cache_dir: str
    ) -> List[ShellCmdResult]:
        """Collect results from completed array job tasks."""
        results = []

        # Read original commands
        with open(commands_file) as f:
            cmds = f.readlines()

        # Collect results from each task
        for i, cmd in enumerate(cmds, 1):
            cmd = cmd.strip()
            task_dir = os.path.join(cache_dir, f"task_{i}")

            # Read output and error files
            output = pathlib.Path(os.path.join(task_dir, "output.txt")).read_text()
            error = pathlib.Path(os.path.join(task_dir, "error.txt")).read_text()
            # Determine return code (0 if output exists and error is empty)
            return_code = 1 if error.strip() else 0

            results.append(
                ShellCmdResult(cmd=cmd, return_code=return_code, output=output + error)
            )

        return results
