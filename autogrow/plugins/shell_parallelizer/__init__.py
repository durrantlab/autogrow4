from abc import abstractmethod
from dataclasses import dataclass
import glob
import os
import random
from typing import Any, Dict, List, Optional, Tuple, Union, cast
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PostDockedCompound, PreDockedCompound, ScoreType
import autogrow.docking.scoring.execute_scoring_mol as Scoring
import autogrow.docking.ranking.ranking_mol as Ranking
from autogrow.utils.logging import log_debug, log_warning


@dataclass
class ShellCmdResult:
    cmd: str
    return_code: int
    output: str


class ShellParallelizerBase(PluginBase):
    def run(self, **kwargs) -> List[ShellCmdResult]:
        """Run the plugin with provided arguments."""
        nprocs = kwargs.get("nprocs", -1)
        return self.run_cmds_in_parallel(cmds=kwargs["cmds"], nprocs=nprocs)

    @abstractmethod
    def run_cmds_in_parallel(
        self, cmds: List[str], nprocs: int = -1
    ) -> List[ShellCmdResult]:
        """
        run_cmds_in_parallel is needs to be implemented in each class.

        Inputs:
        :param List[str] cmds: A list of shell commands to run in parallel.
        :param int nprocs: The number of processors to use. Default is -1 (use
            all available).

        Returns:
        :returns: List[ShellCmdResult]: A list of ShellCmdResult objects, each
            containing the command, return code, and output.
        """
        pass

    def get_nprocs_to_use(self, nprocs: int) -> int:
        """
        Get the number of processors to use.

        Inputs:
        :param int nprocs: The number of processors to use. Default is -1 (use
            all available).

        Returns:
        :returns: int: The number of processors to use.
        """
        if nprocs == -1:
            if os.cpu_count() is None:
                nprocs = 1
                log_warning("Could not determine the number of CPUs. Defaulting to 1.")
            else:
                nprocs = os.cpu_count()  # type: ignore
        return nprocs


class ShellParallelizerPluginManager(PluginManagerBase):
    def execute(self, **kwargs) -> List[ShellCmdResult]:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: List[ShellCmdResult]: A list of ShellCmdResult objects, each
            containing the command, return code, and output.
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
