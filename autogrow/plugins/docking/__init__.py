"""
This module defines base classes for docking plugins and manages their
execution.

It includes abstract base classes for docking plugins and a plugin manager for
handling docking operations. The module also provides functionality for ranking
and saving docked compounds.
"""

from abc import abstractmethod
import os
from typing import List, cast
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PostDockedCompound, PreDockedCompound
from autogrow.utils.logging import log_debug


class DockingBase(PluginBase):
    """
    Abstract base class for docking plugins.

    This class defines the interface for docking plugins and provides a common
    run method that calls the abstract run_docking method.
    """

    def run(self, **kwargs) -> List[PostDockedCompound]:
        """
        Run the docking plugin with provided arguments.

        Args:
            **kwargs: Keyword arguments to be passed to run_docking method.

        Returns:
            List[PostDockedCompound]: A list of PostDockedCompound objects
            containing docking results.
        """
        return self.run_docking(predocked_cmpds=kwargs["predocked_cmpds"])

    @abstractmethod
    def run_docking(
        self, predocked_cmpds: List[PreDockedCompound]
    ) -> List[PostDockedCompound]:
        """
        Abstract method to be implemented by each docking plugin.

        Args:
            predocked_cmpds (List[PreDockedCompound]): A list of PreDockedCompound
                objects to be docked.

        Returns:
            List[PostDockedCompound]: A list of PostDockedCompound objects, each
            containing the score and a docked (posed) SDF file.
        """
        # raise NotImplementedError("run_dock() not implemented")
        pass


class DockingPluginManager(PluginManagerBase):
    """
    Manages the execution of docking plugins.

    This class is responsible for selecting and executing docking plugins,
    as well as ranking and saving the output of docking operations.
    """

    def execute(self, **kwargs) -> List[PostDockedCompound]:
        """
        Execute the selected docking plugin with provided arguments.

        Args:
            **kwargs: A dictionary of arguments to pass to the plugin.

        Returns:
            List[PostDockedCompound]: A list of PostDockedCompound objects, each
            containing the score and a docked (posed) SDF file.

        Raises:
            Exception: If no docking program is specified or if multiple docking
            programs are selected.
        """
        dockings = self.get_selected_plugins_from_params()

        if dockings is None or len(dockings) == 0:
            raise Exception(
                f"You must specify a docking program! Choose from {str(self.plugins.keys())}"
            )
        if len(dockings) > 1:
            raise Exception(
                f"Only one docking program can be selected at a time! You selected {dockings}"
            )

        # Get the selector plugin to use
        docking = cast(DockingBase, self.plugins[dockings[0]])

        post_docked_cmpds = docking.run(**kwargs)

        for post_docked_cmpd in post_docked_cmpds:
            log_debug(
                f"Docked molecule {post_docked_cmpd.smiles}. Score: {post_docked_cmpd.docking_score:.2f}"
            )

        return post_docked_cmpds



