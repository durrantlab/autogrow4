"""
Defines base classes for docking plugins and manages their execution.

It includes abstract base classes for docking plugins and a plugin manager for
handling docking operations. The module also provides functionality for ranking
and saving docked compounds.
"""

from abc import abstractmethod
from typing import List, cast
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound, Compound
from autogrow.utils.logging import log_debug, log_warning


class DockingBase(PluginBase):
    """
    Abstract base class for docking plugins.

    This class defines the interface for docking plugins and provides a common
    run method that calls the abstract run_docking method.
    """

    def run(self, **kwargs) -> List[Compound]:
        """
        Run the docking plugin with provided arguments.

        Args:
            **kwargs: Keyword arguments to be passed to run_docking method.

        Returns:
            List[Compound]: A list of Compound objects
            containing docking results.
        """
        return self.run_docking(predocked_cmpds=kwargs["predocked_cmpds"])

    @abstractmethod
    def run_docking(self, predocked_cmpds: List[Compound]) -> List[Compound]:
        """
        Abstract method to be implemented by each docking plugin.

        Args:
            predocked_cmpds (List[Compound]): A list of Compound
                objects to be docked.

        Returns:
            List[Compound]: A list of Compound objects, each
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

    def execute(self, **kwargs) -> List[Compound]:
        """
        Execute the selected docking plugin with provided arguments.

        Args:
            **kwargs: A dictionary of arguments to pass to the plugin.

        Returns:
            List[Compound]: A list of Compound objects, each
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
            try:
                log_debug(
                    f"Docked molecule {post_docked_cmpd.smiles}. Score: {post_docked_cmpd.docking_score: .2f}"
                )
                post_docked_cmpd.add_history(
                    "DOCKING",
                    f"{post_docked_cmpd.smiles} docked with score {post_docked_cmpd.docking_score: .2f}",
                )
            except:
                log_debug(
                    f"Docked molecule {post_docked_cmpd.smiles} has 'docking_score' attribute as None (Null)"
                )

        # # Sanity check: Make sure each output sdf file exists (should be the
        # # docked pose) and that it belongs to the correct generation.
        # import pdb ;pdb.set_trace()
        # TODO: THIS

        # Sanity check: Make sure each output sdf file contains only one model.
        for sdf_filename in [
            c.sdf_path for c in post_docked_cmpds if c.sdf_path is not None
        ]:
            with open(sdf_filename, "r") as f:
                orig_sdf_content = f.read()

            if orig_sdf_content.count("$$$$") > 1:
                log_warning(
                    f"Docked SDF file {sdf_filename} contains more than one molecule. The docking plugin must output SDF files with a single molecule and a single pose. Keeping only the first model."
                )
                with open(sdf_filename, "w") as f:
                    new_sdf_content = orig_sdf_content.split("$$$$")[0] + "$$$$\n"
                    f.write(new_sdf_content)

        # TODO: Validation needed here. Some docking scores could be the same,
        # but they shouldn't all be the same. SDF paths must be updated. And
        # history must be augmented.

        return post_docked_cmpds
