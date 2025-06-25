"""
Defines base classes for re-scoring plugins and manages their execution.

It includes abstract base classes for re-scoring plugins and a plugin manager for
handling re-scoring operations.
"""

from abc import abstractmethod
from typing import List, cast
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound


class RescoringBase(PluginBase):
    """
    Abstract base class for re-scoring plugins.

    This class defines the interface for re-scoring plugins and provides a common
    run method that calls the abstract run_rescoring method.
    """

    def run(self, **kwargs) -> List[Compound]:
        """
        Run the re-scoring plugin with provided arguments.

        Args:
            **kwargs: Keyword arguments to be passed to run_docking method.

        Returns:
            List[Compound]: A list of Compound objects
            containing the re-scoring result.
        """
        return self.run_rescoring(docked_cmpds=kwargs["docked_cmpds"])

    @abstractmethod
    def run_rescoring(self, docked_cmpds: List[Compound]) -> List[Compound]:
        """
        Abstract method to be implemented by each re-scoring plugin.

        Args:
            docked_cmpds (List[Compound]): A list of Compound objects.

        Returns:
            List[Compound]: A list of Compound objects, each
            containing the re-scoring result.
        """
        pass


class RescoringPluginManager(PluginManagerBase):
    """
    Manages the execution of re-scoring plugins.

    This class is responsible for selecting and executing re-scoring plugins.
    """

    def execute(self, **kwargs) -> List[Compound]:
        """
        Execute the selected re-scoring plugin with provided arguments.

        Args:
            **kwargs: A dictionary of arguments to pass to the plugin.

        Returns:
            List[Compound]: A list of Compound objects, each
            containing the re-scoring result.
        """
        rescoring_methods = self.get_selected_plugins_from_params()

        if len(rescoring_methods) == 0:
            return kwargs["docked_cmpds"]
        elif len(rescoring_methods) > 1:
            raise Exception(
                f"Only one re-scoring method can be selected at a time! You selected {rescoring_methods}"
            )

        # Get the selector plugin to use
        rescoring_method = cast(
            RescoringBase, self.plugins[rescoring_methods[0]]
        )

        rescores = rescoring_method.run(**kwargs)

        for idx, compound in enumerate(kwargs["docked_cmpds"]):
            original_score = compound.docking_score
            new_score = rescores[idx]

            compound.docking_score = new_score
            compound.add_history(
                "LIGAND_EFFICIENCY",
                f"Original docking score: {original_score:.3f}. New score (ligand efficiency): {new_score:.3f}",
            )
        return kwargs["docked_cmpds"]

