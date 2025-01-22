"""Provides base classes and manager for selector plugins in AutoGrow.

This module defines base selector classes and plugin manager for handling compound
selection operations in the genetic algorithm. It includes abstract methods that must
be implemented by specific selector plugins.
"""

from abc import abstractmethod
import random
from typing import Any, List, Optional, Tuple, cast
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound, ScoreType
from autogrow.utils.logging import LogLevel, log_info, log_warning


class SelectorBase(PluginBase):
    """Run the plugin with provided arguments.

    Args:
        **kwargs: Dictionary containing:
            predock_cmpds (List[Compound]): Available compounds to
                select from
            score_type (ScoreType): Type of score to use for selection
            num_to_choose (int): Number of compounds to select

    Returns:
        List[Compound]: Selected compounds
    """

    def run(self, **kwargs) -> List[Compound]:
        """Run the plugin with provided arguments."""
        predock_cmpds: List[Compound] = kwargs["predock_cmpds"]
        score_type: ScoreType = kwargs["score_type"]
        num_to_choose: int = kwargs["num_to_choose"]

        return self.run_selector(
            predock_cmpds=predock_cmpds,
            num_to_choose=num_to_choose,
            score_type=score_type,
        )

    @abstractmethod
    def run_selector(
        self, predock_cmpds: List[Compound], num_to_choose: int, score_type: ScoreType,
    ) -> List[Compound]:
        """
        Abstract method to implement selector-specific compound selection logic.

        Args:
            predock_cmpds (List[Compound]): Available compounds to
                select from
            num_to_choose (int): Number of compounds to select
            score_type (ScoreType): Type of score to use for selection (docking
                or diversity)

        Returns:
            List[Compound]: Selected compounds
        """
        pass


class SelectorPluginManager(PluginManagerBase):
    """Plugin manager for selector plugins in the AutoGrow system."""

    def execute(self, **kwargs) -> List[Compound]:
        """
        Execute selector plugin to choose compounds based on scores.

        Runs the selected plugin to choose compounds based on both docking and
        diversity scores. Combines the selected compounds into a single list.

        Args:
            **kwargs: Dictionary containing:
                predock_cmpds (List[Compound]): Available compounds
                num_seed_dock_fitness (int): Number to choose by docking score
                num_seed_diversity (int): Number to choose by diversity score

        Returns:
            List[Compound]: Combined list of compounds selected by both
                docking and diversity scores

        Raises:
            Exception: If no selector specified or multiple selectors selected
        """
        selectors = self.get_selected_plugins_from_params()

        if selectors is None or len(selectors) == 0:
            raise Exception(
                f"You must specify a selector! Choose from {str(self.plugins.keys())}"
            )
        if len(selectors) > 1:
            raise Exception(
                f"Only one selector can be selected at a time! You selected {selectors}"
            )

        # Get the selector plugin to use
        selector = cast(SelectorBase, self.plugins[selectors[0]])

        docking_fitness_cmpd_list: List[Compound] = []
        diversity_cmpd_list: List[Compound] = []

        # Select the molecules based on the docking score
        num_predock_cmpds = len(kwargs["predock_cmpds"])
        num_to_choose = kwargs["num_seed_dock_fitness"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by docking score from {num_predock_cmpds} available compounds"
            )
            with LogLevel():
                docking_fitness_cmpd_list = selector.run(
                    **{
                        "predock_cmpds": kwargs["predock_cmpds"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DOCKING,
                    }
                )

            for cmpd in docking_fitness_cmpd_list:
                cmpd.add_history(
                    "SELECTOR",
                    f"Selected {cmpd.smiles} due to its docking score: {cmpd.get_score_by_type(ScoreType.DOCKING):.2f}",
                )

        # Select the molecules based on the diversity score
        num_to_choose = kwargs["num_seed_diversity"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by diversity score from {num_predock_cmpds} available compounds"
            )

            with LogLevel():
                diversity_cmpd_list = selector.run(
                    **{
                        "predock_cmpds": kwargs["predock_cmpds"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DIVERSITY,
                    }
                )

            for cmpd in diversity_cmpd_list:
                cmpd.add_history(
                    "SELECTOR",
                    f"Selected {cmpd.smiles} due to its diversity score: {cmpd.get_score_by_type(ScoreType.DIVERSITY):.2f}",
                )

        # Calculate the average docking score of the
        # docking_fitness_smiles_list.
        avg_docking_score = sum(
            x.get_score_by_type(ScoreType.DOCKING) for x in docking_fitness_cmpd_list
        ) / len(docking_fitness_cmpd_list)
        if avg_docking_score > 0:
            log_warning(
                f"Average docking score of selected compounds, +{avg_docking_score:.2f}, was greater than 0. Are you sure the docking plugin is correctly assigning lower (more negative) scores to better compounds?"
            )

        # Combine the two lists
        selected_cmpd_list = list(docking_fitness_cmpd_list)
        selected_cmpd_list.extend(diversity_cmpd_list)

        # Keep only compounds that are unique
        selected_cmpd_list = list({x.smiles: x for x in selected_cmpd_list}.values())

        # Shuffle the list to prevent bias
        random.shuffle(selected_cmpd_list)

        return selected_cmpd_list

        # # Finalize the list
        # return selector.finalize_composite_docking_diversity_list(
        #     docking_diversity_list, kwargs["predock_cmpds"]
        # )
