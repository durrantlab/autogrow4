"""Provides base classes and manager for selector plugins in AutoGrow.

This module defines base selector classes and plugin manager for handling compound
selection operations in the genetic algorithm. It includes abstract methods that must
be implemented by specific selector plugins.
"""

from abc import abstractmethod
import random
from typing import Any, List, Optional, Tuple, cast
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompound, ScoreType
from autogrow.utils.logging import LogLevel, log_info


class SelectorBase(PluginBase):
    """Run the plugin with provided arguments.

    Args:
        **kwargs: Dictionary containing:
            predock_cmpds (List[PreDockedCompound]): Available compounds to
                select from
            score_type (ScoreType): Type of score to use for selection
            num_to_choose (int): Number of compounds to select
            favor_most_negative (bool): If True, lower scores are preferred

    Returns:
        List[PreDockedCompound]: Selected compounds
    """

    def run(self, **kwargs) -> List[PreDockedCompound]:
        """Run the plugin with provided arguments."""
        predock_cmpds: List[PreDockedCompound] = kwargs["predock_cmpds"]
        score_type: ScoreType = kwargs["score_type"]
        num_to_choose: int = kwargs["num_to_choose"]
        favor_most_negative: bool = kwargs["favor_most_negative"]

        return self.run_selector(
            predock_cmpds=predock_cmpds,
            num_to_choose=num_to_choose,
            score_type=score_type,
            favor_most_negative=favor_most_negative,
        )

    @abstractmethod
    def run_selector(
        self,
        predock_cmpds: List[PreDockedCompound],
        num_to_choose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompound]:
        """Abstract method for implementing selector-specific compound
        selection logic.

        Args:
            predock_cmpds (List[PreDockedCompound]): Available compounds to
                select from
            num_to_choose (int): Number of compounds to select
            score_type (ScoreType): Type of score to use for selection (docking
                or diversity)
            favor_most_negative (bool, optional): If True, lower scores are
                preferred. Defaults to True.

        Returns:
            List[PreDockedCompound]: Selected compounds
        """
        pass

class SelectorPluginManager(PluginManagerBase):
    def execute(self, **kwargs) -> List[PreDockedCompound]:
        """Execute selector plugin to choose compounds based on scores.

        Runs the selected plugin to choose compounds based on both docking and
        diversity scores. Combines the selected compounds into a single list.

        Args:
            **kwargs: Dictionary containing:
                predock_cmpds (List[PreDockedCompound]): Available compounds
                num_seed_dock_fitness (int): Number to choose by docking score
                num_seed_diversity (int): Number to choose by diversity score
                favor_most_negative (bool): If True, lower scores are preferred

        Returns:
            List[PreDockedCompound]: Combined list of compounds selected by both
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

        docking_fitness_smiles_list: List[PreDockedCompound] = []
        diversity_smile_list: List[PreDockedCompound] = []

        # Select the molecules based on the docking score
        num_predock_cmpds = len(kwargs["predock_cmpds"])
        num_to_choose = kwargs["num_seed_dock_fitness"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by docking score from {num_predock_cmpds} available compounds"
            )
            with LogLevel():
                docking_fitness_smiles_list = selector.run(
                    **{
                        "predock_cmpds": kwargs["predock_cmpds"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DOCKING,
                        "favor_most_negative": kwargs["favor_most_negative"],
                    }
                )

        # Select the molecules based on the diversity score
        num_to_choose = kwargs["num_seed_diversity"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by diversity score from {num_predock_cmpds} available compounds"
            )

            with LogLevel():
                diversity_smile_list = selector.run(
                    **{
                        "predock_cmpds": kwargs["predock_cmpds"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DIVERSITY,
                        "favor_most_negative": kwargs["favor_most_negative"],
                    }
                )

        # Combine the two lists
        docking_diversity_list = list(docking_fitness_smiles_list)
        docking_diversity_list.extend(diversity_smile_list)

        #Shuffle the list to prevent bias
        random.shuffle(docking_diversity_list)

        return docking_diversity_list

        # # Finalize the list
        # return selector.finalize_composite_docking_diversity_list(
        #     docking_diversity_list, kwargs["predock_cmpds"]
        # )
