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
            usable_smiles (List[PreDockedCompound]): Available compounds to
                select from
            score_type (ScoreType): Type of score to use for selection
            num_to_choose (int): Number of compounds to select
            favor_most_negative (bool): If True, lower scores are preferred

    Returns:
        List[PreDockedCompound]: Selected compounds
    """

    def run(self, **kwargs) -> Any:
        """Run the plugin with provided arguments."""
        usable_smiles: List[PreDockedCompound] = kwargs["usable_smiles"]
        score_type: ScoreType = kwargs["score_type"]
        num_to_choose: int = kwargs["num_to_choose"]
        favor_most_negative: bool = kwargs["favor_most_negative"]

        return self.run_selector(
            usable_smiles=usable_smiles,
            num_to_choose=num_to_choose,
            score_type=score_type,
            favor_most_negative=favor_most_negative,
        )

    @abstractmethod
    def run_selector(
        self,
        usable_smiles: List[PreDockedCompound],
        num_to_choose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompound]:
        """Abstract method for implementing selector-specific compound
        selection logic.

        Args:
            usable_smiles (List[PreDockedCompound]): Available compounds to
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

    @abstractmethod
    def finalize_composite_docking_diversity_list(
        self,
        docking_diversity_list: List[PreDockedCompound],
        usable_smiles: List[PreDockedCompound],
    ) -> List[PreDockedCompound]:
        """Abstract method for finalizing the combined docking and diversity
        selections.

        This method is called after both docking-based and diversity-based
        selections are complete to perform any final processing on the combined
        list.

        Args:
            docking_diversity_list (List[PreDockedCompound]): Combined list of
                compounds selected by both docking and diversity criteria
            usable_smiles (List[PreDockedCompound]): Master list of all
                available compounds with complete data

        Returns:
            List[PreDockedCompound]: Finalized list of selected compounds
        """
        pass

    def get_chosen_mol_full_data_list(
        self,
        chosen_mol_list: List[PreDockedCompound],
        usable_smiles: List[PreDockedCompound],
    ) -> List[PreDockedCompound]:
        """Get complete data for selected molecules from master list.

        Retrieves full compound information for a list of selected molecules by
        matching against a master list. Sorts by docking score, matches by
        SMILES string, and randomly shuffles the final list to prevent order
        bias.

        Args:
            chosen_mol_list (List[PreDockedCompound]): Selected compounds to
                retrieve data for
            usable_smiles (List[PreDockedCompound]): Master list of all
                compounds with complete data

        Returns:
            List[PreDockedCompound]: Randomly shuffled list of selected
                compounds with complete data

        Raises:
            AssertionError: If number of matches found doesn't equal number of
                chosen molecules

        Note:
            - When matching compounds between lists, there may be redundancies
              causing a many-to-many mapping problem. This function enforces
              one-to-one mapping. 
            - Random shuffling prevents biases from compound ordering
        """
        sorted_list = sorted(
            usable_smiles, key=lambda x: x.get_previous_score(ScoreType.DOCKING)
        )
        weighted_order_list: List[PreDockedCompound] = []
        for smile in chosen_mol_list:
            for smile_pair in sorted_list:
                if smile == smile_pair.smiles:
                    weighted_order_list.append(smile_pair)
                    break

        if len(weighted_order_list) != len(chosen_mol_list):
            raise AssertionError(
                "weighted_order_list not the same length as the chosen_mol_list"
            )

        random.shuffle(weighted_order_list)

        return weighted_order_list


class SelectorPluginManager(PluginManagerBase):
    def execute(self, **kwargs) -> List[PreDockedCompound]:
        """Execute selector plugin to choose compounds based on scores.

        Runs the selected plugin to choose compounds based on both docking and
        diversity scores. Combines the selected compounds into a single list.

        Args:
            **kwargs: Dictionary containing:
                usable_smiles (List[PreDockedCompound]): Available compounds
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
        num_usable_smiles = len(kwargs["usable_smiles"])
        num_to_choose = kwargs["num_seed_dock_fitness"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by docking score from {num_usable_smiles} available compounds"
            )
            with LogLevel():
                docking_fitness_smiles_list = selector.run(
                    **{
                        "usable_smiles": kwargs["usable_smiles"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DOCKING,
                        "favor_most_negative": kwargs["favor_most_negative"],
                    }
                )

        # Select the molecules based on the diversity score
        num_to_choose = kwargs["num_seed_diversity"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by diversity score from {num_usable_smiles} available compounds"
            )

            with LogLevel():
                diversity_smile_list = selector.run(
                    **{
                        "usable_smiles": kwargs["usable_smiles"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DIVERSITY,
                        "favor_most_negative": kwargs["favor_most_negative"],
                    }
                )

        # Combine the two lists
        docking_diversity_list = list(docking_fitness_smiles_list)
        docking_diversity_list.extend(diversity_smile_list)

        # Finalize the list
        return selector.finalize_composite_docking_diversity_list(
            docking_diversity_list, kwargs["usable_smiles"]
        )
