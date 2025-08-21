"""
Implements a roulette selector for choosing compounds in AutoGrow.

This module provides a RouletteSelector class that selects compounds to advance
to the next generation using a weighted roulette selection method.
"""

import __future__

from autogrow.plugins.selectors import SelectorBase
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.types import Compound, ScoreType
import numpy.random as rn

from autogrow.utils.logging import log_debug


class RouletteSelector(SelectorBase):
    """
    A selector that uses weighted roulette selection to choose compounds.

    This selector chooses compounds to advance to the next generation using a
    weighted roulette selection method. The selection is stochastic and done
    without replacement.
    """

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments specific to the Roulette Selector.

        This method defines the command-line arguments that can be used to
        configure the Roulette Selector.

        Returns:
            List[ArgumentVars]: A list with one ArgumentVars object defining
            the argument to enable the Roulette Selector.
        Note:
            The Roulette Selector doesn't require additional parameters beyond
            its activation flag.
        """
        return [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable weighted roulette wheel selection. This stochastic method selects compounds based on their scores, where higher-scoring compounds have a proportionally higher chance of being chosen. Selection is performed without replacement.",
            )
        ]

    def validate(self, params: dict):
        """
        Validate the arguments provided for the Roulette Selector.

        This method is a placeholder for argument validation. Currently, the
        Roulette Selector doesn't require any additional validation beyond its
        activation.

        Args:
            params (dict): A dictionary of parameters provided to the selector.
                Not used in the current implementation.
        """
        pass

    def run_selector(
        self, predock_cmpds: List[Compound], num_to_choose: int, score_type: ScoreType,
    ) -> List[Compound]:
        """
        Select compounds using weighted roulette selection without replacement.

        Args:
            predock_cmpds (List[Compound]): A list of all compounds
                from the previous generation.
            num_to_choose (int): The number of compounds to select based on
                their score.
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores for weighting.

        Returns:
            List[Compound]: List of selected compounds.

        Raises:
            Exception: If predock_cmpds is not a list or is empty.

        Note:
            The selection is weighted by the compounds' scores, with the
                weighting adjusted based on the score_type.
        """
        if type(predock_cmpds) is not type([]):
            raise Exception("predock_cmpds Must be a list, wrong data type")

        num_ligands = len(predock_cmpds)
        if num_ligands == 0:
            raise Exception(
                "predock_cmpds is an empty list. There is nothing to chose from."
            )

        if num_to_choose <= 0:
            return []

        adjusted_scores = self._adjust_scores(predock_cmpds, score_type)

        total = sum(adjusted_scores)
        probabilities = [x / total for x in adjusted_scores]
        # smiles_list = [x.smiles for x in predock_cmpds]

        # NOTE: I'm not sure why type: ignore is needed below. It seems like it should be
        # inferred correctly.
        chosen_predock_cmpds: List[Compound] = rn.choice(
            predock_cmpds, size=num_to_choose, replace=False, p=probabilities  # type: ignore
        ).tolist()  # type: ignore

        for chosen_predock_cmpd in chosen_predock_cmpds:
            log_debug(chosen_predock_cmpd.smiles)

        return chosen_predock_cmpds

    def _adjust_scores(
        self, predock_cmpds: List[Compound], score_type: ScoreType
    ) -> List[float]:
        """
        Adjust scores for weighting in the selection process.

        This method adjusts the scores to emphasize smaller differences and
        account for the directionality of the scores (i.e., whether lower or
        higher scores are better).

        Args:
            predock_cmpds (List[Compound]): A list of all compounds
                from the previous generation.
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores.

        Returns:
            List[float]: A list of adjusted scores (positive weights).
        Raises:
            Exception: If an invalid score_type is provided.

        Note:
            For DIVERSITY, lower scores are better. We use 1/score^2 to create
            positive weights where smaller scores get larger weights.
            For FITNESS, lower (more negative) scores are better. We transform them
            into positive weights using: max_score - score + epsilon, so the
            best scores get the highest weights.
        """
        if score_type == ScoreType.DIVERSITY:
            weight_scores = [
                x.diversity_score
                for x in predock_cmpds
                if x.diversity_score is not None
            ]
            # Use 1/x^2 to make smaller (more diverse) scores have larger weights
            # Add a small epsilon to avoid division by zero
            adjusted = [(1 / (x**2 + 1e-6)) for x in weight_scores]
        elif score_type == ScoreType.FITNESS:
            weight_scores = [
                x.fitness_score for x in predock_cmpds if x.fitness_score is not None
            ]
            # To handle negative scores where lower is better, we transform them.
            # A common method is to subtract scores from the max score.
            if not weight_scores:
                return []
            max_score = max(weight_scores)
            # Add a small epsilon to ensure all weights are non-zero
            epsilon = 1e-6
            adjusted = [(max_score - x + epsilon) for x in weight_scores]
        else:
            raise Exception("docking_or_diversity choice not an option")

        return adjusted
