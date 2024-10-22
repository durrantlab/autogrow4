"""
Implements a roulette selector for choosing compounds in AutoGrow.

This module provides a RouletteSelector class that selects compounds to advance
to the next generation using a weighted roulette selection method.
"""

import __future__

from autogrow.plugins.selectors import SelectorBase
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.types import PreDockedCompound, ScoreType
import numpy.random as rn

from autogrow.utils.logging import log_debug


class RouletteSelector(SelectorBase):
    """
    A selector that uses weighted roulette selection to choose compounds.

    This selector chooses compounds to advance to the next generation using a
    weighted roulette selection method. The selection is stochastic and done
    without replacement.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the Roulette Selector.

        This method defines the command-line arguments that can be used to
        configure the Roulette Selector.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("Selectors")
                - A list with one ArgumentVars object defining the argument
                  to enable the Roulette Selector

        Note:
            The Roulette Selector doesn't require additional parameters beyond
            its activation flag.
        """
        return (
            "Selectors",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Select compounds to advance to the next generation per a weighted roulette selector. The roulette selector chooses without replacement and is stoichastic.",  # TODO: Add more detail here.
                )
            ],
        )

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
        self,
        usable_smiles: List[PreDockedCompound],
        num_to_choose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompound]:
        """
        Select compounds using weighted roulette selection without replacement.

        Args:
            usable_smiles (List[PreDockedCompound]): A list of all compounds
                from the previous generation.
            num_to_choose (int): The number of compounds to select based on
                their score.
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores for weighting.
            favor_most_negative (bool): If True, lower scores are considered
                better. Defaults to True.

        Returns:
            List[PreDockedCompound]: List of selected compounds.

        Raises:
            Exception: If usable_smiles is not a list or is empty.

        Note:
            The selection is weighted by the compounds' scores, with the
                weighting adjusted based on the score_type.
        """
        if type(usable_smiles) is not type([]):
            raise Exception("usable_smiles Must be a list, wrong data type")

        num_ligands = len(usable_smiles)
        if num_ligands == 0:
            raise Exception(
                "usable_smiles is an empty list. There is nothing to chose from."
            )

        if num_to_choose <= 0:
            return []

        adjusted = self._adjust_scores(usable_smiles, score_type, favor_most_negative)

        total = sum(adjusted)
        probability = [x / total for x in adjusted]
        smiles_list = [x.smiles for x in usable_smiles]

        chosen_smis = rn.choice(
            smiles_list, size=num_to_choose, replace=False, p=probability
        ).tolist()

        for chosen_smi in chosen_smis:
            log_debug(chosen_smi)

        return chosen_smis

    def _adjust_scores(
        self,
        usable_smiles: List[PreDockedCompound],
        score_type: ScoreType,
        favor_most_negative: bool,
    ) -> List[float]:
        """
        Adjust scores for weighting in the selection process.

        This method adjusts the scores to emphasize smaller differences and
        account for the directionality of the scores (i.e., whether lower or
        higher scores are better).

        Args:
            usable_smiles (List[PreDockedCompound]): A list of all compounds
                from the previous generation.
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores.
            favor_most_negative (bool): If True, lower scores are considered
                better.

        Returns:
            List[float]: A list of adjusted scores.

        Raises:
            Exception: If an invalid score_type is provided.

        Note:
            For diversity scores, the adjustment is (1/x^2) to make smaller
            (more diverse) scores more prominent. For docking scores, the
            adjustment depends on the favor_most_negative parameter.
        """
        if score_type == ScoreType.DIVERSITY:
            weight_scores = [
                x.previous_diversity_score
                for x in usable_smiles
                if x.previous_diversity_score is not None
            ]
            # adjust by squaring the number to make the discrpency larger and
            # invert by dividing 1/x^2 (because the more diverse a mol is the
            # smaller the number)
            adjusted = [(x ** -2) for x in weight_scores]

        elif ScoreType.DOCKING:
            weight_scores = [
                x.previous_docking_score
                for x in usable_smiles
                if x.previous_docking_score is not None
            ]
            # minimum is the most positive value from usable_smiles the
            # more negative the docking score the better the dock

            # TODO: See comment above? Need to account for the fact that the docking
            # score could be reversed depending on docking program. Need to use
            # favor_most_negative

            minimum = max(weight_scores) + 0.1
            minimum = max(minimum, 0)
            adjusted = [(x ** 10) + minimum for x in weight_scores]

        else:
            raise Exception("docking_or_diversity choice not an option")

        return adjusted

    def finalize_composite_docking_diversity_list(
        self,
        docking_diversity_list: List[PreDockedCompound],
        usable_smiles: List[PreDockedCompound],
    ) -> List[PreDockedCompound]:
        """
        Retrieve full data for selected compounds.

        Args:
            docking_diversity_list (List[PreDockedCompound]): List of selected
                compounds.
            usable_smiles (List[PreDockedCompound]): List of all available
                compounds with their full data.

        Returns:
            List[PreDockedCompound]: A list of selected compounds with their
                complete information.
        """
        # Get all the information about the chosen molecules. chosen_mol_list is
        # 1D list of all chosen ligands chosen_mol_full_data_list is a 1D list
        # with each item of the list having multiple pieces of information such
        # as the ligand name/id, the smiles string, the diversity and docking
        # score...
        return self.get_chosen_mol_full_data_list(docking_diversity_list, usable_smiles)
