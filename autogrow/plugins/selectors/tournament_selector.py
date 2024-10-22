"""
Implements a tournament selector for choosing compounds in AutoGrow.

This module provides a TournamentSelector class that selects compounds to
advance to the next generation using a tournament selection method.
"""

import __future__
import copy
import math
import random

from autogrow.plugins.selectors import SelectorBase
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.types import PreDockedCompound, ScoreType
from autogrow.utils.logging import log_debug


class TournamentSelector(SelectorBase):
    """
    A selector that uses tournament selection to choose compounds.

    This selector chooses compounds to advance to the next generation using a
    tournament selection method. The selection is stochastic and done without
    replacement.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the Tournament Selector.

        This method defines the command-line arguments that can be used to
        configure the Tournament Selector.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("Selectors")
                - A list of ArgumentVars objects defining the arguments:
                    1. An argument to enable the Tournament Selector
                    2. The 'tourn_size' parameter to set the tournament size

        Note:
            The 'tourn_size' parameter determines the fraction of the total
            population to include in each tournament.
        """
        return (
            "Selectors",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Select compounds to advance to the next generation per a tournament selector. The tournament selector chooses without replacement and is stoichastic.",  # TODO: Add more detail here.
                ),
                ArgumentVars(
                    name="tourn_size",
                    type=float,
                    default=0.1,
                    help="If using the Tournament_Selector, this determines the size of each tournament. The number of ligands used for each tournament will (tourn_size) * (the number of considered ligands).",
                ),
            ],
        )

    def validate(self, params: dict):
        """
        Validate the arguments provided for the Tournament Selector.

        This method checks if the required 'tourn_size' parameter is present in
        the provided parameters.

        Args:
            params (dict): A dictionary of parameters provided to the selector.

        Raises:
            Exception: If the 'tourn_size' parameter is not specified when using
                the Tournament Selector.
        """
        if "tourn_size" not in self.params:
            raise Exception(
                "You are using the tournament selector, but you have not specified the tourn_size parameter."
            )

    def run_selector(
        self,
        usable_smiles: List[PreDockedCompound],
        num_to_choose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompound]:
        """
        Select compounds using tournament selection without replacement.

        This method runs multiple tournaments to select compounds. In each
        tournament, a subset of compounds is randomly chosen, and the
        best-scoring compound wins the tournament.

        Args:
            usable_smiles (List[PreDockedCompound]): A list of all compounds
                from the previous generation.
            num_to_choose (int): The number of compounds to select (also the
                number of tournaments to run).
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores for selection.
            favor_most_negative (bool): If True, lower scores are considered
                better. Defaults to True.

        Returns:
            List[PreDockedCompound]: List of selected compounds.

        Raises:
            Exception: If usable_smiles is not a list or is empty.

        Note:
            The tournament size is determined by the 'tourn_size' parameter,
            which is a fraction of the total number of compounds.
        """
        tourn_size: float = self.params["tourn_size"]

        if type(usable_smiles) is not type([]):
            raise Exception("list_of_ligands Must be a list, wrong data type")

        num_ligands = len(usable_smiles)
        if num_ligands == 0:
            raise Exception(
                "list_of_ligands is an empty list. There is nothing to chose from."
            )

        # if score_type != -1 and len(list_of_ligands[0].to_list()) < score_type:
        #     raise Exception(
        #         "The idx to select by does not exist in the provided list_of_ligand."
        #     )

        num_per_tourn = int(math.ceil(num_ligands * tourn_size))

        chosen_ligands = []
        list_of_ligands_reduced = copy.deepcopy(usable_smiles)
        for _ in range(num_to_choose):
            chosen_ligand = self._run_one_tournament(
                usable_smiles, num_per_tourn, score_type, favor_most_negative
            )
            list_of_ligands_reduced = [
                x for x in list_of_ligands_reduced if x != chosen_ligand
            ]
            chosen_ligands.append(chosen_ligand)

            scre = chosen_ligand.get_previous_score(score_type)
            log_debug(
                f"{chosen_ligand.smiles} ({chosen_ligand.name}): score {scre:.2f}"
            )

        return chosen_ligands

    def _run_one_tournament(
        self,
        list_of_ligands: List[PreDockedCompound],
        num_per_tourn: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> PreDockedCompound:
        """
        Run a single tournament to select one compound.

        This method randomly selects a subset of compounds for the tournament
        and returns the best-scoring compound from this subset.

        Args:
            list_of_ligands (List[PreDockedCompound]): The list of all compounds
                to select from.
            num_per_tourn (int): The number of compounds to include in the
                tournament.
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores for selection.
            favor_most_negative (bool): If True, lower scores are considered
                better. Defaults to True.

        Returns:
            PreDockedCompound: The winning compound from the tournament.
        """
        num_ligands = len(list_of_ligands)

        chosen_option = PreDockedCompound(smiles="", name="")  # init
        temp = []
        for i in range(num_per_tourn):
            temp.append(i)
            if i == 0:
                chosen_option = list_of_ligands[random.randint(0, num_ligands - 1)]
            else:
                choice = list_of_ligands[random.randint(0, num_ligands - 1)]
                if favor_most_negative:
                    if float(chosen_option.get_previous_score(score_type)) > float(
                        choice.get_previous_score(score_type)
                    ):
                        chosen_option = choice
                elif float(chosen_option.get_previous_score(score_type)) < float(
                    choice.get_previous_score(score_type)
                ):
                    chosen_option = choice
                else:
                    continue

        return chosen_option

    def finalize_composite_docking_diversity_list(
        self,
        docking_diversity_list: List[PreDockedCompound],
        usable_smiles: List[PreDockedCompound],
    ) -> List[PreDockedCompound]:
        """
        Finalize the list of selected compounds.

        For the Tournament Selector, this method simply returns the input list
        as the selection process is already complete.

        Args:
            docking_diversity_list (List[PreDockedCompound]): List of selected
                compounds from the tournament selection process.
            usable_smiles (List[PreDockedCompound]): List of all available
                compounds (not used in this method).

        Returns:
            List[PreDockedCompound]: The input docking_diversity_list,
                unmodified.
        """

        # TODO: What is the point of this?

        # Tournament_Selector returns an already full list of ligands so you can
        # skip the get_chosen_mol_full_data_list step
        return docking_diversity_list
