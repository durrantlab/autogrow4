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
from autogrow.config.argument_vars import ArgumentVars
from autogrow.types import Compound, ScoreType
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
        self, predock_cmpds: List[Compound], num_to_choose: int, score_type: ScoreType,
    ) -> List[Compound]:
        """
        Select compounds using tournament selection without replacement.

        This method runs multiple tournaments to select compounds. In each
        tournament, a subset of compounds is randomly chosen, and the
        best-scoring compound wins the tournament.

        Args:
            predock_cmpds (List[Compound]): A list of all compounds
                from the previous generation.
            num_to_choose (int): The number of compounds to select (also the
                number of tournaments to run).
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores for selection.

        Returns:
            List[Compound]: List of selected compounds.

        Raises:
            Exception: If predock_cmpds is not a list or is empty.

        Note:
            The tournament size is determined by the 'tourn_size' parameter,
            which is a fraction of the total number of compounds.
        """
        tourn_size: float = self.params["tourn_size"]

        if type(predock_cmpds) is not type([]):
            raise Exception("list_of_ligands Must be a list, wrong data type")

        num_ligands = len(predock_cmpds)
        if num_ligands == 0:
            raise Exception(
                "list_of_ligands is an empty list. There is nothing to chose from."
            )

        # if score_type != -1 and len(list_of_ligands[0].to_list()) < score_type:
        #     raise Exception(
        #         "The idx to select by does not exist in the provided list_of_ligand."
        #     )

        num_per_tourn = int(math.ceil(num_ligands * tourn_size))

        chosen_predock_cmpds: List[Compound] = []
        list_of_ligands_reduced = copy.deepcopy(predock_cmpds)
        for _ in range(num_to_choose):
            chosen_predock_cmpd = self._run_one_tournament(
                predock_cmpds, num_per_tourn, score_type
            )
            list_of_ligands_reduced = [
                x for x in list_of_ligands_reduced if x != chosen_predock_cmpd
            ]
            chosen_predock_cmpds.append(chosen_predock_cmpd)

            scre = chosen_predock_cmpd.get_score_by_type(score_type)
            log_debug(
                f"{chosen_predock_cmpd.smiles} ({chosen_predock_cmpd.id}): score {scre:.2f}"
            )

        return chosen_predock_cmpds

    def _run_one_tournament(
        self,
        list_of_ligands: List[Compound],
        num_per_tourn: int,
        score_type: ScoreType,
    ) -> Compound:
        """
        Run a single tournament to select one compound.

        This method randomly selects a subset of compounds for the tournament
        and returns the best-scoring compound from this subset.

        Args:
            list_of_ligands (List[Compound]): The list of all compounds
                to select from.
            num_per_tourn (int): The number of compounds to include in the
                tournament.
            score_type (ScoreType): Specifies whether to use "docking" or
                "diversity" scores for selection.

        Returns:
            Compound: The winning compound from the tournament.
        """
        num_ligands = len(list_of_ligands)

        chosen_option = Compound(smiles="", id="")  # init
        temp = []
        for i in range(num_per_tourn):
            temp.append(i)
            if i == 0:
                chosen_option = list_of_ligands[random.randint(0, num_ligands - 1)]
            else:
                choice = list_of_ligands[random.randint(0, num_ligands - 1)]
                if float(chosen_option.get_score_by_type(score_type)) > float(
                    choice.get_score_by_type(score_type)
                ):
                    chosen_option = choice

        return chosen_option
