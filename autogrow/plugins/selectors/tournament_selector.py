import __future__
import copy
import math
import random

from autogrow.plugins.selectors import SelectorBase
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.types import PreDockedCompoundInfo, ScoreType


class TournamentSelector(SelectorBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
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
        """Validate the provided arguments."""
        if "tourn_size" not in self.params:
            raise Exception(
                "You are using the tournament selector, but you have not specified the tourn_size parameter."
            )

    def run_selector(
        self,
        usable_smiles: List[PreDockedCompoundInfo],
        num_to_chose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompoundInfo]:
        """
        This runs a tournament style selector given a list of ligands and
        specified metric. It will randomly select ligands for tournaments. The
        best scoring ligand for each of these groups will end up in the
        chosen_ligands list.

        This is done WITHOUT REPLACEMENT. This does provide an opportunity for any
        ligand to make it into the chosen_ligands list even if it doesn't have a
        high score, but that chance of random incorporation decreases as
        tourn_size increases.
            -ie tourn_size=1.0 will randomly pick N number of ligands equal to the
                total number of ligands in the list this means theres a high chance that
                the top ligand will be chosen enter every tournament and will win
                everytime. This could result in a very homogenous choice.

        Inputs:
        :param list list_of_ligands: The list of lists containing info about
            ligands with scores to select from.
        :param int num_to_chose: the number of ligands to be chosen total this
            also is the number of tournaments that will be conducted.
        :param float tourn_size: percentage of the total pool of ligands to be
            tested in each tournament.
        :param int idx_to_sel: the idx within each sublist which will serve as
            the metric for each tournament.
        :param bol favor_most_negative: True if the most negative number is
            the best solution. False if the most positive number is the best
            solution default to True.

        Returns:
        :returns: list chosen_ligands: a list of chosen ligands containing all the
            info for each ligand with potential for redundancy
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
        for _ in range(num_to_chose):
            chosen_ligand = self._run_one_tournament(
                usable_smiles, num_per_tourn, score_type, favor_most_negative
            )
            list_of_ligands_reduced = [
                x for x in list_of_ligands_reduced if x != chosen_ligand
            ]
            chosen_ligands.append(chosen_ligand)

        return chosen_ligands

    def _run_one_tournament(
        self,
        list_of_ligands: List[PreDockedCompoundInfo],
        num_per_tourn: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> PreDockedCompoundInfo:
        """
        This runs a single tournament style selection given a list of ligands and
        specified metric. It will randomly select ligands for the tournament. The
        best scoring ligand from the tournament will be returned.

        This is done WITHOUT REPLACEMENT. This does provide an opportunity for any
        ligand to make it into the chosen_ligands list even if it doesn't have a
        high score, but that chance of random incorporation decreases as
        tourn_size increases.
            -ie tourn_size=1.0 will randomly pick N number of ligands equal to the
                total number of ligands in the list this means theres a high chance
                that the top ligand will be chosen enter every tournament and will
                win everytime. This could result in a very homogenous choice.

            -num_per_tourn is the int(math.ceil(num_ligands * tourn_size)) so that
                it rounds to the nearest int with a minimum values of 1

        Inputs:
        :param list list_of_ligands: The list of lists containing info about
            ligands with scores to select from.
        :param int num_per_tourn: the number of ligands to be tested in each
            tournament.
        :param int idx_to_sel: the idx within each sublist which will serve as
            the metric for each tournament.
        :param bol favor_most_negative: True if the most negative number is
            the best solution. False if the most positive number is the best
            solution default to True.

        Returns:
        :returns: list chosen_option: a list with a single ligand chosen from a
            single tournament
        """

        num_ligands = len(list_of_ligands)

        chosen_option = PreDockedCompoundInfo(smiles="", name="")  # init
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
        docking_diversity_list: List[PreDockedCompoundInfo],
        usable_smiles: List[PreDockedCompoundInfo],
    ) -> List[PreDockedCompoundInfo]:
        # Tournament_Selector returns an already full list of ligands so you can
        # skip the get_chosen_mol_full_data_list step
        return docking_diversity_list