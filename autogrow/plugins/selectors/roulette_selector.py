import __future__

from autogrow.plugins.selectors import SelectorBase
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.types import PreDockedCompoundInfo, ScoreType
import numpy.random as rn


class RouletteSelector(SelectorBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
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
        """Validate the provided arguments."""
        pass

    def run_selector(
        self,
        usable_smiles: List[PreDockedCompoundInfo],
        num_to_choose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompoundInfo]:
        """
        Make a list of ligands chosen by a random weighted roulette selection,
        without replacement, weighted by its docking score

        Inputs:
        :param list usable_smiles: a list with all the information of all
            the mols in the previous generation
        :param int number_to_chose: the number of molecules to chose based on
            docking score
        :param str score_type: an string describing either "docking" or
            "diversity" this tells the function how to adjust the weighted scores

        Returns:
        :returns: list top_choice_smile_order: list of ligands chosen by a random
            weighted selection, without replacement, -weighted by its docking score
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

        return rn.choice(
            smiles_list, size=num_to_choose, replace=False, p=probability
        ).tolist()

    def _adjust_scores(
        self,
        usable_smiles: List[PreDockedCompoundInfo],
        score_type: ScoreType,
        favor_most_negative: bool,
    ) -> List[float]:
        """
        This function adjusts the scores appropriately. This is where we weight
        the scores so smaller differences are more pronounced and where we adjust
        for the fact that docking score is better with a more negative number
        while diversity score is the smallest positive number is the most unique.

        Inputs:
        :param list usable_smiles: a list with all the information of all
            the mols in the previous generation
        :param str score_type: an string describing either "docking"
            or "diversity" this tells the function how to adjust the weighted
            scores

        Returns:
        :returns: list adjusted: list of ligand scores which have been weighted
            and adjusted
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
        docking_diversity_list: List[PreDockedCompoundInfo],
        usable_smiles: List[PreDockedCompoundInfo],
    ) -> List[PreDockedCompoundInfo]:
        # Get all the information about the chosen molecules. chosen_mol_list is
        # 1D list of all chosen ligands chosen_mol_full_data_list is a 1D list
        # with each item of the list having multiple pieces of information such
        # as the ligand name/id, the smiles string, the diversity and docking
        # score...
        return self.get_chosen_mol_full_data_list(docking_diversity_list, usable_smiles)
