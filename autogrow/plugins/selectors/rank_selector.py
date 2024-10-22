"""
Rank-based selector plugin for AutoGrow.

This module provides the RankSelector class, which implements a selection
strategy based on ranking compounds according to their scores. 

The RankSelector chooses compounds based on their rank in a specified score type
(such as docking or diversity scores), selecting the top-ranked compounds up to
a specified number. This selector is non-redundant, meaning it avoids selecting
duplicate compounds.

Key features:
- Rank-based selection of compounds
- Support for both docking and diversity scores
- Non-redundant selection to ensure diversity
- Configurable to favor either lower or higher scores

Note: This selector is not recommended for small runs where the number of
desired ligands might exceed the number of available ligands to choose from.
"""

import __future__

from autogrow.plugins.selectors import SelectorBase
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.types import PreDockedCompound, ScoreType
from autogrow.utils.logging import log_debug


class RankSelector(SelectorBase):
    """
    A selector plugin that chooses compounds based on their rank in a specified
    score.

    This selector is non-redundant and selects the top-ranked compounds up to
    the specified number. It's not recommended for small runs where the number
    of desired ligands might exceed the number of ligands to choose from.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the plugin
                category name and a list of ArgumentVars for the plugin's specific
                arguments.
        """
        return (
            "Selectors",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Select compounds to advance to the next generation per a ranked selector. The Rank option is a non-redundant selector. Do not use Rank_Selector for small runs as there is potential that the number of desired ligands exceed the number of ligands to chose from.",  # TODO: Add more detail here.
                )
            ],
        )

    def validate(self, params: dict):
        """
        Validate the provided arguments.

        Args:
            params (dict): A dictionary of parameters to validate.

        Note:
            This method is currently a placeholder and needs to be implemented.
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
        Select compounds based on their rank in the specified score.

        This method selects the top-ranked compounds up to the specified number,
        removing any redundancies in the process.

        Args:
            usable_smiles (List[PreDockedCompound]): A list of all compounds
                from the previous generation.
            num_to_choose (int): The number of compounds to select.
            score_type (ScoreType): The type of score to use for ranking (e.g.,
                docking or diversity).
            favor_most_negative (bool): If True, lower scores are considered
                better. Default is True.

        Returns:
            List[PreDockedCompound]: A list of selected compounds.

        Raises:
            Exception: If usable_smiles is not a list, is empty, or if there are
                fewer unique compounds than requested.
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

        # Sort by chosen idx property
        sorted_list = sorted(
            usable_smiles,
            key=lambda x: x.get_previous_score(score_type),
            reverse=not favor_most_negative,
        )

        # sorted_list = sorted(
        #     usable_smiles,
        #     key=lambda x: x.score_by_index_lookup(column_idx_to_select),
        #     reverse=reverse_sort,
        # )

        # remove any redundants
        new_sorted_list: List[PreDockedCompound] = []
        temp_list_info: List[str] = []
        for i in range(len(sorted_list)):
            info = sorted_list[i]
            if "\t".join(info.to_list()) in temp_list_info:
                continue

            temp_list_info.append("\t".join(info.to_list()))
            new_sorted_list.append(info)

        del sorted_list
        del temp_list_info
        if len(new_sorted_list) < num_to_choose:

            raise Exception(
                "Asked for {} but only {} availabe to chose from \
                There are more ligands to chose to seed the list than ligands to select from. \
                Please lower the top_mols_to_seed_next_generation and/or \
                diversity_mols_to_seed_first_generation".format(
                    num_to_choose, len(new_sorted_list)
                )
            )

        new_sorted_list = sorted(
            new_sorted_list,
            key=lambda x: x.get_previous_score(score_type),
            reverse=not favor_most_negative,
        )

        if len(list({x.smiles for x in new_sorted_list})) >= num_to_choose:
            sorted_list = []
            smiles_list = []
            for mol_info in new_sorted_list:
                if mol_info.smiles in smiles_list:
                    continue

                sorted_list.append(mol_info)
                smiles_list.append(mol_info.smiles)
        else:
            sorted_list = new_sorted_list

        top_choice_smiles_in_order = []
        for i in range(num_to_choose):
            smiles = sorted_list[i]
            top_choice_smiles_in_order.append(smiles.smiles)
            scre = smiles.get_previous_score(score_type)
            log_debug(f"{smiles.smiles} ({smiles.name}): score {scre:.2f}")

        return top_choice_smiles_in_order

    def finalize_composite_docking_diversity_list(
        self,
        docking_diversity_list: List[PreDockedCompound],
        usable_smiles: List[PreDockedCompound],
    ) -> List[PreDockedCompound]:
        """
        Finalize the list of selected compounds by retrieving full data for
        each.

        This method takes the list of selected compounds and retrieves their
        full data from the usable_smiles list.

        Args:
            docking_diversity_list (List[PreDockedCompound]): The list of
                selected compounds.
            usable_smiles (List[PreDockedCompound]): The full list of available
                compounds.

        Returns:
            List[PreDockedCompound]: A list of selected compounds with their
                full data.
        """

        # Get all the information about the chosen molecules. chosen_mol_list is
        # 1D list of all chosen ligands chosen_mol_full_data_list is a 1D list
        # with each item of the list having multiple pieces of information such
        # as the ligand name/id, the smiles string, the diversity and docking
        # score...
        return self.get_chosen_mol_full_data_list(docking_diversity_list, usable_smiles)
