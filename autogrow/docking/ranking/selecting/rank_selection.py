"""
This script is use to select molecules using a ranked selector
"""
import __future__


import os
import random
from typing import List, Tuple

from autogrow.types import PreDockedCompoundInfo, ScoreType


def run_rank_selector(
    usable_smiles: List[PreDockedCompoundInfo],
    num_to_chose: int,
    score_type: ScoreType,
    favor_most_negative: bool = True,
) -> List[PreDockedCompoundInfo]:
    """
    Given a data set and an idx number to select based on it will select the
    top rank scores for that critera. The number is choses is defined by
    number_to_chose.

    This is an alternative to the weight roulette style selectors.

    Inputs:
    :param list usable_smiles: a list with all the information of all
        the mols in the previous generation
    :param int number_to_chose: the number of molecules to chose based on
        diversity score
    :param int score_type: Whether to consider docking or diversity scores.
    :param bol favor_most_negative:    Set to False if you want to select the most
        positive number is the best choice Set to False if you want to select the
        most negative number

    Returns:
    :returns: list top_choice_smile_order: list of ligands chosen by a elitism
        selection, without replacement,
    """

    if type(usable_smiles) is not type([]):
        raise Exception("usable_smiles Must be a list, wrong data type")

    num_ligands = len(usable_smiles)
    if num_ligands == 0:
        raise Exception(
            "usable_smiles is an empty list. There is nothing to chose from."
        )

    if num_to_chose <= 0:
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
    new_sorted_list: List[PreDockedCompoundInfo] = []
    temp_list_info: List[str] = []
    for i in range(len(sorted_list)):
        info = sorted_list[i]
        if "\t".join(info.to_list()) in temp_list_info:
            continue

        temp_list_info.append("\t".join(info.to_list()))
        new_sorted_list.append(info)

    del sorted_list
    del temp_list_info
    if len(new_sorted_list) < num_to_chose:

        raise Exception(
            "Asked for {} but only {} availabe to chose from \
            There are more ligands to chose to seed the list than ligands to select from. \
            Please lower the top_mols_to_seed_next_generation and/or \
            diversity_mols_to_seed_first_generation".format(
                num_to_chose, len(new_sorted_list)
            )
        )

    new_sorted_list = sorted(
        new_sorted_list,
        key=lambda x: x.get_previous_score(score_type),
        reverse=not favor_most_negative,
    )

    if len(list({x.smiles for x in new_sorted_list})) >= num_to_chose:
        sorted_list = []
        smiles_list = []
        for mol_info in new_sorted_list:
            if mol_info.smiles in smiles_list:
                continue

            sorted_list.append(mol_info)
            smiles_list.append(mol_info.smiles)
    else:
        sorted_list = new_sorted_list

    top_choice_smile_order = []
    for i in range(num_to_chose):
        smiles = sorted_list[i]
        top_choice_smile_order.append(smiles.smiles)

    return top_choice_smile_order
