"""
This script runs the ranking and selection of ligands.
"""
import __future__

import os
import random
from typing import Dict, List, Optional, Tuple, Union

from autogrow.types import PreDockedCompoundInfo, PostDockedCompoundInfo, ScoreType
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint  # type: ignore
from rdkit import DataStructs  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH
import autogrow.docking.ranking.selecting.rank_selection as Rank_Sel
import autogrow.docking.ranking.selecting.roulette_selection as Roulette_Sel
import autogrow.docking.ranking.selecting.tournament_selection as Tournament_Sel


def create_seed_list(
    usable_smiles: List[PreDockedCompoundInfo],
    num_seed_diversity: int,
    num_seed_dock_fitness: int,
    selector_choice: str,
    tourn_size: float,
) -> List[PreDockedCompoundInfo]:
    """
    this function will take ausable_list_of_smiles which can be derived from
    either the previous generation or a source_compounds_file. Then it will
    select a set of smiles chosen by a weighted function to their
    docking-fitness (docking score) Then it will select a set of smiles chosen
    by a weighted function to their diversity-fitness (docking score) Then it
    will merge these sets of smiles into a single list

    Using the merged list it will make a list of all the smiles in the
    merged-list (chosen_mol_list) with all the other information about each of
    the chosen mols from the usable_smiles

    It will return this list with the complete information of each chosen mol
    (weighted_order_list)

    Inputs:
    :param list usable_smiles: a list with SMILES strings, names, and
        information about the smiles from either the previous generation or the
        source compound list
    :param int num_seed_diversity: the number of seed molecules which come
        from diversity selection
    :param int num_seed_dock_fitness: the number of seed molecules which come
        from eite selection by docking score
    :param int selector_choice: the choice of selector method. Choices are
        Roulette_Selector, Rank_Selector, or Tournament_Selector
    :param float tourn_size: percentage of the total pool of ligands to be
        tested in each tournament.

    Returns:
    :returns: list chosen_mol_full_data_list: a list of all the smiles in a
        weighted ranking ie ["CCCC"  "zinc123"   1    -0.1]
    """

    if selector_choice == "Rank_Selector":
        print("Rank_Selector")
        # This assumes the most negative number is the best option which is
        # true for both This is true for both the diversity score and the
        # docking score. This may need to be adjusted for different scoring
        # functions.

        # Get seed molecules based on docking scores
        docking_fitness_smiles_list = Rank_Sel.run_rank_selector(
            usable_smiles,
            num_seed_dock_fitness,
            ScoreType.DOCKING,
            True,  # TODO: True is hardcoded? Why?
        )

        # Get seed molecules based on diversity scores
        diversity_smile_list = Rank_Sel.run_rank_selector(
            usable_smiles,
            num_seed_diversity,
            ScoreType.DIVERSITY,
            True,  # TODO: True is hardcoded? Why?
        )

    elif selector_choice == "Roulette_Selector":
        print("Roulette_Selector")
        # Get seed molecules based on docking scores
        docking_fitness_smiles_list = Roulette_Sel.spin_roulette_selector(
            usable_smiles, num_seed_dock_fitness, ScoreType.DOCKING
        )

        # Get seed molecules based on diversity scores
        diversity_smile_list = Roulette_Sel.spin_roulette_selector(
            usable_smiles, num_seed_diversity, ScoreType.DIVERSITY
        )

    elif selector_choice == "Tournament_Selector":
        print("Tournament_Selector")
        # This assumes the most negative number is the best option which is
        # true for both This is true for both the diversity score and the
        # docking score. This may need to be adjusted for different scoring
        # functions.

        # Get seed molecules based on docking scores
        docking_fitness_smiles_list = Tournament_Sel.run_Tournament_Selector(
            usable_smiles,
            num_seed_dock_fitness,
            tourn_size,
            ScoreType.DOCKING,
            True,  # TODO: True is hardcoded? Why?
        )

        # Get seed molecules based on diversity scores
        diversity_smile_list = Tournament_Sel.run_Tournament_Selector(
            usable_smiles,
            num_seed_diversity,
            tourn_size,
            ScoreType.DIVERSITY,
            True,  # TODO: True is hardcoded? Why?
        )

    else:
        print(selector_choice)
        raise Exception(
            "selector_choice value is not Roulette_Selector, Rank_Selector, nor Tournament_Selector"
        )

    chosen_mol_list = list(docking_fitness_smiles_list)
    chosen_mol_list.extend(diversity_smile_list)

    if selector_choice in {"Rank_Selector", "Roulette_Selector"}:
        # Get all the information about the chosen molecules. chosen_mol_list
        # is 1D list of all chosen ligands chosen_mol_full_data_list is a 1D
        # list with each item of the list having multiple pieces of
        # information such as the ligand name/id, the smiles string, the
        # diversity and docking score...
        chosen_mol_full_data_list = get_chosen_mol_full_data_list(
            chosen_mol_list, usable_smiles
        )

    elif selector_choice == "Tournament_Selector":
        # Tournament_Selector returns an already full list of ligands so you
        # can skip the get_chosen_mol_full_data_list step
        chosen_mol_full_data_list = chosen_mol_list

    else:
        print(selector_choice)
        raise Exception(
            "selector_choice value is not Roulette_Selector, Rank_Selector, nor Tournament_Selector"
        )

    return chosen_mol_full_data_list


def get_chosen_mol_full_data_list(
    chosen_mol_list: List[PreDockedCompoundInfo],
    usable_smiles: List[PreDockedCompoundInfo],
) -> List[PreDockedCompoundInfo]:
    """
    This function will take a list of chosen molecules and a list of all the
    SMILES which could have been chosen and all of the information about those
    SMILES (ie. ligand name, SMILES string, docking score, diversity score...)

    It will iterated through the list of chosen mols (chosen_mol_list), get
    all the information from the usable_smiles Then it appends the
    corresponding item in usable_smiles to a new list
    weighted_order_list

    --- an issue to be aware of is that there may be redundancies in both
        chosen_mol_list and usable_smiles this causes a many-to-many
        problem so if manipulating this section you need to solve for
        one-to-many
    ---for this reason if this gets altered it will raise an
        AssertionError if the one-to-many is violated.

    It then shuffles the order of the list which to prevent biasing by the
    order of the ligands.

    It will return that list of the chosen molecules in a randomly shuffled
    order.

    Inputs:
    :param list chosen_mol_list: a list of chosen molecules
    :param list usable_smiles: List of all the possibly chosen ligs
        and all the of the info about it (ie. ligand name, SMILES string, docking
        score, diversity score...) ["CCCC"  "zinc123"   1    -0.1  -0.1]

    Returns:
    :returns: list weighted_order_list: a list of all the SMILES with all of
        the associated information in a random order
    """

    sorted_list = sorted(
        usable_smiles, key=lambda x: x.get_previous_score(ScoreType.DOCKING)
    )
    weighted_order_list: List[PreDockedCompoundInfo] = []
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


def get_usable_format(infile: str) -> List[PreDockedCompoundInfo]:
    """
    This code takes a string for an file which is formatted as an .smi file. It
    opens the file and reads in the components into a usable list.

    The .smi must follow the following format for each line:
        MANDATORY INFO
            part 1 is the SMILES string
            part 2 is the SMILES name/ID

        Optional info
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating
                - For default setting this is the Docking score
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)

    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """

    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_smiles: List[PreDockedCompoundInfo] = []

    if os.path.exists(infile) is False:
        print(f"\nFile of Source compounds does not exist: {infile}\n")
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts separated by 4-spaces
            if len(parts) == 1:
                # split line into parts separated by 4-spaces
                parts = line.split("    ")

            # choice_list = [parts[i] for i in range(len(parts))]
            compoundInfo = PreDockedCompoundInfo(smiles=parts[0], name=parts[1])
            if len(parts) > 2:
                compoundInfo.previous_diversity_score = float(parts[-1])
                compoundInfo.previous_docking_score = float(parts[-2])
            usable_smiles.append(compoundInfo)

    return usable_smiles


def convert_usable_list_to_lig_dict(
    usable_smiles: List[PreDockedCompoundInfo],
) -> Optional[Dict[str, PreDockedCompoundInfo]]:
    """
    This will convert a list created by get_usable_format() to a dictionary
    using the ligand smile+lig_id as the key. This makes for faster searching
    in larger Autogrow runs.

    Inputs:
    :param list usable_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow

    Returns:
    :returns: list usable_dict_of_smiles: dict of all the ligand info with a
        key containing both the SMILES string and the unique lig ID
    """

    if type(usable_smiles) is not type([]):
        return None

    usable_dict_of_smiles: Dict[str, PreDockedCompoundInfo] = {}
    for item in usable_smiles:
        key = item.smiles + item.name
        if key in usable_dict_of_smiles and usable_dict_of_smiles[
            key
        ].get_previous_score(ScoreType.DIVERSITY) < item.get_previous_score(
            ScoreType.DOCKING
        ):
            # TODO: Why DIVERSITY vs. DOCKING? Not understanding the comparison here.
            continue
        usable_dict_of_smiles[key] = item
    return usable_dict_of_smiles


##### Called in the docking class ######
def score_and_calc_diversity_scores(
    postDockedCompoundInfos: List[PostDockedCompoundInfo],
) -> List[PostDockedCompoundInfo]:
    """
    This function will take list of molecules which makes up a population. It
    will then create a diversity score for each molecules:
    It creates the diversity score by determining the Morgan Fingerprint for
        each molecule in the population.
    It then compares the fingerprints for every molecule against every
        molecule in a pairwise manner.
        Based on the approach provided on
            http://www.rdkit.org/docs/GettingStartedInPython.html section: "Picking
            Diverse Molecules Using Fingerprints"
        It determines a score of similarity using the RDKit function
            DataStructs.DiceSimilarity
            -The higher the similarity the higher the similarity score
                -ie) if you compare two identical SMILES the similarity score
                    is 1.0. I.e., if you compare 4 identical SMILES the
                    similarity score for each is 4.0.
                -ie) if you compare two completely different SMILES, the score
                    is 0.0

        It sums the similarity score for each pairwise comparison.
            -ie) if there are 15 ligands the max score is 15 the minimum is 0.0
                    with 15.0 if all ligands are identical

        It then appends the diversity score to the molecule list which it
        returns.

        It can raise an AssertionError if there are ligs which fail to
            sanitize or deprotanate.
                -this prevents future errors from occuring in later steps and
                    makes this funciton usable for multiple codes
        It will remove any Nones from the input list

    Inputs:
    :param list molecules_list: list of all molecules in the populations with
    the respective info

    Returns:
    :returns: list molecules_list: list of all molecules in the populations
        with the respective info and append diversity score
    """

    postDockedCompoundInfosToKeep: List[PostDockedCompoundInfo] = []

    for postDockedCompoundInfo in postDockedCompoundInfos:
        if postDockedCompoundInfo is not None:
            smile = postDockedCompoundInfo.smiles
            # name = pair[1]
            try:
                mol = Chem.MolFromSmiles(smile, sanitize=False)
            except Exception:
                mol = None

            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                    def score_and_calc_diversity_scores"
                )

            mol = MOH.check_sanitization(mol)
            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                        def score_and_calc_diversity_scores"
                )

            mol = MOH.try_deprotanation(mol)
            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                        def score_and_calc_diversity_scores"
                )

            postDockedCompoundInfo.mol = mol
            if mol is None:
                print(postDockedCompoundInfo)
                print("None in temp list, skip this one")
                continue
            # mol_list.append(temp)

            postDockedCompoundInfosToKeep.append(postDockedCompoundInfo)
        else:
            print("noneitem in molecules_list in score_and_calc_diversity_scores")

    for postDockedCompoundInfo in postDockedCompoundInfosToKeep:
        fp = GetMorganFingerprint(postDockedCompoundInfo.mol, 10, useFeatures=True)
        postDockedCompoundInfo.fp = fp

    for i in range(len(postDockedCompoundInfosToKeep)):
        diversity_score = 0
        for j in range(len(postDockedCompoundInfosToKeep)):
            if i != j:
                # if DiceSimilarity=1.0 its a perfect match, the smaller the
                # number the more diverse it is. The sum of all of these gives
                # the distance from the normal. The smaller the number means
                # the more distant
                diversity_score = diversity_score + DataStructs.DiceSimilarity(
                    postDockedCompoundInfosToKeep[i].fp,
                    postDockedCompoundInfosToKeep[j].fp,
                )
        postDockedCompoundInfosToKeep[i].diversity_score = diversity_score

    return postDockedCompoundInfos
