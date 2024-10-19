"""
Handles ranking and selection of ligands based on docking scores and diversity.

This module provides functions to process and evaluate ligands from .smi files,
calculate diversity scores, and prepare data structures for further processing
in AutoGrow.
"""
import __future__

import os
from typing import Dict, List, Optional

from autogrow.types import PreDockedCompound, PostDockedCompound, ScoreType
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint  # type: ignore
from rdkit import DataStructs  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.utils.mol_object_handling as MOH


def get_usable_format(infile: str) -> List[PreDockedCompound]:
    """
    Reads and processes a .smi file into a list of PreDockedCompound objects.

    The .smi file must follow this format for each line:
    SMILES<tab>ID[<tab>optional_info...]<tab>docking_score<tab>diversity_score

    Args:
        infile (str): Path to the formatted .smi file to be read.

    Returns:
        List[PreDockedCompound]: List of PreDockedCompound objects with
        information from the .smi file.

    Raises:
        Exception: If the input file does not exist.

    Note:
        - The last two fields are assumed to be docking_score and 
          diversity_score if present.
        - Any fields between ID and docking_score are ignored but allowed.
    """
    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_smiles: List[PreDockedCompound] = []

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
            postDockedCompound = PostDockedCompound.from_list(parts)
            compoundInfo = PreDockedCompound(
                smiles=postDockedCompound.smiles, name=postDockedCompound.id
            )
            if len(parts) > 2:
                compoundInfo.previous_docking_score = postDockedCompound.docking_score
                compoundInfo.previous_diversity_score = (
                    postDockedCompound.diversity_score
                )
            usable_smiles.append(compoundInfo)

    return usable_smiles


def convert_usable_list_to_lig_dict(
    usable_smiles: List[PreDockedCompound],
) -> Optional[Dict[str, PreDockedCompound]]:
    """
    Converts a list of PreDockedCompound objects to a dictionary.

    Args:
        usable_smiles (List[PreDockedCompound]): List of PreDockedCompound
        objects.

    Returns:
        Optional[Dict[str, PreDockedCompound]]: Dictionary with keys as
        'SMILES+ID' and values as PreDockedCompound objects. Returns None if
        input is not a list.

    Note:
        If duplicate keys exist, it keeps the entry with the higher docking
        score (lower is better).
        # TODO: Should this be hardcoded?
    """
    if type(usable_smiles) is not type([]):
        return None

    usable_dict_of_smiles: Dict[str, PreDockedCompound] = {}
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
    postDockedCompoundInfos: List[PostDockedCompound],
) -> List[PostDockedCompound]:
    """
    Calculates diversity scores for a list of PostDockedCompound objects.

    This function computes Morgan Fingerprints for each molecule and calculates
    pairwise Dice Similarity scores. The sum of these scores becomes the
    diversity score for each molecule.

    Args:
        postDockedCompoundInfos (List[PostDockedCompound]): List of
        PostDockedCompound objects to process.

    Returns:
        List[PostDockedCompound]: Input list with updated diversity scores and
        fingerprints.

    Raises:
        AssertionError: If any molecule fails to sanitize or deprotonate.

    Note:
        - Lower diversity scores indicate more diverse molecules. TODO: Verify this is true.
        - Removes any None entries from the input list.
        - Uses Morgan Fingerprints with radius 10 and feature-based encoding.
    """
    postDockedCompoundInfosToKeep: List[PostDockedCompound] = []

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
