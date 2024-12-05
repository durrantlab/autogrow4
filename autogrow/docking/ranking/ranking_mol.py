"""
Handles ranking and selection of ligands based on docking scores and diversity.

This module provides functions to process and evaluate ligands from .smi files,
calculate diversity scores, and prepare data structures for further processing
in AutoGrow.
"""
import __future__

import os
from typing import Any, Dict, List, Optional

from autogrow.types import Compound, Compound, ScoreType
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint  # type: ignore
from rdkit import DataStructs  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.utils.mol_object_handling as MOH
import autogrow.docking.ranking.ranking_mol as Ranking


def rank_and_save_output_smi(
    current_generation_dir: str,
    current_gen_int: int,
    smiles_file: str,
    postDockedCompoundInfos: List[Compound],
    params: Dict[str, Any],
) -> str:
    """
    Rank and save docked compounds based on docking score.

    This method ranks all the SMILES based on docking score (low to high),
    formats them into a .smi file, and saves the file.

    Args:
        current_generation_dir (str): Path of directory of current generation.
        current_gen_int (int): The integer of the current generation
            (zero-indexed).
        smiles_file (str): File path for the file with the ligands for the
            generation (a .smi file).
        postDockedCompoundInfos (List[Compound]): List of
            Compound objects containing docking results.
        params (Dict[str, Optional[str]]): Dictionary of parameters.

    Returns:
        str: The path of the output ranked .smi file.

    Note:
        This method handles pass-through ligands from the previous generation
        and the current generation number.
    """
    # Get directory string of PDB files for Ligands
    # folder_with_pdbqts = f"{current_generation_dir}PDBs{os.sep}"

    # Run any compatible Scoring Function
    # postDockedCompoundInfos = Scoring.run_scoring_common(
    #     params, smiles_file, folder_with_pdbqts
    # )

    # Before ranking these we need to handle Pass-Through ligands from the last
    # generation
    #
    # current_gen_int==1: dock all ligands from the last generation so all of
    # the pass-through lig are already in the PDB's folder thus they should be
    # accounted for in smiles_list
    #
    # current_gen_int != 1: We need to append the scores form the last gen to
    # smiles_list

    # Only add these when we haven't already redocked the ligand
    if current_gen_int != 0:
        # Go to previous generation folder
        prev_gen_num = str(current_gen_int - 1)
        run_folder = params["output_directory"]
        previous_gen_folder = f"{run_folder}generation_{prev_gen_num}{os.sep}"
        ranked_smi_file_prev_gen = (
            f"{previous_gen_folder}generation_{prev_gen_num}_ranked.smi"
        )

        # Also check sometimes Generation 1 won't have a previous
        # generation to do this with and sometimes it will
        if (
            current_gen_int != 1
            or os.path.exists(ranked_smi_file_prev_gen) is not False
        ):
            # Shouldn't happen but to be safe.
            _process_ligand_scores_from_prev_gen(
                ranked_smi_file_prev_gen,
                current_generation_dir,
                current_gen_int,
                postDockedCompoundInfos,
            )
    # Output format of the .smi file will be: SMILES    Full_lig_name
    # shorthandname   ...AnyCustominfo... Fitness_metric  diversity
    # Normally the docking score is the fitness metric but if we use a
    # Custom metric than dock score gets moved to index -3 and the new
    # fitness metric gets -2

    # sort list by the affinity of each sublist (which is the last index
    # of sublist)
    postDockedCompoundInfos.sort(key=lambda x: x.docking_score, reverse=False)

    # score the diversity of each ligand compared to the rest of the
    # ligands in the group this adds on a float in the last column for the
    # sum of pairwise comparisons the lower the diversity score the more
    # unique a molecule is from the other mols in the same generation
    postDockedCompoundInfos = Ranking.calc_diversity_scores(postDockedCompoundInfos)

    # name for the output file
    output_ranked_smile_file = smiles_file.replace(".smi", "") + "_ranked.smi"

    # save to a new output smiles file. ie. save to ranked_smiles_file

    with open(output_ranked_smile_file, "w") as output:
        for cmpdInf in postDockedCompoundInfos:
            output.write(cmpdInf.tsv_line)
            # sdf_path = "" if cmpdInf.sdf_path is None else cmpdInf.sdf_path
            # sdf_basename = os.path.basename(sdf_path)
            # output_line = f"{cmpdInf.smiles}\t{cmpdInf.id}\t{cmpdInf.additional_info}\t{cmpdInf.docking_score}\t{cmpdInf.diversity_score}\t{sdf_basename}\n"
            # output.write(output_line)

    return output_ranked_smile_file


def _process_ligand_scores_from_prev_gen(
    ranked_smi_file_prev_gen: str,
    current_generation_dir: str,
    current_gen_int: int,
    smiles_list: List[Compound],
):
    """
    Process ligand scores from the previous generation.

    This method retrieves ligand scores from the previous generation and
    updates the current generation's smiles_list with pass-through ligands.

    Args:
        ranked_smi_file_prev_gen (str): Path to the ranked .smi file from the
            previous generation.
        current_generation_dir (str): Path to the current generation directory.
        current_gen_int (int): The current generation number.
        smiles_list (List[Compound]): List of Compound
            objects to be updated.

    Raises:
        Exception: If the previous generation ranked .smi file does not exist.

    Note:
        This method modifies the smiles_list in place.
    """

    # CHECKED: smiles_list is of type List[Compound] here.

    print("Getting ligand scores from the previous generation")

    # Shouldn't happen but to be safe.
    if os.path.exists(ranked_smi_file_prev_gen) is False:
        raise Exception(
            "Previous generation ranked .smi file does not exist. "
            + "Check if output folder has been moved"
        )

    # Get the data for all ligands from previous generation ranked
    # file
    prev_gen_data_list = Ranking.get_predockcmpds_from_smi_file(
        ranked_smi_file_prev_gen
    )
    # CHECKED: prev_gen_data_list of type List[Compound] here.

    # Get the list of pass through ligands
    current_gen_pass_through_smi = (
        current_generation_dir
        + f"SeedFolder{os.sep}Chosen_Elite_To_advance_Gen_{current_gen_int}.smi"
    )
    pass_through_list = Ranking.get_predockcmpds_from_smi_file(
        current_gen_pass_through_smi
    )
    # CHECKED: pass_through_list is of type List[Compound] here.

    # Convert lists to searchable Dictionaries.
    prev_gen_data_dict = Ranking.convert_usable_list_to_lig_dict(prev_gen_data_list)
    # CHECKED: prev_gen_data_dict is of type Dict[str, Compound] here.

    assert prev_gen_data_dict is not None, "prev_gen_data_dict is None"

    pass_through_data: List[Compound] = []
    for lig in pass_through_list:
        # CHECKED: lig is of type Compound here.

        lig_data = prev_gen_data_dict[str(lig.smiles + lig.id)]
        # CHECKED: lig_data is of type Compound here.

        # NOTE: Here it must be converted to a Compound
        assert (
            lig_data.docking_score is not None
        ), "lig_data.previous_docking_score is None"

        # TODO: Nervous that additional_info = "". Not sure what to put there.
        lig_info_remove_diversity_info = Compound(
            smiles=lig.smiles,
            id=lig.id,
            additional_info="",
            docking_score=lig_data.docking_score,
            diversity_score=None,
        )
        pass_through_data.append(lig_info_remove_diversity_info)

    smiles_list.extend(pass_through_data)


def get_predockcmpds_from_smi_file(infile: str) -> List[Compound]:
    """
    Reads and processes a .smi file into a list of Compound objects.

    The .smi file must follow this format for each line:
    SMILES<tab>ID[<tab>optional_info...]<tab>docking_score<tab>diversity_score

    Args:
        infile (str): Path to the formatted .smi file to be read.

    Returns:
        List[Compound]: List of Compound objects with
        information from the .smi file.

    Raises:
        Exception: If the input file does not exist.

    Note:
        - The last two fields are assumed to be docking_score and 
          diversity_score if present.
        - Any fields between ID and docking_score are ignored but allowed.
    """
    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    predock_cmpds: List[Compound] = []

    if os.path.exists(infile) is False:
        print(f"\nFile of Source compounds does not exist: {infile}\n")
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            compoundInfo = Compound.from_tsv_line(line)
            predock_cmpds.append(compoundInfo)

    return predock_cmpds


def convert_usable_list_to_lig_dict(
    predock_cmpds: List[Compound],
) -> Optional[Dict[str, Compound]]:
    """
    Converts a list of Compound objects to a dictionary.

    Args:
        predock_cmpds (List[Compound]): List of Compound
        objects.

    Returns:
        Optional[Dict[str, Compound]]: Dictionary with keys as
        'SMILES+ID' and values as Compound objects. Returns None if
        input is not a list.

    Note:
        If duplicate keys exist, it keeps the entry with the higher docking
        score (lower is better).
    """
    if type(predock_cmpds) is not type([]):
        return None

    usable_dict_of_predock_cmpds: Dict[str, Compound] = {}
    for cmpd in predock_cmpds:
        key = cmpd.smiles + cmpd.id
        if key in usable_dict_of_predock_cmpds and usable_dict_of_predock_cmpds[
            key
        ].get_score_by_type(ScoreType.DOCKING) < cmpd.get_score_by_type(
            ScoreType.DOCKING
        ):
            continue
        usable_dict_of_predock_cmpds[key] = cmpd
    return usable_dict_of_predock_cmpds


##### Called in the docking class ######
def calc_diversity_scores(postDockedCompoundInfos: List[Compound],) -> List[Compound]:
    """
    Calculates diversity scores for a list of Compound objects.

    This function computes Morgan Fingerprints for each molecule and calculates
    pairwise Dice Similarity scores (1.0 meanas a perfect match, 0.0 means no
    match at all). The sum of these scores becomes the diversity score for each
    molecule. Lower scores indicate better diversity, because you want low
    similarity between molecules.

    Args:
        postDockedCompoundInfos (List[Compound]): List of
        Compound objects to process.

    Returns:
        List[Compound]: Input list with updated diversity scores and
        fingerprints.

    Raises:
        AssertionError: If any molecule fails to sanitize or deprotonate.

    Note:
        - Lower diversity scores indicate more diverse molecules.
        - Removes any None entries from the input list.
        - Uses Morgan Fingerprints with radius 10 and feature-based encoding.
    """
    postDockedCompoundInfosToKeep: List[Compound] = []

    for postDockedCompoundInfo in postDockedCompoundInfos:
        if postDockedCompoundInfo is not None:
            smile = postDockedCompoundInfo.smiles
            # name = pair[1]
            try:
                mol = Chem.MolFromSmiles(smile, sanitize=False)
            except Exception:
                mol = None

            # TODO: Should not be assertion here.

            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                    def calc_diversity_scores"
                )

            mol = MOH.check_sanitization(mol)
            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                        def calc_diversity_scores"
                )

            mol = MOH.try_deprotanation(mol)
            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                        def calc_diversity_scores"
                )

            postDockedCompoundInfo.mol = mol
            if mol is None:
                print(postDockedCompoundInfo)
                print("None in temp list, skip this one")
                continue
            # mol_list.append(temp)

            postDockedCompoundInfosToKeep.append(postDockedCompoundInfo)
        else:
            print("noneitem in molecules_list in calc_diversity_scores")

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
