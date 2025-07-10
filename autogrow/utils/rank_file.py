# Files for writing to and loading from a given generation's rank file.

import glob
from typing import List
import os

from autogrow.types import Compound

def _get_rank_file_path(gen_folder_name: str) -> str:
    """
    Get the path to the rank file for a given generation.

    Args:
        gen_folder_name (str): The name of the generation folder.

    Returns:
        str: The path to the rank file.
    """
    return glob.glob(f"{gen_folder_name}{os.sep}*_ranked.smi")[0]

def _get_gen_number_from_rank_filename(rank_file: str) -> str:
    """
    Extract the generation number from the rank file name.
    Args:
        rank_file (str): The path to the rank file.
    Returns:
        str: The generation number extracted from the rank file name.
    """
    return os.path.basename(rank_file).split("_")[1]

def get_gen_number_from_folder_name(gen_folder_name: str) -> str:
    """
    Get the generation number from the folder name.

    Args:
        gen_folder_name (str): The name of the generation folder.

    Returns:
        str: The generation number extracted from the folder name.
    """
    rank_file = _get_rank_file_path(gen_folder_name)
    return _get_gen_number_from_rank_filename(rank_file)

def load_rank_file(gen_folder_name: str) -> List[Compound]:
    """
    Load the rank file for a given generation.

    Args:
        gen_folder_name (str): The name of the generation folder.

    Returns:
        list: A list of lines from the rank file.
    """
    rank_file_path = _get_rank_file_path(gen_folder_name)

    if not os.path.exists(rank_file_path):
        raise FileNotFoundError(f"Rank file not found: {rank_file_path}")

    with open(rank_file_path, "r") as f:
        lines = [l.strip() for l in f.readlines()] if f else []

    cmpd_infos = []
    cmpd_infos.extend(Compound.from_tsv_line(tsv_line) for tsv_line in lines)
    return cmpd_infos

def save_rank_file(rank_file_path: str, cmpds: Compound) -> None:
    """
    Save the rank file for a given generation.

    Args:
        rank_file_path (str): The path to the rank file.
    """
    if not os.path.exists(os.path.dirname(rank_file_path)):
        os.makedirs(os.path.dirname(rank_file_path))

    with open(rank_file_path, "w") as output:
        for cmpdInf in cmpds:
            output.write(cmpdInf.tsv_line)
            # sdf_path = "" if cmpdInf.sdf_path is None else cmpdInf.sdf_path
            # sdf_basename = os.path.basename(sdf_path)
            # output_line = f"{cmpdInf.smiles}\t{cmpdInf.id}\t{cmpdInf.additional_info}\t{cmpdInf.docking_score}\t{cmpdInf.diversity_score}\t{sdf_basename}\n"
            # output.write(output_line)
