"""
Handles docking and file conversion for molecular docking simulations.

This module provides functions to run docking simulations using various docking
programs and process the results. It includes common operations for different
docking software and manages the workflow from pre-docked compounds to
post-docked, ranked compounds.
"""
import __future__

from typing import Any, Dict, List, cast
from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.docking import DockingPluginManager
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import LogLevel, log_info

def run_docking_common(
    params: Dict[str, Any],  # TODO: Not used.
    current_gen_int: int,
    cur_gen_dir: str,
    smiles_file_new_gen: str,
    new_gen_predock_cmpds: List[PreDockedCompound],
) -> str:
    """
    Runs common docking operations for all docking programs.

    This function handles the docking process, including running the docking
    simulation, processing results, and ranking the docked compounds.

    Args:
        params (Dict[str, Any]): User variables governing program execution.
        current_gen_int (int): Current generation number.
        cur_gen_dir (str): Directory for the current generation.
        smiles_file_new_gen (str): Filename containing new population molecules.
        new_gen_predock_cmpds (List[PreDockedCompound]): List of PreDockedCompound
            objects for the new generation.

    Returns:
        str: Filename of the unweighted-ranked SMILES with their docking scores.
    """

    docking_plugin_manager = cast(DockingPluginManager, plugin_managers.Docking)

    log_info("Starting docking")
    with LogLevel():
        post_docked_compounds = docking_plugin_manager.run(
            predocked_cmpds=new_gen_predock_cmpds, cache_dir=cur_gen_dir
        )

    # Remove those that failed to convert
    post_docked_compounds = [x for x in post_docked_compounds if x is not None]

    # Remove those not associated with a docked sdf file
    post_docked_compounds = [
        x for x in post_docked_compounds if x.docked_sdf_path is not None
    ]

    print("\nBegin Ranking and Saving results")
    unweighted_ranked_smile_file = docking_plugin_manager.rank_and_save_output_smi(
        cur_gen_dir,
        current_gen_int,
        smiles_file_new_gen,
        post_docked_compounds,
    )
    print("\nCompleted Ranking and Saving results\n")
    return unweighted_ranked_smile_file
