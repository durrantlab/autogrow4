"""
Handles docking and file conversion for molecular docking simulations.

This module provides functions to run docking simulations using various docking
programs and process the results. It includes common operations for different
docking software and manages the workflow from pre-docked compounds to
post-docked, ranked compounds.
"""
import __future__

from typing import Any, Dict, List, cast
from autogrow.docking.ranking.ranking_mol import rank_and_save_output_smi
from autogrow.plugins.registry_base import plugin_managers
from autogrow.plugins.docking import DockingPluginManager
from autogrow.types import Compound
from autogrow.utils.logging import LogLevel, log_info


def run_docking_common(
    cur_gen_dir: str,
    new_gen_predock_cmpds: List[Compound],
) -> str:
    """
    Run common docking operations for all docking programs.

    This function handles the docking process, including running the docking
    simulation, processing results, and ranking the docked compounds.

    Args:
        cur_gen_dir (str): Directory for the current generation.
        new_gen_predock_cmpds (List[Compound]): List of Compound
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
        post_docked_compounds = post_docked_compounds

    # Remove those that failed to convert
    post_docked_compounds = [x for x in post_docked_compounds if x is not None]

    # Remove those not associated with a docked sdf file
    return [x for x in post_docked_compounds if x.sdf_path is not None and x.docking_score is not None]


