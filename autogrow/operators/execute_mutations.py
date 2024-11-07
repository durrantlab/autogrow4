"""
Handles mutation operations for molecular structures.

This module contains functions for creating molecular mutations using various
reaction libraries. It includes functionality for parallel processing of
mutations and handling of mutation results.
"""
import __future__

import random
import copy
from typing import Any, Dict, List, Optional, Tuple, Union


from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.mutation import MutationBase, MutationPluginManager
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import LogLevel, log_debug, log_warning


#######################################
# Functions for creating molecular models
##########################################
def make_mutants(
    params: Dict[str, Any],
    generation_num: int,
    number_of_processors: int,
    num_mutants_to_make: int,
    predock_cmpds: List[PreDockedCompound],
) -> List[PreDockedCompound]:
    """
    Creates mutant compounds based on a list of ligands.

    Args:
        params (Dict[str, Any]): User parameters governing the mutation process.
        generation_num (int): Current generation number.
        number_of_processors (int): Number of processors for parallelization.
        num_mutants_to_make (int): Number of mutants to generate.
        ligands_list (List[PreDockedCompound]): List of ligands to mutate.

    Returns:
        Optional[List[PreDockedCompound]]: List of new mutant ligands, or None if
        generation fails.

    Note:
        Uses multiprocessing to generate mutations efficiently.
        Attempts to create unique mutants, avoiding duplicates.
    """
    new_predock_cmpds: List[PreDockedCompound] = []

    # initialize the smileclickclass
    mutation_plugin_manager = plugin_managers.Mutation
    mutation_plugin_manager.setup_plugins()

    log_debug("Creating new compounds from selected compounds via mutation")

    with LogLevel():
        predock_cmpds_queue = copy.deepcopy(predock_cmpds)
        smiles_already_generated = set([])
        ids_already_generated = set([])
        attempts_to_fill_queue = 0  # to prevent infinite loop

        while len(new_predock_cmpds) < num_mutants_to_make and attempts_to_fill_queue < 5:
            attempts_to_fill_queue += 1
            random.shuffle(predock_cmpds_queue)

            job_input_list = []
            buffer_num = attempts_to_fill_queue * 2  # To avoid zeno's paradox
            for i in range(num_mutants_to_make - len(new_predock_cmpds) + buffer_num):
                imod = i % len(predock_cmpds_queue)
                job_input_list.append(
                    (predock_cmpds_queue[imod], mutation_plugin_manager)
                )
            job_input = tuple(job_input_list)

            results = params["parallelizer"].run(
                job_input, _run_mutation_for_multithread
            )

            # Remove mutants that have already been generated (so unique)
            for i, result in enumerate(results):
                if result is None:
                    continue
                if result[0] in smiles_already_generated:
                    results[i] = None
                smiles_already_generated.add(result[0])

            # Remove None mutant compounds
            results = [x for x in results if x is not None]

            # Add the remaining mutants to the list of new compounds

            for i in results:
                # Get the new molecule's (aka the Child lig) Smile string
                child_lig_smile, reaction_id_number, zinc_id_comp_mol, parent_lig_id = i

                new_lig_id = ""
                while new_lig_id in ids_already_generated or not new_lig_id:
                    # make unique ID with the 1st number being the
                    # parent_lig_id for the derived mol, Followed by
                    # Mutant, folowed by the generationnumber,
                    # followed by a unique.

                    # get the unique ID (last few diget ID of the
                    # parent mol
                    parent_lig_id = parent_lig_id.split(")")[-1]

                    random_id_num = random.randint(100, 1000000)
                    if zinc_id_comp_mol is None:
                        new_lig_id = f"({parent_lig_id})Gen_{generation_num}_Mutant_{reaction_id_number}_{random_id_num}"
                    else:
                        new_lig_id = f"({parent_lig_id}+{zinc_id_comp_mol})Gen_{generation_num}_Mutant_{reaction_id_number}_{random_id_num}"

                ids_already_generated.add(new_lig_id)

                # make a temporary list containing the smiles string
                # of the new product and the unique ID
                ligand_info = PreDockedCompound(
                    smiles=child_lig_smile, name=new_lig_id
                )

                # append the new ligand smile and ID to the list of
                # all newly made ligands
                new_predock_cmpds.append(ligand_info)

    if len(new_predock_cmpds) < num_mutants_to_make:
        log_warning(
            f"Only able to create {len(new_predock_cmpds)} of {num_mutants_to_make} requested mutants."
        )

    # once the number of mutants we need is generated return the list
    return new_predock_cmpds


def _run_mutation_for_multithread(
    predock_cmpd: PreDockedCompound, mutation_obj: MutationPluginManager
) -> Optional[Tuple[str, int, Union[str, None]]]:
    """
    Performs a single mutation operation on a PreDockedCompound.

    This function is designed to be used in a multithreaded context, allowing
    for parallel processing of mutations.

    Args:
        smile (PreDockedCompound): PreDockedCompound of the molecule to mutate.
        mutation_obj (MutationPluginManager): Mutation object to perform the mutation.

    Returns:
        Optional[Tuple[str, int, Union[str, None]]]: Tuple containing the mutated SMILES
        string, reaction ID, and complementary molecule ID (if any), or None if
        the mutation fails.

    Note:
        This function is a wrapper around the mutation object's run method,
        making it suitable for use in multiprocessing contexts.
    """
    resp = mutation_obj.run(predock_cmpd=predock_cmpd)
    return None if resp is None else tuple(list(resp) + [predock_cmpd.name])
