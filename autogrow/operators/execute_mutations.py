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
from autogrow.plugins.mutation import MutationBase
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import LogLevel, log_debug


#######################################
# Functions for creating molecular models
##########################################
def make_mutants(
    params: Dict[str, Any],
    generation_num: int,
    number_of_processors: int,
    num_mutants_to_make: int,
    ligands_list: List[PreDockedCompound],
    new_mutation_smiles_list: List[PreDockedCompound],
) -> Optional[List[PreDockedCompound]]:
    """
    Creates mutant compounds based on a list of ligands.

    Args:
        params (Dict[str, Any]): User parameters governing the mutation process.
        generation_num (int): Current generation number.
        number_of_processors (int): Number of processors for parallelization.
        num_mutants_to_make (int): Number of mutants to generate.
        ligands_list (List[PreDockedCompound]): List of ligands to mutate.
        new_mutation_smiles_list (List[PreDockedCompound]): Previously generated mutants.

    Returns:
        Optional[List[PreDockedCompound]]: List of new mutant ligands, or None if
        generation fails.

    Note:
        Uses multiprocessing to generate mutations efficiently.
        Attempts to create unique mutants, avoiding duplicates.
    """
    new_ligands_list: List[PreDockedCompound] = new_mutation_smiles_list or []
    loop_counter = 0

    # TODO: What is this???
    number_of_processors = int(params["parallelizer"].return_node())

    # initialize the smileclickclass
    mutation_plugin_manager = plugin_managers.Mutation
    mutation_plugin_manager.setup_plugins()

    log_debug("Creating new compounds from selected compounds via mutation")

    with LogLevel():

        # SmileClickClass.SmilesClickChem(
        #     rxn_library_variables, new_mutation_smiles_list
        # )

        while loop_counter < 2000 and len(new_ligands_list) < num_mutants_to_make:

            react_list = copy.deepcopy(ligands_list)

            while len(new_ligands_list) < num_mutants_to_make and react_list:

                # mutation_plugin_manager.add_mutant_smiles(new_ligands_list)

                num_to_grab = num_mutants_to_make - len(new_ligands_list)
                num_to_make = num_to_grab

                # to minimize a big loop of running a single mutation at a time we
                # will make 1 new lig/processor. This will help to prevent wasting
                # reasources and time.
                num_to_make = max(num_to_make, number_of_processors)
                smile_pairs = [
                    react_list.pop() for _ in range(num_to_make) if len(react_list) > 0
                ]

                smile_inputs = [x.smiles for x in smile_pairs]
                smile_names = [x.name for x in smile_pairs]

                job_input = tuple(
                    (smile, mutation_plugin_manager) for smile in smile_inputs
                )

                results = params["parallelizer"].run(
                    job_input, _run_smiles_click_for_multithread
                )

                for index, i in enumerate(results):
                    if i is None:
                        continue

                    # Get the new molecule's (aka the Child lig) Smile string
                    child_lig_smile = i[0]

                    # get the ID for the parent of a child mol and the
                    # complementary parent mol. comp mol could be None or a
                    # zinc database ID
                    parent_lig_id = smile_names[index]
                    # Make a list of all smiles and smile_id's of all
                    # previously made smiles in this generation
                    list_of_already_made_smiles = []
                    list_of_already_made_id = []

                    # fill lists of all smiles and smile_id's of all
                    # previously made smiles in this generation
                    for x in new_ligands_list:
                        list_of_already_made_smiles.append(x.smiles)
                        list_of_already_made_id.append(x.name)

                    if child_lig_smile not in list_of_already_made_smiles:
                        # if the smiles string is unique to the list of
                        # previous smile strings in this round of reactions
                        # then we append it to the list of newly created
                        # ligands we append it with a unique ID, which also
                        # tracks the progress of the reactant
                        is_name_unique = False
                        new_lig_id = ""

                        # get the reaction id number
                        reaction_id_number = i[1]

                        zinc_id_comp_mol = i[2]

                        while not is_name_unique:
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

                            # check name is unique
                            if new_lig_id not in list_of_already_made_id:
                                is_name_unique = True

                        # make a temporary list containing the smiles string
                        # of the new product and the unique ID
                        ligand_info = PreDockedCompound(
                            smiles=child_lig_smile, name=new_lig_id
                        )

                        # append the new ligand smile and ID to the list of
                        # all newly made ligands
                        new_ligands_list.append(ligand_info)

            loop_counter += 1

    if len(new_ligands_list) < num_mutants_to_make:
        return None

    # once the number of mutants we need is generated return the list
    return new_ligands_list


def _run_smiles_click_for_multithread(
    smile: str, mutation_obj: MutationBase
) -> Optional[List[Union[str, int, None]]]:
    """
    Performs a single mutation operation on a SMILES string.

    This function is designed to be used in a multithreaded context, allowing
    for parallel processing of mutations.

    Args:
        smile (str): SMILES string of the molecule to mutate.
        mutation_obj (MutationBase): Mutation object to perform the mutation.

    Returns:
        Optional[List[Union[str, int, None]]]: List containing the mutated SMILES
        string, reaction ID, and complementary molecule ID (if any), or None if
        the mutation fails.

    Note:
        This function is a wrapper around the mutation object's run method,
        making it suitable for use in multiprocessing contexts.
    """
    return mutation_obj.run(parent_smiles=smile)
