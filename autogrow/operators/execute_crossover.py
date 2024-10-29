"""
Performs crossover operations on molecular structures using SMILES.

This module implements a crossover algorithm for molecular structures, using the
Most Common Substructure (MCS) method to identify suitable pairs of molecules
for crossover. It includes functions for selecting molecules, finding MCS, and
performing the actual crossover operation. TODO: Shouldn't be MCS specific here.
"""
import __future__

import random
import copy
from typing import Any, Dict, List, Optional, Tuple, Union, cast

from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.crossover import CrossoverPluginManager
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import LogLevel, log_debug
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import rdFMCS  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


import autogrow.utils.mol_object_handling as MOH

# TODO: Lots of this code is specific to the MCS crossover. But that code should
# all be in the appropriate plugin, not here in this generic function.


def _test_for_mcs(
    params: Dict[str, Any], mol_1: rdkit.Chem.rdchem.Mol, mol_2: rdkit.Chem.rdchem.Mol
) -> Optional[rdkit.Chem.rdFMCS.MCSResult]:
    """
    Finds the Most Common Substructure (MCS) between two molecules.

    Args:
        params (Dict[str, Any]): User parameters governing the MCS search.
        mol_1 (rdkit.Chem.rdchem.Mol): The first RDKit molecule.
        mol_2 (rdkit.Chem.rdchem.Mol): The second RDKit molecule.

    Returns:
        Optional[rdkit.Chem.rdFMCS.MCSResult]: MCS result object if found,
        None otherwise.

    Note:
        - Recommended to use with molecules that have H's removed.
        - Implicit H's are recognized as part of MCS.
    """
    mols = [mol_1, mol_2]
    time_timeout = params["max_time_mcs_prescreen"]
    min_number_atoms_matched = params["min_atom_match_mcs"]

    try:
        result = rdFMCS.FindMCS(
            mols,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            timeout=time_timeout,
        )
    except Exception:
        return None

    # could be used for a theoretical timeout prefilter was to be implement
    # (but this isn't implemented) (ie. if it times out the prefilter dont use
    # in thorough MCS ligmerge) canceled: if True, the MCS calculation did not
    # finish

    # filter by UserDefined minimum number of atoms found. The higher the
    # number the more similar 2 ligands are but the more restrictive for
    # finding mergable ligands number of atoms in common found
    if result.numAtoms < min_number_atoms_matched:
        return None

    return None if result.canceled else result


def _find_sufficiently_similar_cmpd(
    params: Dict[str, Any],
    predock_cmpds: List[PreDockedCompound],
    query_predock_cmpd: PreDockedCompound,
) -> Optional[PreDockedCompound]:
    """
    Selects a random molecule with satisfactory MCS to the given ligand.

    Args:
        params (Dict[str, Any]): User parameters governing the selection.
        predock_cmpds (List[PreDockedCompound]): List of ligands to choose from.
        query_predock_cmpd (PreDockedCompound): The reference (query) ligand.

    Returns:
        Optional[PreDockedCompound]: A suitable second ligand if found,
        None otherwise.
    """
    count = 0
    shuffled_num_list = list(range(len(predock_cmpds) - 1))
    random.shuffle(shuffled_num_list)

    # Convert lig1 into an RDkit mol
    lig1_string = query_predock_cmpd.smiles
    lig1_mol = _convert_mol_from_smiles(lig1_string)

    while count < len(predock_cmpds) - 1:
        rand_num = shuffled_num_list[count]
        mol2_pair = predock_cmpds[rand_num]

        if mol2_pair.smiles == lig1_string:
            count += 1
            continue

        # Convert lig1 into an RDkit mol
        lig_2_string = mol2_pair.smiles
        lig2_mol = _convert_mol_from_smiles(lig_2_string)

        if lig2_mol is None:
            count += 1
            continue

        # it converts and it is not Ligand1. now lets test for a common
        # substructure
        if _test_for_mcs(params, lig1_mol, lig2_mol) is None:
            count += 1
            continue

        # We found a good pair of Ligands
        return mol2_pair

    return None


def _convert_mol_from_smiles(smiles: str) -> Union[rdkit.Chem.rdchem.Mol, bool, None]:
    """
    Converts a SMILES string to an RDKit molecule object.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        Union[rdkit.Chem.rdchem.Mol, bool, None]: RDKit molecule object if 
        conversion is successful, False otherwise.

    Note:
        The function also sanitizes and deprotonates the molecule.
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    except Exception:
        return None

    mol = MOH.check_sanitization(mol)
    if mol is None:
        return None

    mol = MOH.try_deprotanation(mol)
    return False if mol is None else mol


#########################
#### RUN MAIN PARTS #####
#########################


def make_crossovers(
    params: Dict[str, Any],
    generation_num: int,
    number_of_processors: int,
    num_crossovers_to_make: int,
    list_previous_gen_smiles: List[PreDockedCompound],
    new_crossover_smiles_list: List[PreDockedCompound],
) -> Optional[List[PreDockedCompound]]:
    """
    Creates crossover compounds based on previous generation molecules.

    Args:
        params (Dict[str, Any]): User parameters governing the crossover process.
        generation_num (int): Current generation number.
        number_of_processors (int): Number of processors for multithreading.
        num_crossovers_to_make (int): Number of crossovers to generate.
        list_previous_gen_smiles (List[PreDockedCompound]): Molecules from previous generation.
        new_crossover_smiles_list (List[PreDockedCompound]): Previously generated crossovers.

    Returns:
        Optional[List[PreDockedCompound]]: List of new unique ligands, or None if
        generation fails.

    Note:
        Uses multiprocessing to generate crossovers efficiently.
    """
    # TODO: This gets defined again below, so what's the point of defining it
    # here?
    if not new_crossover_smiles_list:
        new_ligands: List[PreDockedCompound] = []
    else:
        new_ligands = copy.deepcopy(new_crossover_smiles_list)

    # Use a temp params dict so you don't put mpi multiprocess info through
    # itself...
    temp_params = {
        key: params[key] for key in list(params.keys()) if key != "parallelizer"
    }
    new_ligands: List[PreDockedCompound] = []
    number_of_processors = int(params["parallelizer"].return_node())

    log_debug("Creating new compounds from selected compounds via crossover")

    with LogLevel():

        loop_counter = 0
        while loop_counter < 2000 and len(new_ligands) < num_crossovers_to_make:

            react_list = copy.deepcopy(list_previous_gen_smiles)

            while len(new_ligands) < num_crossovers_to_make and react_list:

                num_to_grab = num_crossovers_to_make - len(new_ligands)
                num_to_make = num_to_grab

                # to minimize a big loop of running a single crossover at a time
                # we will make 1 new lig/processor. This will help to prevent
                # wasting reasources and time.
                num_to_make = max(num_to_make, number_of_processors)
                smile_pairs = [
                    react_list.pop() for _ in range(num_to_make) if len(react_list) > 0
                ]

                # smile_inputs = [x[0] for x in smile_pairs]
                # smile_names = [x[1] for x in smile_pairs]

                # make a list of tuples for multi-processing Crossover
                job_input: List[
                    Tuple[Dict[str, Any], PreDockedCompound, List[PreDockedCompound],]
                ] = []
                for i in smile_pairs:
                    temp = temp_params, i, list_previous_gen_smiles
                    job_input.append(temp)

                # Example information:
                # result is a list of lists
                # result = [[ligand_new_smiles, lig1_smile_pair,lig2_pair],...]
                # ligand_new_smiles is the smiles string of a new ligand from crossover
                # lig1_smile_pair = ["NCCCCCC","zinc123"]
                # Lig2_smile_pair = ["NCCCO","zinc456"]
                # Lig1 and lig 2 were used to generate the ligand_new_smiles

                results: List[
                    Tuple[str, PreDockedCompound, PreDockedCompound]
                ] = params["parallelizer"].run(
                    tuple(job_input), _do_crossovers_smiles_merge
                )
                results = [x for x in results if x is not None]

                for i in results:
                    if i is None:
                        continue

                    # Get the new molecule's (aka the Child lig) Smile string
                    child_lig_smile = i[0]

                    # get the ID for the parent of a child mol
                    parent_lig1_id = i[1].name
                    parent_lig_2_id = i[2].name

                    # get the unique ID (last few diget ID of the parent mol)
                    parent_lig1_id = parent_lig1_id.split(")")[-1]
                    parent_lig_2_id = parent_lig_2_id.split(")")[-1]

                    # Make a list of all smiles and smile_id's of all previously
                    # made smiles in this generation
                    list_of_already_made_smiles = []
                    list_of_already_made_id = []

                    # fill lists of all smiles and smile_id's of all previously made
                    # smiles in this generation
                    for x in new_ligands:
                        list_of_already_made_smiles.append(x.smiles)
                        list_of_already_made_id.append(x.name)

                    if child_lig_smile not in list_of_already_made_smiles:
                        # if the smiles string is unique to the list of previous
                        # smile strings in this round of reactions then we append it
                        # to the list of newly created ligands we append it with a
                        # unique ID, which also tracks the progress of the reactant
                        is_name_unique = False
                        new_lig_id = None
                        while not is_name_unique:

                            # make unique ID with the 1st number being the
                            # ligand_id_Name for the derived mol. second being the
                            # lig2 number. Followed by Cross. folowed by the
                            # generation number. followed by a  unique.

                            random_id_num = random.randint(100, 1000000)
                            new_lig_id = f"({parent_lig1_id}+{parent_lig_2_id})Gen_{generation_num}_Cross_{random_id_num}"

                            # check name is unique
                            if new_lig_id not in list_of_already_made_id:
                                is_name_unique = True

                        # make a temporary list containing the smiles string of
                        # the new product and the unique ID
                        assert new_lig_id is not None, "new_lig_id is None"
                        ligand_info = PreDockedCompound(
                            smiles=child_lig_smile, name=new_lig_id
                        )

                        # append the new ligand smile and ID to the list of all
                        # newly made ligands
                        new_ligands.append(ligand_info)

            loop_counter += 1

    if len(new_ligands) < num_crossovers_to_make:
        return None

    # once the number of mutants we need is generated return the list
    return new_ligands


def _find_similar_cmpd(
    params: Dict[str, Any],
    ligands_list: List[PreDockedCompound],
    lig1_smile_pair: PreDockedCompound,
) -> Optional[PreDockedCompound]:
    """
    Finds a molecule with sufficient shared structure to the given ligand.

    Args:
        params (Dict[str, Any]): User parameters governing the process.
        ligands_list (List[PreDockedCompound]): List of ligands to choose from.
        ligand1_pair (PreDockedCompound): The reference ligand.

    Returns:
        Optional[PreDockedCompound]: A suitable second ligand if found, None otherwise.
    """
    ligand_1_string = lig1_smile_pair.smiles

    # check if ligand_1 can be converted to an rdkit mol
    lig1 = _convert_mol_from_smiles(ligand_1_string)
    if lig1 is False:
        # Ligand1_string failed to be converted to rdkit mol format
        return None

    # GET TWO UNIQUE LIGANDS TO WITH A SHARED SUBSTRUCTURE
    return _find_sufficiently_similar_cmpd(params, ligands_list, lig1_smile_pair)


def _do_crossovers_smiles_merge(
    params: Dict[str, Any],
    lig1_predock_cmpd: PreDockedCompound,
    all_predock_cmpds: List[PreDockedCompound],
) -> Optional[Tuple[str, PreDockedCompound, PreDockedCompound]]:
    """
    Performs a crossover operation between two ligands.

    Args:
        params (Dict[str, Any]): User parameters governing the process.
        all_predock_cmpds (PreDockedCompound): Information for the first ligand.
        predock_cmpds (List[PreDockedCompound]): List of all seed ligands.

    Returns:
        Optional[Tuple[str, PreDockedCompound, PreDockedCompound]]: 
        Tuple containing new ligand SMILES and parent ligand information,
        or None if crossover fails.

    Note:
        Attempts crossover up to 3 times before giving up.
    """
    # Run the run_smiles_merge_prescreen of the ligand. This gets a new a lig2
    # which passed the prescreen.
    # check if ligand_1 can be converted to an rdkit mol
    lig1_rdkit_mol = _convert_mol_from_smiles(lig1_predock_cmpd.smiles)
    if lig1_rdkit_mol is False:
        # Ligand1_string failed to be converted to rdkit mol format
        lig2_predock_cmpd = None

    # GET A UNIQUE LIGANDS TO WITH A SHARED SUBSTRUCTURE
    lig2_predock_cmpd = _find_sufficiently_similar_cmpd(
        params, all_predock_cmpds, lig1_predock_cmpd
    )

    if lig2_predock_cmpd is None:
        return None

    crossover_manager = cast(CrossoverPluginManager, plugin_managers.Crossover)

    counter = 0
    while counter < 3:
        # run SmilesMerge
        ligand_new_smiles = crossover_manager.run(
            predock_cmpd1=lig1_predock_cmpd, predock_cmpd2=lig2_predock_cmpd
        )

        if ligand_new_smiles is None:
            counter += 1
        else:
            # TODO: Filter accepts a list of PreDockedCompounds. So we need to
            # convert smiles string to that just for the purpose of filtering.
            # This is because crossover doesn't return a PreDockCompound object
            # yet. Need to refactor so that happens. As is, smiles getting
            # converted to PreDockCompound twice.
            tmp_predock_cmpd = PreDockedCompound(smiles=ligand_new_smiles, name="tmp")

            # Filter Here
            pass_or_not = (
                len(plugin_managers.SmilesFilter.run(predock_cmpds=[tmp_predock_cmpd])) > 0
            )

            if not pass_or_not:
                counter += 1
            else:
                return (ligand_new_smiles, lig1_predock_cmpd, lig2_predock_cmpd)
    return None
