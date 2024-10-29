"""
Populates an AutoGrow generation via mutation, crossover, and elitism.

This module handles the creation of a new generation of compounds in the 
AutoGrow evolutionary algorithm. It includes functions for generating mutations
and crossovers, selecting elite compounds, and managing the overall population
generation process. It also handles filtering and conversion of SMILES to 3D SDFs.
"""
import __future__

import os
import random
import copy
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.types import PreDockedCompound
from autogrow.utils.caching import CacheManager
from autogrow.utils.logging import LogLevel, log_info, log_warning
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.docking.ranking.ranking_mol as Ranking
import autogrow.operators.execute_mutations as Mutation
import autogrow.operators.execute_crossover as execute_crossover
import autogrow.utils.mol_object_handling as MOH

#############
# Main run Autogrow operators to make a generation
#############


def populate_generation(
    params: Dict[str, Any], generation_num: int, cur_gen_dir: str
) -> Tuple[str, List[PreDockedCompound]]:
    """
    Populates a new generation of ligands through mutation, crossover, and elitism.

    This function orchestrates the entire process of creating a new generation,
    including mutation, crossover, elite selection, and file management.

    Args:
        params (Dict[str, Any]): Dictionary of all user variables.
        generation_num (int): The current generation number.
        cur_gen_dir (str): Directory for the current generation.

    Returns:
        Tuple[str, List[PreDockedCompound]]: A tuple containing:
            - The name of the .smi file containing the new population.
            - A list of PreDockedCompound objects representing the new population.

    Raises:
        AssertionError: If the population fails to make enough compounds.
    """
    number_of_processors = int(params["number_of_processors"])

    # Determine which generation it is and how many mutations and crossovers to make
    if generation_num == 1:
        # If 1st generation
        num_crossovers = params["number_of_crossovers_first_generation"]
        num_mutations = params["number_of_mutants_first_generation"]

        # How many advance from previous generation to the next generation directly
        # This will be done later but we are unpacking params here
        num_elite_prev_gen = params[
            "number_elitism_advance_from_previous_gen_first_generation"
        ]
    else:
        # Later generations
        num_crossovers = params["number_of_crossovers"]
        num_mutations = params["number_of_mutants"]
        num_elite_prev_gen = params["number_elitism_advance_from_previous_gen"]

    # Get the source compound list
    src_cmpds = _get_cmpds_prev_gen(params, generation_num)

    num_seed_diversity, num_seed_dock_fitness = _get_seed_pop_sizes(
        params, generation_num
    )

    # Total population size of this generation
    total_num_desired_new_ligands = num_crossovers + num_mutations + num_elite_prev_gen

    # Generate mutations
    log_info("Creating mutant compounds")
    with LogLevel():
        with CacheManager("MutantCompounds", cur_gen_dir) as cache:
            if cache.exists:
                mut_predock_cmpds = cache.data
            else:
                if num_mutations > 0:
                    mut_predock_cmpds = _generate_mutations(
                        params,
                        generation_num,
                        num_mutations,
                        num_seed_diversity,
                        num_seed_dock_fitness,
                        src_cmpds,
                        number_of_processors,
                        cur_gen_dir,
                    )
                else:
                    log_warning("No mutations made, per user settings")
                    mut_predock_cmpds: List[PreDockedCompound] = []
                cache.data = mut_predock_cmpds
        log_info(f"Created {len(mut_predock_cmpds)} mutant compounds")

    # Generate crossovers
    log_info("Creating crossover compounds")

    with LogLevel():
        with CacheManager("CrossoverCompounds", cur_gen_dir) as cache:
            if cache.exists:
                cross_predock_cmpds = cache.data
            else:
                if num_crossovers > 0:
                    cross_predock_cmpds = _generate_crossovers(
                        params,
                        generation_num,
                        num_crossovers,
                        num_seed_diversity,
                        num_seed_dock_fitness,
                        src_cmpds,
                        number_of_processors,
                    )
                else:
                    log_warning("No crossovers made, per user settings")
                    cross_predock_cmpds: List[PreDockedCompound] = []
                cache.data = cross_predock_cmpds

        log_info(f"Created {len(cross_predock_cmpds)} crossover compounds")

    # Get ligands from previous generation
    log_info("Identifying elite compounds to advance, without mutation or crossover")
    with LogLevel():
        with CacheManager("EliteCompounds", cur_gen_dir) as cache:
            if cache.exists:
                elite_predock_cmpds = cache.data
            else:
                # No cache, so we need to generate the elite compounds
                if num_elite_prev_gen > 0:
                    elite_predock_cmpds = _get_elite_cmpds_prev_gen(
                        params, src_cmpds, num_elite_prev_gen, generation_num,
                    )
                else:
                    log_warning("No elite ligands advanced, per user settings")
                    elite_predock_cmpds: List[PreDockedCompound] = []

                # Saves to cache
                cache.data = elite_predock_cmpds
        log_info(f"Identified {len(elite_predock_cmpds)} elite compounds")

    # Build new_gen_smis and full_gen_smis
    # TODO: Need to understand why these two are separate.
    new_gen_predock_cmpds: List[
        PreDockedCompound
    ] = mut_predock_cmpds + cross_predock_cmpds

    full_gen_predock_cmpds: List[
        PreDockedCompound
    ] = mut_predock_cmpds + cross_predock_cmpds

    # TODO: Consider depreciating redock_elite_from_previous_gen
    if params["redock_elite_from_previous_gen"] is False and generation_num != 1:
        # Doesn't append to the new_generation_smiles_list
        full_gen_predock_cmpds.extend(iter(elite_predock_cmpds))
    else:
        for i in elite_predock_cmpds:
            new_gen_predock_cmpds.append(i)
            full_gen_predock_cmpds.append(i)

    if len(full_gen_predock_cmpds) < total_num_desired_new_ligands:
        print("We needed ", total_num_desired_new_ligands)
        print("We made ", len(full_gen_predock_cmpds))
        print(
            "Population failed to make enough mutants or crossovers... "
            "Errors could include not enough diversity, too few seeds to "
            "the generation, the seed mols are unable to cross-over due "
            "to lack of similarity, or all of the seed lack functional groups "
            "for performing reactions"
        )
        assert False

    # Save the full generation and the SMILES to convert
    (
        full_gen_smi_file,
        smiles_to_convert_file,
        new_gen_folder_path,
    ) = _save_smiles_files(
        params,
        generation_num,
        full_gen_predock_cmpds,
        new_gen_predock_cmpds,
        "_to_convert",
    )

    # Convert SMILES to .sdf using Gypsum and convert .sdf to .pdb with RDKit

    # Note that smiles_to_convert_file is a single with with many smi files
    # in it.
    log_info("Converting SMILES to 3D SDF files")
    with LogLevel():
        full_gen_predock_cmpds = plugin_managers.SmiTo3DSdf.run(
            predock_cmpds=full_gen_predock_cmpds,
            pwd=new_gen_folder_path,
            cache_dir=cur_gen_dir,
        )

    # Remove those that failed to convert
    full_gen_predock_cmpds = [
        x for x in full_gen_predock_cmpds if x.sdf_3d_path is not None
    ]

    return full_gen_smi_file, full_gen_predock_cmpds


def _generate_mutations(
    params: Dict[str, Any],
    generation_num: int,
    num_mutations: int,
    num_seed_diversity: int,
    num_seed_dock_fitness: int,
    src_cmpds: List[PreDockedCompound],
    number_of_processors: int,
    cur_gen_dir: str,
) -> List[PreDockedCompound]:
    """
    Generates mutations for the current generation.

    This function creates mutant compounds based on the seed list from the
    previous generation, using the specified reaction library.

    Args:
        params (Dict[str, Any]): User parameters governing the mutation process.
        generation_num (int): Current generation number.
        num_mutations (int): Number of mutations to generate.
        num_seed_diversity (int): Number of seed molecules chosen for diversity.
        num_seed_dock_fitness (int): Number of seed molecules chosen for docking fitness.
        src_cmpds (List[PreDockedCompound]): Source compounds from previous generation.
        number_of_processors (int): Number of processors for parallel processing.
        cur_gen_dir (str): Current generation directory.

    Returns:
        List[PreDockedCompound]: List of newly generated mutant compounds.

    Raises:
        Exception: If insufficient mutants are generated.
    """
    # Get starting compounds for Mutations
    seed_list_mutations = _make_seed_list(
        src_cmpds, generation_num, num_seed_diversity, num_seed_dock_fitness,
    )

    # Save seed list for Mutations
    _save_ligand_list(
        params["output_directory"],
        generation_num,
        seed_list_mutations,
        "Mutation_Seed_List",
    )

    # Package user params specifying the Reaction library to use for mutation
    rxn_library_variables = [
        params["rxn_library_path"],
    ]

    # List of SMILES from mutation
    new_mutation_smiles_list: List[PreDockedCompound] = []

    # Make all the required ligands by mutations
    while len(new_mutation_smiles_list) < num_mutations:
        num_mutants_to_make = num_mutations - len(new_mutation_smiles_list)

        # Make all mutants
        new_mutants = Mutation.make_mutants(
            params,
            generation_num,
            number_of_processors,
            num_mutants_to_make,
            seed_list_mutations,
            new_mutation_smiles_list,
        )

        if new_mutants is None:
            break

        # Remove Nones
        new_mutants = [x for x in new_mutants if x is not None]

        for i in new_mutants:
            new_mutation_smiles_list.append(i)
            if len(new_mutation_smiles_list) == num_mutations:
                break

        # Save new_mutation_smiles_list
        _save_ligand_list(
            params["output_directory"],
            generation_num,
            new_mutation_smiles_list,
            "Chosen_Mutants",
        )

        if (
            new_mutation_smiles_list is None
            or len(new_mutation_smiles_list) < num_mutations
        ):
            _throw_error_not_enough_compounds_made(
                num_mutations,
                " ligands through Mutation",
                new_mutation_smiles_list,
                "Mutation failed to make enough new ligands.",
            )
        # print("FINISHED MAKING MUTATIONS")

    return new_mutation_smiles_list


def _generate_crossovers(
    params: Dict[str, Any],
    generation_num: int,
    num_crossovers: int,
    num_seed_diversity: int,
    num_seed_dock_fitness: int,
    src_cmpds: List[PreDockedCompound],
    number_of_processors: int,
) -> List[PreDockedCompound]:
    """
    Generates crossovers for the current generation.

    This function creates crossover compounds based on the seed list from the
    previous generation.

    Args:
        params (Dict[str, Any]): User parameters governing the crossover process.
        generation_num (int): Current generation number.
        num_crossovers (int): Number of crossovers to generate.
        num_seed_diversity (int): Number of seed molecules chosen for diversity.
        num_seed_dock_fitness (int): Number of seed molecules chosen for docking fitness.
        src_cmpds (List[PreDockedCompound]): Source compounds from previous generation.
        number_of_processors (int): Number of processors for parallel processing.

    Returns:
        List[PreDockedCompound]: List of newly generated crossover compounds.

    Raises:
        Exception: If insufficient crossovers are generated.
    """
    # Get starting compounds to seed Crossovers
    seed_list_crossovers = _make_seed_list(
        src_cmpds, generation_num, num_seed_diversity, num_seed_dock_fitness,
    )

    # Save seed list for Crossovers
    _save_ligand_list(
        params["output_directory"],
        generation_num,
        seed_list_crossovers,
        "Crossover_Seed_List",
    )

    # print("MAKE CROSSOVERS")

    # Making Crossovers
    # List of smiles from crossover
    new_crossover_smiles_list: List[PreDockedCompound] = []

    # Make all the required ligands by Crossover
    while len(new_crossover_smiles_list) < num_crossovers:
        num_crossovers_to_make = num_crossovers - len(new_crossover_smiles_list)

        # Make all crossovers
        new_crossovers = execute_crossover.make_crossovers(
            params,
            generation_num,
            number_of_processors,
            num_crossovers_to_make,
            seed_list_crossovers,
            new_crossover_smiles_list,
        )

        if new_crossovers is None:
            break

        # Remove Nones
        new_crossovers = [x for x in new_crossovers if x is not None]

        # Append those which passed the filter
        for i in new_crossovers:
            new_crossover_smiles_list.append(i)
            if len(new_crossover_smiles_list) == num_crossovers:
                break

    # Save new_crossover_smiles_list
    _save_ligand_list(
        params["output_directory"],
        generation_num,
        new_crossover_smiles_list,
        "Chosen_Crossovers",
    )

    if (
        new_crossover_smiles_list is None
        or len(new_crossover_smiles_list) < num_crossovers
    ):
        _throw_error_not_enough_compounds_made(
            num_crossovers,
            " ligands through Crossover",
            new_crossover_smiles_list,
            "Crossover failed to make enough new ligands.",
        )
    # print("FINISHED MAKING CROSSOVERS")

    return new_crossover_smiles_list


def _get_elite_cmpds_prev_gen(
    params: Dict[str, Any],
    src_cmpds: List[PreDockedCompound],
    num_elite_prev_gen: int,
    generation_num: int,
) -> List[PreDockedCompound]:
    """
    Selects elite compounds from the previous generation to advance.

    This function handles the selection of top-performing compounds from the
    previous generation to be carried forward without modification.

    Args:
        params (Dict[str, Any]): User parameters.
        src_cmpds (List[PreDockedCompound]): Source compounds from previous generation.
        num_elite_prev_gen (int): Number of elite compounds to select.
        generation_num (int): Current generation number.

    Returns:
        List[PreDockedCompound]: List of selected elite compounds.

    Raises:
        Exception: If selection process fails or insufficient compounds are available.
    """
    chosen_mol_to_pass_through_list = _make_pass_through_list(
        params, src_cmpds, num_elite_prev_gen, generation_num
    )

    if isinstance(chosen_mol_to_pass_through_list, str):
        printout = (
            chosen_mol_to_pass_through_list
            + "\nIf this is the 1st generation, it may be due to the starting "
            + "library having SMILES which could not be converted to sanitizable "
            + "RDKit Molecules"
        )
        raise Exception(printout)

    assert isinstance(
        chosen_mol_to_pass_through_list, list
    ), "Chosen mol to pass through list is not a list"

    # Save chosen_mol_to_pass_through_list
    _save_ligand_list(
        params["output_directory"],
        generation_num,
        chosen_mol_to_pass_through_list,
        "Chosen_Elite_To_advance",
    )

    return chosen_mol_to_pass_through_list


def _save_smiles_files(
    params: Dict[str, Any],
    generation_num: int,
    full_gen_smis: List[PreDockedCompound],
    new_gen_smis: List[PreDockedCompound],
    suffix: Optional[str],
) -> Tuple[str, str, str]:
    """
    Save SMILES files for the current generation.

    This function saves the full generation and the SMILES file to be 
    converted to 3D structures.

    Args:
        params (Dict[str, Any]): User parameters for output paths.
        generation_num (int): Current generation number.
        full_gen_smis (List[PreDockedCompound]): Full list of compounds for 
            the generation.
        new_gen_smis (List[PreDockedCompound]): List of new compounds for 
            conversion.
        suffix (Optional[str]): Suffix to append to the file name.

    Returns:
        Tuple[str, str, str]: A tuple containing:
            - The full generation SMILES file path.
            - The SMILES-to-convert file path.
            - The folder path for the current generation.
    """
    # Save the full generation
    full_generation_smiles_file, new_gen_folder_path = _save_generation_smi(
        params["output_directory"], generation_num, full_gen_smis, None
    )

    # Save the file to convert to 3D
    smiles_to_convert_file, _ = _save_generation_smi(
        params["output_directory"], generation_num, new_gen_smis, suffix,
    )

    return full_generation_smiles_file, smiles_to_convert_file, new_gen_folder_path


def _throw_error_not_enough_compounds_made(arg0, arg1, arg2, arg3):
    """
    Raise an error if insufficient compounds are generated.

    This function prints details about the number of required and generated 
    compounds, then raises an exception if the generated compounds are 
    fewer than required.

    Args:
        arg0 (int): Number of compounds that should have been generated.
        arg1 (str): Descriptor of the compounds (e.g., 'mutants' or 'crossovers').
        arg2 (List[PreDockedCompound]): List of generated compounds.
        arg3 (str): Error message to raise if there are insufficient compounds.

    Raises:
        Exception: Raised with the provided error message.
    """
    print("")
    print("")
    print(f"We needed to make {arg0}{arg1}")
    print(f"We only made {len(arg2)}{arg1}")
    print("")
    print("")
    raise Exception(arg3)


#############
# Get seeds
#############
def _test_source_smiles_convert(
    smile_info: PreDockedCompound,
) -> Union[PreDockedCompound, str]:
    """
    Attempts to convert a SMILES string to an rdkit.Chem.rdchem.Mol object.

    This function is done in a try statement to handle bad SMILES strings that
    are incapable of being converted. It also checks that the SMILES string is
    able to be sanitized.

    Args:
        smile_info (PreDockedCompound): A PreDockedCompound object containing
            the SMILES of a ligand, its ID, and potentially additional
            information about the ligand.

    Returns:
        Union[PreDockedCompound, str]: If it passed the test, it returns the
            PreDockedCompound object. If it failed to convert, it returns an
            error message string. This passes out to prevent MPI print issues.
    """
    if smile_info is None or not smile_info:
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: Blank "
            + "entry in source compound list.\n"
        )
        return f"{printout}\tRemoving: {smile_info}"

    # TODO: JDD: What was this? Never understood, so commented out
    # if len(smile_info) == 1:
    #     printout = (
    #         "REMOVING SMILES FROM SOURCE LIST: Unformatted or blank "
    #         + "entry in source compound list.\n"
    #     )
    #     printout += f"\tRemoving: {smile_info}"
    #     return printout

    # separate out SMILES str and ID
    smile_str = smile_info.smiles
    smile_id = smile_info.name

    if type(smile_str) is not type(""):
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: SMILES string is not a "
            + "String. Check for formatting errors. \n"
        )
        printout += f"\tIgnored SMILES is: {smile_str}"
        return printout

    # Try importing it into RDKit with Sanitization off. Tests for errors in
    # having the wrong data type
    try:
        mol = Chem.MolFromSmiles(str(smile_str), sanitize=False)
    except Exception:
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: SMILES string failed "
            + "to import into RDKit.\n\t "
        )
        printout += f"Removed SMILE string is: {smile_str} \n"
        printout += f"\t Removed SMILE ID is: {smile_id}"
        return printout

    # This will fail if there are valence errors. We won't try to correct
    # someones source compound list Although the MOH.check_sanitization will
    # do that. try sanitizing, which is necessary later
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: SMILES "
            + "string failed to Sanitize in RDKit.\n"
        )
        printout += f"\t Removed SMILE string is: {smile_str} \n"
        printout += f"\t Removed SMILE ID is: {smile_id}"
        return printout

    # Make the mol again fresh and try running it through MOH.handleHs() This
    # will try protanating and Deprotanating the mol. If it can't handle that
    # We reject it as many functions will require this sort of manipulation.
    # More advanced sanitization issues will also be removed in this step
    mol = Chem.MolFromSmiles(str(smile_str), sanitize=False)
    mol = MOH.handleHs(mol, True)

    if mol is None:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string failed \
                    to be protanated or deprotanated.\n"
        printout = (
            printout
            + "\t This is often an issue with valence and sanitization "
            + "issues with the SMILES string."
        )
        return _report_removed_compound_info(smile_str, printout, smile_id)
    # Check there are no * which are atoms with atomic number=0
    mol = MOH.check_for_unassigned_atom(mol)
    if mol is None:
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: SMILES string contained "
            + "an unassigned atom type labeled as *.\n"
        )
        return _report_removed_compound_info(smile_str, printout, smile_id)
    # Check for fragments.
    if len(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)) != 1:

        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string was fragmented.\n"
        return _report_removed_compound_info(smile_str, printout, smile_id)
    # the ligand is good enough to use throughout the program!
    return smile_info


def _report_removed_compound_info(smile_str, printout, smile_id):
    """
    Create a report for a removed compound.

    This function formats the information about a removed SMILES string 
    and its ID for logging purposes.

    Args:
        smile_str (str): SMILES string of the compound.
        printout (str): Initial message to include in the report.
        smile_id (str): Identifier for the compound.

    Returns:
        str: Formatted string describing the removed compound.
    """
    printout += f"\t Removed SMILE string is: {smile_str} \n"
    printout += f"\t Removed SMILE ID is: {smile_id}"
    return printout


def _get_cmpds_prev_gen(
    params: Dict[str, Any], generation_num: int
) -> List[PreDockedCompound]:
    """
    Get the source compounds list from the previous generation of the source
    compound list

    This also filters the list to ensure mols can be imported to RDKit and
    that they pass the drug-likliness filters.

    If generation = 1 use the User specified starting compounds If generation
    is >1 than use the previous generations top ligands. This takes an .smi
    file

    Inputs:
    :param dict params: a dictionary of all user variables
    :param int generation_num: the interger of the current generation

    Returns:
    :returns: list predock_cmpds: a list with SMILES strings, names,
        and information about the smiles from the previous generation or the
        source compound list
    """
    source_file_gen_0 = (
        f"{params['output_directory']}generation_0{os.sep}generation_0_ranked.smi"
    )
    if generation_num == 0:
        predock_cmpds = _get_source_compounds_or_raise(params)
    elif generation_num == 1 and os.path.exists(source_file_gen_0) is False:
        predock_cmpds = _get_source_compounds_or_raise(params)
    else:
        source_file = (
            params["output_directory"]
            + f"generation_{generation_num - 1}{os.sep}generation_{generation_num - 1}_ranked.smi"
        )
        if os.path.exists(source_file) is False:
            _handle_no_ligands_found("\tCheck formatting or if file has been moved.\n")
        predock_cmpds = Ranking.get_predockcmpds_from_smi_file(source_file)

        if len(predock_cmpds) == 0:
            _handle_no_ligands_found("\tCheck formatting or if file has been moved. \n")
    # Test that every SMILES in the predock_cmpds is a valid SMILES
    # which will import and Sanitize in RDKit. SMILES will be excluded if they
    # are fragmented, contain atoms with no atomic number (*), or do not
    # sanitize
    job_input = tuple((i,) for i in predock_cmpds)

    predock_cmpds: List[PreDockedCompound] = params["parallelizer"].run(
        job_input, _test_source_smiles_convert
    )
    predock_cmpds = [x for x in predock_cmpds if x is not None]
    print_errors = [x for x in predock_cmpds if type(x) is str]
    predock_cmpds = [x for x in predock_cmpds if type(x) is PreDockedCompound]
    for x in print_errors:
        print(x)

    if not predock_cmpds:
        _raise_exception_with_message(
            "\nThere were no ligands in source compound or previous \
            generation which could sanitize.\n"
        )

    random.shuffle(predock_cmpds)

    return predock_cmpds


def _raise_exception_with_message(arg0):
    """
    Raise an exception with a specific message.

    This function prints the provided error message and then raises it as an 
    exception.

    Args:
        arg0 (str): Error message to display and raise.

    Raises:
        Exception: Raised with the provided error message.
    """
    printout = arg0
    print(printout)
    raise Exception(printout)


def _get_source_compounds_or_raise(params) -> List[PreDockedCompound]:
    """
    Retrieve source compounds or raise an error if none are found.

    This function loads the source compounds specified in the user parameters 
    and raises an exception if the list is empty.

    Args:
        params (Dict[str, Any]): User parameters containing the source 
            compound file path.

    Returns:
        List[PreDockedCompound]: List of source compounds.

    Raises:
        Exception: If no ligands are found in the source compound file.
    """
    # This will be the full length list of starting molecules as the seed
    source_file = str(params["source_compound_file"])
    result = Ranking.get_predockcmpds_from_smi_file(source_file)

    if len(result) == 0:
        print(
            "\nThere were no available ligands in source compound. Check formatting\n"
        )
        raise Exception(
            "There were no available ligands in source compound. Check formatting"
        )

    return result


def _handle_no_ligands_found(arg0):
    """
    Handle the case where no ligands are found.

    This function logs and raises an error when no ligands are available in 
    the previous generation's ranked file.

    Args:
        arg0 (str): Additional information about the error.

    Raises:
        Exception: Raised with a message indicating no ligands were found.
    """
    printout = (
        "\n"
        + "There were no available ligands in previous"
        + " generation ranked ligand file.\n"
    )
    printout = printout + arg0
    print(printout)
    raise Exception(printout)


def _make_seed_list(
    source_compounds_list: List[PreDockedCompound],
    generation_num: int,
    num_seed_diversity: int,
    num_seed_dock_fitness: int,
) -> List[PreDockedCompound]:
    """
    Get the starting compound list for running the Mutation and Crossovers.

    Args:
        source_compounds_list (List[PreDockedCompound]): A list with SMILES
            strings, names, and information about the smiles from either the
            previous generation or the source compound list.
        generation_num (int): The integer of the current generation.
        num_seed_diversity (int): The number of seed molecules which come
            from diversity selection.
        num_seed_dock_fitness (int): The number of seed molecules which come
            from elite selection by docking score.

    Returns:
        List[PreDockedCompound]: A list with SMILES strings, names, and
            information about the smiles which will be used to seed the next
            generation.
    """
    predock_cmpds: List[PreDockedCompound] = copy.deepcopy(source_compounds_list)

    # Code no longer uses generation_num == 0, but just to make sure...
    assert generation_num >= 0, "Generation number must be greater than or equal to 0"

    if generation_num == 1:
        # If generation 1, then we don't need to do anything. All source compounds
        # are available to participate in mutations and crossovers.
        log_info(f"Selecting all {len(predock_cmpds)} source compounds (generation 1)")
    else:
        # selector_choice = params["selector_choice"]
        # tourn_size = params["tourn_size"]
        # Get subset of the source_file based on diversity scores and docking
        # scores

        # This assumes the most negative number is the best option which is true
        # for both. This is true for both the diversity score and the docking
        # score, because the docking plugin should return more negative scores
        # for better compounds.
        predock_cmpds = plugin_managers.Selector.run(
            predock_cmpds=predock_cmpds,
            num_seed_dock_fitness=num_seed_dock_fitness,
            num_seed_diversity=num_seed_diversity,
        )

    random.shuffle(predock_cmpds)

    return predock_cmpds


def _get_seed_pop_sizes(params: Dict[str, Any], gen_num: int) -> Tuple[int, int]:
    """
    Determines how many molecules will be chosen to seed a generation based on
    their docking score and diversity score.

    Args:
        params (Dict[str, Any]): A dictionary of all user variables.
        gen_num (int): The integer of the current generation.

    Returns:
        Tuple[int, int]: A tuple containing:
            - The number of seed molecules which come from diversity selection.
            - The number of seed molecules which come from elite selection by
              docking score.
    """
    # How many fewer seed mols are chosen from diversity compared to the 1st
    # generation This is also how many more get chosen from elitist selection
    diversity_depreciation = (
        int(gen_num - 1) * params["diversity_seed_depreciation_per_gen"]
    )

    if gen_num == 1:
        top_mols_to_seed_next_gen = params[
            "top_mols_to_seed_next_generation_first_generation"
        ]

    else:
        top_mols_to_seed_next_gen = params["top_mols_to_seed_next_generation"]

    # Number of mols chosen because of their diversity score
    num_seed_diversity = (
        params["diversity_mols_to_seed_first_generation"] - diversity_depreciation
    )

    # Number of mols chosen because of their docking score. Elitist style
    # selection but this will be chosen using a weight roulette selector later
    # on.
    if num_seed_diversity <= 0:
        num_seed_dock_fitness = (
            top_mols_to_seed_next_gen
            + params["diversity_mols_to_seed_first_generation"]
        )
        num_seed_diversity = 0
    else:
        num_seed_dock_fitness = top_mols_to_seed_next_gen + diversity_depreciation

    return num_seed_diversity, num_seed_dock_fitness


def _make_pass_through_list(
    params: Dict[str, Any],
    smis_from_prev_gen: List[PreDockedCompound],
    num_elite_prev_gen: int,
    gen_num: int,
) -> Union[List[PreDockedCompound], str]:
    """
    Determines the elite ligands to advance from the previous generation without
    being altered into the next generation.

    Args:
        params (Dict[str, Any]): A dictionary of all user variables.
        smis_from_prev_gen (List[PreDockedCompound]): List of SMILES from the
            last generation chosen to seed the list of molecules to advance to
            the next generation without modification via elitism.
        num_elite_prev_gen (int): The number of molecules to advance from the
            last generation without modifications.
        gen_num (int): The integer of the current generation.

    Returns:
        Union[List[PreDockedCompound], str]: A list of ligands which should
            advance into the new generation without modifications, via elitism
            from the last generation. Returns a printout of why it failed if it
            fails.
    """
    # this will be a list of lists. Each sublist will be  [SMILES_string, ID]
    ligs_to_advance = []

    # If not enough of your previous generation sanitize to make the list
    # Return None and trigger an Error
    if len(smis_from_prev_gen) < num_elite_prev_gen:
        return (
            "Not enough ligands in initial list the filter to progress"
            + "\n len(smiles_from_previous_gen_list): {} ; \
                num_elite_prev_gen: {}".format(
                len(smis_from_prev_gen), num_elite_prev_gen,
            )
        )
    smis_from_prev_gen = [x for x in smis_from_prev_gen if type(x) == PreDockedCompound]

    ligs_that_passed_filters = [x for x in smis_from_prev_gen if x is not None]

    # If not enough of your previous generation sanitize to make the list
    # Return None and trigger an Error
    if len(ligs_that_passed_filters) < num_elite_prev_gen:
        return "Not enough ligands passed the filter to progress"
    # Save seed list of all ligands which passed which will serve as the seed
    # list.
    _save_ligand_list(
        params["output_directory"],
        gen_num,
        ligs_that_passed_filters,
        "Previous_Gen_Elite_Seed_List",
    )

    # check if ligands_which_passed_filters have docking scores
    has_dock_score = all(
        x.previous_docking_score is not None for x in ligs_that_passed_filters
    )
    # try:
    #     temp = [float(x[-2]) for x in ligands_which_passed_filters]
    #     has_dock_score = True
    # except Exception:
    #     has_dock_score = False

    # In past, generation 0 was permitted. This is no longer the case, but
    # sanity check here.
    assert gen_num > 0, "Generation number must be greater than 0"

    if not has_dock_score:
        # Take the 1st num_elite_prev_gen number of
        # molecules from ligands_which_passed_filters
        log_info(
            f"Randomly selecting {num_elite_prev_gen} compounds from the {len(ligs_that_passed_filters)} source compounds"
        )

        random.shuffle(ligs_that_passed_filters)
        ligs_to_advance = [
            ligs_that_passed_filters[x] for x in range(num_elite_prev_gen)
        ]
    else:
        # Use the make_seed_list function to select the list to advance. This
        # list will be chosen strictly by
        ligs_to_advance = _make_seed_list(
            ligs_that_passed_filters,
            gen_num,
            0,  # 0 means no diversity
            num_elite_prev_gen,  # But yes based on docking score
        )

    if len(ligs_to_advance) >= num_elite_prev_gen:
        return ligs_to_advance

    return "Not enough ligands were chosen to advance to the next generation."


#############
# Saving Output files for generations and seeds
#############
def _save_generation_smi(
    output_directory: str,
    generation_num: int,
    formatted_smile_list: List[PreDockedCompound],
    nomenclature_tag: Optional[str],
) -> Tuple[str, str]:
    """
    Save a list of newly generated compounds as an .smi file.

    This function saves the SMILES and IDs of the current generation's 
    compounds to a tab-delimited .smi file. The .smi file contains two 
    columns:
        - Column 1: SMILES string of the compound.
        - Column 2: ID of the compound.

    Args:
        output_directory (str): Directory for saving the generation files.
        generation_num (int): Current generation number.
        formatted_smile_list (List[PreDockedCompound]): List of compounds to save.
        nomenclature_tag (Optional[str]): Tag to append to the file name. If 
            None, no tag is added.

    Returns:
        Tuple[str, str]: A tuple containing:
            - The output .smi file name.
            - The path to the folder containing the generation files.
    """
    # folder for this new generation
    new_gen_folder_path = f"{output_directory}generation_{generation_num}{os.sep}"
    output_file_name = (
        f"{new_gen_folder_path}generation_{generation_num}.smi"
        if nomenclature_tag is None
        else f"{new_gen_folder_path}generation_{generation_num}{nomenclature_tag}.smi"
    )
    # write as a tab delineated .smi file
    with open(output_file_name, "w") as f:
        for smile in formatted_smile_list:
            # smile_string = smile[0]
            # smile_id = smile[1]
            x = smile.smiles + "\t" + smile.name + "\n"
            f.write(x)

    sys.stdout.flush()
    return output_file_name, new_gen_folder_path


def _save_ligand_list(
    output_directory: str,
    generation_num: int,
    chosen_ligands: List[PreDockedCompound],
    nomenclature_tag: str,
) -> None:
    """
    Save a list of ligands to an .smi file.

    This function saves a specified list of ligands (e.g., 'Mutation', 
    'Crossover', or 'Previous_Gen_choice') to an .smi file. The nomenclature_tag 
    describes the purpose of the ligand list. If the tag indicates 'seeding', 
    it represents ligands carried over from the previous generation to seed the 
    next one.

    The .smi file is tab-delimited with two columns:
        - Column 1: SMILES string of the ligand.
        - Column 2: ID of the ligand.

    Args:
        output_directory (str): Directory to save the generation files.
        generation_num (int): The current generation number.
        chosen_ligands (List[PreDockedCompound]): The list of ligands to save 
            for the next generation.
        nomenclature_tag (str): Description of the ligand list, such as 
            'Mutation', 'Crossover', 'Previous_Gen_choice', or 'seeding'.
    """
    # make a folder for the new generation
    new_gen_folder_path = f"{output_directory}generation_{generation_num}{os.sep}"

    # make a folder for the Seed files
    seed_folder_path = f"{new_gen_folder_path}SeedFolder{os.sep}"

    # check if folders exist, if not make them
    if not os.path.isdir(new_gen_folder_path):
        os.makedirs(new_gen_folder_path)
    if not os.path.isdir(seed_folder_path):
        os.makedirs(seed_folder_path)

    output_file_name = f"{seed_folder_path}{nomenclature_tag}_Gen_{generation_num}.smi"

    # save to a new output smiles file. ie. save to ranked_smiles_file
    with open(output_file_name, "w") as output:
        for chosen_ligand in chosen_ligands:
            output_line = ""
            try:
                output_line = "\t".join(chosen_ligand.to_list()) + "\n"
            except Exception:
                import pdb

                pdb.set_trace()
            output.write(output_line)

    sys.stdout.flush()
