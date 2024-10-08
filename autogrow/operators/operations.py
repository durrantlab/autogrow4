"""
Populates an AutoGrow generation via mutation, crossover, and elitism.
Also filters and converts SMILES to 3d SDFS.
"""
import __future__

import os
import random
import copy
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

from autogrow.plugins.plugin_manager_base import get_plugin_manager
from autogrow.types import PreDockedCompoundInfo
from autogrow.utils.logging import log_info
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.docking.ranking.ranking_mol as Ranking
import autogrow.operators.mutation.execute_mutations as Mutation
import autogrow.operators.crossover.execute_crossover as execute_crossover
import autogrow.operators.convert_files.conversion_to_3d as conversion_to_3d
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


#############
# Main run Autogrow operators to make a generation
#############
def populate_generation(
    params: Dict[str, Any], generation_num: int
) -> Tuple[str, List[List[str]]]:
    """
    This will run all of the mutations, crossovers, and filters for a single
        generation. Populates a new generation of ligands.

    Inputs:
    :param dict params: a dictionary of all user variables
    :param int generation_num: the generation number

    Returns:
    :returns: str full_generation_smiles_file: the name of the .smi file
        containing the new population
    :returns: list full_gen_smis: list with the new population
        of ligands
    :returns: bool None: returns None twice if any step failed. This will
        result in the program ending
    """
    number_of_processors = int(params["number_of_processors"])

    # Determine which generation it is and how many mutations and crossovers
    # to make
    if generation_num == 1:
        # If 1st generation
        num_crossovers = params["number_of_crossovers_first_generation"]
        num_mutations = params["number_of_mutants_first_generation"]

        # How many advance from previous generation to the next generation
        # directly This will be done later but we are unpacking params here
        num_elite_to_advance_from_previous_gen = params[
            "number_elitism_advance_from_previous_gen_first_generation"
        ]
    else:
        # Later generations
        num_crossovers = params["number_of_crossovers"]
        num_mutations = params["number_of_mutants"]
        num_elite_to_advance_from_previous_gen = params[
            "number_elitism_advance_from_previous_gen"
        ]

    # Get the Source compound list. This list is the full population from
    # either the previous generations or if its Generation 1 than the its the
    # entire User specified Source compound list If either has a SMILES that
    # does not sanitize in RDKit it will be excluded and a printout of its
    # Name and SMILES string will be printed.
    src_cmpds = get_complete_list_prev_gen_or_source_compounds(params, generation_num)

    num_seed_diversity, num_seed_dock_fitness = determine_seed_population_sizes(
        params, generation_num
    )

    # Total Population size of this generation
    total_num_desired_new_ligands = (
        num_crossovers + num_mutations + num_elite_to_advance_from_previous_gen
    )

    # Get starting compounds for Mutations
    seed_list_mutations = make_seed_list(
        params, src_cmpds, generation_num, num_seed_diversity, num_seed_dock_fitness,
    )

    # Save seed list for Mutations
    save_ligand_list(
        params["output_directory"],
        generation_num,
        seed_list_mutations,
        "Mutation_Seed_List",
    )
    sys.stdout.flush()

    print("MAKE MUTATIONS")
    # Making Mutations

    # Package user params specifying the Reaction library to use for mutation
    rxn_library_variables = [
        params["rxn_library"],
        params["rxn_library_file"],
        params["function_group_library"],
        params["complementary_mol_directory"],
    ]

    # List of SMILES from mutation
    new_mutation_smiles_list: List[PreDockedCompoundInfo] = []

    # Make all the required ligands by mutations
    while len(new_mutation_smiles_list) < num_mutations:
        sys.stdout.flush()

        num_mutants_to_make = num_mutations - len(new_mutation_smiles_list)

        # Make all mutants
        new_mutants = Mutation.make_mutants(
            params,
            generation_num,
            number_of_processors,
            num_mutants_to_make,
            seed_list_mutations,
            new_mutation_smiles_list,
            rxn_library_variables,
        )
        if new_mutants is None:
            # try once more
            new_mutants = Mutation.make_mutants(
                params,
                generation_num,
                number_of_processors,
                num_mutants_to_make,
                seed_list_mutations,
                new_mutation_smiles_list,
                rxn_library_variables,
            )

        if new_mutants is None:
            break

        # Remove Nones:
        new_mutants = [x for x in new_mutants if x is not None]

        for i in new_mutants:
            new_mutation_smiles_list.append(i)
            if len(new_mutation_smiles_list) == num_mutations:
                break
    sys.stdout.flush()

    # save new_mutation_smiles_list
    save_ligand_list(
        params["output_directory"],
        generation_num,
        new_mutation_smiles_list,
        "Chosen_Mutants",
    )

    if (
        new_mutation_smiles_list is None
        or len(new_mutation_smiles_list) < num_mutations
    ):
        _extracted_from_populate_generation_142(
            num_mutations,
            " ligands through Mutation",
            new_mutation_smiles_list,
            "Mutation failed to make enough new ligands.",
        )
    print("FINISHED MAKING MUTATIONS")

    # Get starting compounds to seed Crossovers
    seed_list_crossovers = make_seed_list(
        params, src_cmpds, generation_num, num_seed_diversity, num_seed_dock_fitness,
    )

    # Save seed list for Crossovers
    save_ligand_list(
        params["output_directory"],
        generation_num,
        seed_list_crossovers,
        "Crossover_Seed_List",
    )

    print("MAKE CROSSOVERS")
    sys.stdout.flush()

    # Making Crossovers
    # List of smiles from crossover
    new_crossover_smiles_list: List[PreDockedCompoundInfo] = []

    # Make all the required ligands by Crossover
    while len(new_crossover_smiles_list) < num_crossovers:
        sys.stdout.flush()
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
            # try once more
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

        # Remove Nones:
        new_crossovers = [x for x in new_crossovers if x is not None]

        # append those which passed the filter
        for i in new_crossovers:
            new_crossover_smiles_list.append(i)
            if len(new_crossover_smiles_list) == num_crossovers:
                break

    # save new_crossover_smiles_list
    save_ligand_list(
        params["output_directory"],
        generation_num,
        new_crossover_smiles_list,
        "Chosen_Crossovers",
    )

    if (
        new_crossover_smiles_list is None
        or len(new_crossover_smiles_list) < num_crossovers
    ):
        _extracted_from_populate_generation_142(
            num_crossovers,
            " ligands through Crossover",
            new_crossover_smiles_list,
            "Crossover failed to make enough new ligands.",
        )
    print("FINISHED MAKING CROSSOVERS")

    # Get unaltered samples from the previous generation
    print("GET SOME LIGANDS FROM THE LAST GENERATION")
    sys.stdout.flush()

    # Make a list of the ligands chosen to pass through to the next generation
    # via Elitism This handles creating a seed list and defining the advance
    # to next generation final selection

    chosen_mol_to_pass_through_list = make_pass_through_list(
        params, src_cmpds, num_elite_to_advance_from_previous_gen, generation_num,
    )

    if type(chosen_mol_to_pass_through_list) == str:
        printout = (
            chosen_mol_to_pass_through_list
            + "\nIf this is the 1st generation, it may be due to the starting "
            + "library has SMILES which could not be converted to sanitizable "
            + "RDKit Molecules"
        )

        raise Exception(printout)
    sys.stdout.flush()

    assert (
        type(chosen_mol_to_pass_through_list) is list
    ), "Chosen mol to pass through list is not a list"

    # save chosen_mol_to_pass_through_list
    save_ligand_list(
        params["output_directory"],
        generation_num,
        chosen_mol_to_pass_through_list,
        "Chosen_Elite_To_advance",
    )

    print("GOT LIGANDS FROM THE LAST GENERATION")

    # make a list of all the ligands from mutations, crossovers, and from the
    # last generation
    new_gen_smis = []
    full_gen_smis = []
    for i in new_mutation_smiles_list:
        new_gen_smis.append(i)
        full_gen_smis.append(i)

    for i in new_crossover_smiles_list:
        new_gen_smis.append(i)
        full_gen_smis.append(i)

    if params["redock_elite_from_previous_gen"] is False and generation_num != 1:
        for i in chosen_mol_to_pass_through_list:
            # Doesn't append to the new_generation_smiles_list
            full_gen_smis.append(i)

    # Generation 0 pass through gets added to the convert and dock list
    # because it has no docking score to compare with This is independent of
    # the params['redock_elite_from_previous_gen']
    else:
        for i in chosen_mol_to_pass_through_list:
            new_gen_smis.append(i)
            full_gen_smis.append(i)

    if len(full_gen_smis) < total_num_desired_new_ligands:
        print("We needed ", total_num_desired_new_ligands)
        print("We made ", len(full_gen_smis))
        print(
            "population failed to make enough mutants or crossovers... \
            Errors could include not enough diversity, too few seeds to \
            the generation, the seed mols are unable to cross-over due \
            to lack of similariy, or all of the seed lack functional groups \
            for performing reactions"
        )
        assert False

    # Save the Full Generation
    full_generation_smiles_file, new_gen_folder_path = save_generation_smi(
        params["output_directory"], generation_num, full_gen_smis, None
    )

    # Save the File to convert to 3d
    smiles_to_convert_file, new_gen_folder_path = save_generation_smi(
        params["output_directory"], generation_num, new_gen_smis, "_to_convert",
    )

    sys.stdout.flush()
    # CONVERT SMILES TO .sdf USING GYPSUM and convert .sdf to .pdb with rdkit
    # This will output sdf files into a folder. The .smi.0.sdf file is not a
    # valid mol, but all the others will be valid the 1st Smiles in the
    # original .smi file is saved as .smi.1.sdf and 2nd file is saved as
    # .smi.2.sdf
    conversion_to_3d.convert_to_3d(params, smiles_to_convert_file, new_gen_folder_path)
    sys.stdout.flush()

    return full_generation_smiles_file, full_gen_smis


# TODO: Rename this here and in `populate_generation`
def _extracted_from_populate_generation_142(arg0, arg1, arg2, arg3):
    print("")
    print("")
    print(f"We needed to make {arg0}{arg1}")
    print(f"We only made {len(arg2)}{arg1}")
    print("")
    print("")
    raise Exception(arg3)


# TODO: There's a lot of overlap between populate_generation_zero() and
# populate_generation()
def populate_generation_zero(
    params: Dict[str, Any], gen_num: int = 0
) -> Tuple[bool, str, List[PreDockedCompoundInfo]]:
    """
    This will handle all that is required for generation 0redock and handle
    the generation 0

    Inputs:
    :param dict params: a dictionary of all user variables
    :param int generation_num: the generation number

    Returns:
    :returns: str full_generation_smiles_file: the name of the .smi file
        containing the new population
    :returns: list full_gen_smis: list with the new population
        of ligands.
    :returns: bool already_docked: if true we won't redock the source ligands.
        If False we will dock the source ligands.
    """
    number_of_processors = int(params["number_of_processors"])

    num_crossovers = 0
    num_mutations = 0

    # Get the Source compound list. This list is the full population from
    # either the previous generations or if its Generation 1 than the its the
    # entire User specified Source compound list If either has a SMILES that
    # does not sanitize in RDKit it will be excluded and a printout of its
    # Name and SMILES string will be printed.
    src_cmpds = get_complete_list_prev_gen_or_source_compounds(params, gen_num)
    num_elite_to_advance_from_previous_gen = len(src_cmpds)

    num_seed_diversity, num_seed_dock_fitness = determine_seed_population_sizes(
        params, gen_num
    )

    # Total Population size of this generation
    total_num_desired_new_ligands = num_crossovers + num_mutations + 1

    # Get unaltered samples from the previous generation
    print("GET SOME LIGANDS FROM THE LAST GENERATION")

    # Make a list of the ligands chosen to pass through to the next generation
    # This handles creating a seed list and defining the advance to next
    # generation final selection
    mols_to_pass_through = make_pass_through_list(params, src_cmpds, 1, 0)

    if type(mols_to_pass_through) == str:
        printout = (
            mols_to_pass_through
            + "\nIf this is the 1st generation, it may be due to the starting "
            + "library has SMILES which could not be converted to "
            + "sanitizable RDKit Molecules"
        )

        raise Exception(printout)

    assert (
        type(mols_to_pass_through) is list
    ), "Chosen mol to pass through list is not a list"

    # save chosen_mol_to_pass_through_list
    save_ligand_list(
        params["output_directory"],
        gen_num,
        mols_to_pass_through,
        "Chosen_Elite_To_advance",
    )

    print("GOT LIGANDS FROM THE LAST GENERATION")

    # make a list of all the ligands from mutations, crossovers, and from the
    # last generation
    new_gen_smis: List[PreDockedCompoundInfo] = []
    full_gen_smis: List[PreDockedCompoundInfo] = []

    # These will be docked and scored for generation 0
    for i in mols_to_pass_through:
        new_gen_smis.append(i)
        full_gen_smis.append(i)

    assert (
        full_gen_smis
    ), "population failed to import any molecules from the source_compounds_list."

    # Save the Full Generation
    full_gen_smis_file, new_gen_folder_path = save_generation_smi(
        params["output_directory"], gen_num, full_gen_smis, None
    )

    # Save the File to convert to 3d
    smis_to_convert_file, new_gen_folder_path = save_generation_smi(
        params["output_directory"], gen_num, new_gen_smis, "_to_convert",
    )

    # order files by -2 (docking score) of each lig
    full_gen_smis_printout = []
    try:

        def _get_score(x: PreDockedCompoundInfo) -> float:
            assert x.previous_docking_score is not None, "Docking score is None"
            return float(x.previous_docking_score)

        full_gen_smis.sort(key=lambda x: _get_score(x), reverse=False)
        full_gen_smis_printout = ["\t".join(x.smiles) for x in full_gen_smis]
        already_docked = True
    except Exception:
        print(
            "Not all ligands in source compound list are scored. "
            + "We will convert and redock them all."
        )
        already_docked = False

    if already_docked:
        # Write all the ligands to a ranked file
        full_gen_smis_printout = "\n".join(full_gen_smis_printout)
        ranked_file = f"{new_gen_folder_path}{os.sep}generation_0_ranked.smi"
        with open(ranked_file, "w") as f:
            f.write(full_gen_smis_printout)
        return already_docked, full_gen_smis_file, full_gen_smis

    # If you are to redock and convert the generation zero you will also need
    # to do the following:

    # CONVERT SMILES TO .sdf USING GYPSUM and convert .sdf to .pdb with
    # rdkit This will output sdf files into a folder. The .smi.0.sdf file
    # is not a valid mol, but all the others will be valid the 1st Smiles
    # in the original .smi file is saved as .smi.1.sdf and 2nd file is
    # saved as .smi.2.sdf
    conversion_to_3d.convert_to_3d(params, smis_to_convert_file, new_gen_folder_path)

    return already_docked, full_gen_smis_file, full_gen_smis


#############
# Get seeds
#############
def test_source_smiles_convert(
    smile_info: PreDockedCompoundInfo,
) -> Union[PreDockedCompoundInfo, str]:
    """
    This attempts to convert a SMILES string to an rdkit.Chem.rdchem.Mol
    object
        -done in a try statement so that it is some bad SMILES string which is
        incapable of being converted
        - it also checks that the SMILES string is able to be sanitized

    Inputs:
    :param list smile_info: a list containing the SMILES of a ligand, its ID
        and potentially additional information about the ligand

    Returns:
    :returns: list smile_info: If it passed the test, it returns the list
        containing the SMILES of a ligand, its ID and potentially additional
        information about the ligand
    :returns: str printout: If it failed to convert it returns the error
        message. This passess out to prevent MPI print issues
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
        return _extracted_from_test_source_smiles_convert_90(
            smile_str, printout, smile_id
        )
    # Check there are no * which are atoms with atomic number=0
    mol = MOH.check_for_unassigned_atom(mol)
    if mol is None:
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: SMILES string contained "
            + "an unassigned atom type labeled as *.\n"
        )
        return _extracted_from_test_source_smiles_convert_90(
            smile_str, printout, smile_id
        )
    # Check for fragments.
    if len(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)) != 1:

        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string was fragmented.\n"
        return _extracted_from_test_source_smiles_convert_90(
            smile_str, printout, smile_id
        )
    # the ligand is good enough to use throughout the program!
    return smile_info


# TODO Rename this here and in `test_source_smiles_convert`
def _extracted_from_test_source_smiles_convert_90(smile_str, printout, smile_id):
    printout += f"\t Removed SMILE string is: {smile_str} \n"
    printout += f"\t Removed SMILE ID is: {smile_id}"
    return printout


def get_complete_list_prev_gen_or_source_compounds(
    params: Dict[str, Any], generation_num: int
) -> List[PreDockedCompoundInfo]:
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
    :returns: list usable_smiles: a list with SMILES strings, names,
        and information about the smiles from the previous generation or the
        source compound list
    """
    source_file_gen_0 = (
        f"{params['output_directory']}generation_0{os.sep}generation_0_ranked.smi"
    )
    if generation_num == 0:
        usable_smiles = _get_source_compounds_or_raise(params)
    elif generation_num == 1 and os.path.exists(source_file_gen_0) is False:
        usable_smiles = _get_source_compounds_or_raise(params)
    else:
        source_file = (
            params["output_directory"]
            + f"generation_{generation_num - 1}{os.sep}generation_{generation_num - 1}_ranked.smi"
        )
        if os.path.exists(source_file) is False:
            _handle_no_ligands_found("\tCheck formatting or if file has been moved.\n")
        usable_smiles = Ranking.get_usable_format(source_file)

        if len(usable_smiles) == 0:
            _handle_no_ligands_found("\tCheck formatting or if file has been moved. \n")
    # Test that every SMILES in the usable_smiles is a valid SMILES
    # which will import and Sanitize in RDKit. SMILES will be excluded if they
    # are fragmented, contain atoms with no atomic number (*), or do not
    # sanitize
    job_input = tuple((i,) for i in usable_smiles)

    usable_smiles: List[PreDockedCompoundInfo] = params["parallelizer"].run(
        job_input, test_source_smiles_convert
    )
    usable_smiles = [x for x in usable_smiles if x is not None]
    print_errors = [x for x in usable_smiles if type(x) is str]
    usable_smiles = [x for x in usable_smiles if type(x) is PreDockedCompoundInfo]
    for x in print_errors:
        print(x)

    if not usable_smiles:
        _raise_exception_with_message(
            "\nThere were no ligands in source compound or previous \
            generation which could sanitize.\n"
        )

    random.shuffle(usable_smiles)

    return usable_smiles


def _raise_exception_with_message(arg0):
    printout = arg0
    print(printout)
    raise Exception(printout)


def _get_source_compounds_or_raise(params) -> List[PreDockedCompoundInfo]:
    # This will be the full length list of starting molecules as the seed
    source_file = str(params["source_compound_file"])
    result = Ranking.get_usable_format(source_file)

    if len(result) == 0:
        print(
            "\nThere were no available ligands in source compound. Check formatting\n"
        )
        raise Exception(
            "There were no available ligands in source compound. Check formatting"
        )

    return result


def _handle_no_ligands_found(arg0):
    printout = (
        "\n"
        + "There were no available ligands in previous"
        + " generation ranked ligand file.\n"
    )
    printout = printout + arg0
    print(printout)
    raise Exception(printout)


def make_seed_list(
    params: Dict[str, Any],
    source_compounds_list: List[PreDockedCompoundInfo],
    generation_num: int,
    num_seed_diversity: int,
    num_seed_dock_fitness: int,
) -> List[PreDockedCompoundInfo]:
    """
    Get the starting compound list for running the Mutation and Crossovers

    If generation = 0 use the User specified starting compounds If generation
    is >0 than use the previous generations top ligands. This takes an .smi
    file

    Inputs:
    :param dict params: a dictionary of all user variables
    :param list source_compounds_list: a list with SMILES strings, names, and
        information about the smiles from either the previous generation or the
        source compound list
    :param int generation_num: the interger of the current generation
    :param int num_seed_diversity: the number of seed molecules which come
        from diversity selection
    :param int num_seed_dock_fitness: the number of seed molecules which come
        from eite selection by docking score

    Returns:
    :returns: list usable_smiles: a list with SMILES strings, names,
        and information about the smiles which will be used to seed the next
        generation
    """
    usable_smiles = copy.deepcopy(source_compounds_list)

    full_length = False
    if generation_num == 0:
        # Get starting compounds for Mutations
        full_length = True
    elif generation_num == 1:
        if params["dock_source_compounds_first"] is False:
            # Get starting compounds for Mutations
            full_length = True
        else:
            source_file_gen_0 = params[
                "output_directory"
            ] + "generation_{}{}generation_{}_ranked.smi".format(0, os.sep, 0)
            if not os.path.exists(source_file_gen_0):
                full_length = True

            else:
                # generation_num 1 may run into problems if the source
                # compounds are smaller than the seed pool required to seed
                # generation 1. Because the seeding options are connected to
                # the generation number (due to the depreciation of diversity
                # option) Because of this we may need to ignore the ranking
                # for the seeds of generation 1 to accomidate the smaller
                # source size. This is especially important with
                # lead-optimization in which the source pool may be much
                # smaller For this reason we will override the seeding of
                # generation 1 if the number to seed is greater than exists
                # but will provide a warning message.
                if (
                    len(usable_smiles) < num_seed_diversity
                    or len(usable_smiles) < num_seed_diversity
                ):
                    # This is problematic so just use what is available
                    printout = "\n\nNot enough ligands in source compound \
                        list to seed generation 1. We will use the entire \
                        list of every ligand in the source compound list \
                        to seed generation 1. This means there is no \
                        selection in generation 1's seeding process.\n\n"
                    print(printout)
                    full_length = True
                else:
                    full_length = False
    else:
        full_length = False

    if full_length is True or generation_num == 0:
        # This will be the full length list of starting molecules as the seed
        random.shuffle(usable_smiles)

    else:
        # selector_choice = params["selector_choice"]
        # tourn_size = params["tourn_size"]
        # Get subset of the source_file based on diversity scores and docking
        # scores
        usable_smiles = Ranking.create_seed_list(
            usable_smiles, num_seed_diversity, num_seed_dock_fitness,
        )

    random.shuffle(usable_smiles)

    return usable_smiles


def determine_seed_population_sizes(
    params: Dict[str, Any], gen_num: int
) -> Tuple[int, int]:
    """
    This function determines how many molecules will be chosen to seed a
    generation because of their docking score and how many will assess because
    of their diversity score.

    Inputs:
    :param dict params: a dictionary of all user variables
    :param int generation_num: the interger of the current generation

    Returns:
    :returns:  int num_seed_diversity: the number of seed molecules which come
        from diversity selection
    :returns:  int num_seed_dock_fitness: the number of seed molecules which
        come from eite selection by docking score

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


def make_pass_through_list(
    params: Dict[str, Any],
    smis_from_prev_gen: List[PreDockedCompoundInfo],
    num_elite_from_prev_gen: int,
    gen_num: int,
) -> Union[List[PreDockedCompoundInfo], str]:
    """
    This function determines the molecules which elite ligands will advance
    from the previous generation without being altered into the next
    generation.

    Inputs:
    :param dict params: a dictionary of all user variables
    :param list smis_from_prev_gen: List of SMILES from the last
        generation chosen to seed the list of molecules to advance to the next
        generation without modification via elitism.
    :param int num_elite_from_prev_gen: the number of molecules
        to advance from the last generation without modifications.
    :param int gen_num: the interger of the current generation

    Returns:
    :returns: list ligs_to_advance: a list of ligands which should
        advance into the new generation without modifications, via elitism from
        the last generation. Returns a printout of why it failed if it fails
    """
    # this will be a list of lists. Each sublist will be  [SMILES_string, ID]
    ligs_to_advance = []

    # If not enough of your previous generation sanitize to make the list
    # Return None and trigger an Error
    if gen_num != 0 and len(smis_from_prev_gen) < num_elite_from_prev_gen:
        return (
            "Not enough ligands in initial list the filter to progress"
            + "\n len(smiles_from_previous_gen_list): {} ; \
                num_elite_to_advance_from_previous_gen: {}".format(
                len(smis_from_prev_gen), num_elite_from_prev_gen,
            )
        )
    smis_from_prev_gen = [
        x for x in smis_from_prev_gen if type(x) == PreDockedCompoundInfo
    ]

    ligs_that_passed_filters = [x for x in smis_from_prev_gen if x is not None]
    
    # If not enough of your previous generation sanitize to make the list
    # Return None and trigger an Error
    if gen_num != 0 and len(ligs_that_passed_filters) < num_elite_from_prev_gen:
        return "Not enough ligands passed the filter to progress"
    # Save seed list of all ligands which passed which will serve as the seed
    # list.
    save_ligand_list(
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

    if gen_num == 0 and not has_dock_score:
        # Take the 1st num_elite_to_advance_from_previous_gen number of
        # molecules from ligands_which_passed_filters
        random.shuffle(ligs_that_passed_filters)
        ligs_to_advance = [
            ligs_that_passed_filters[x] for x in range(len(ligs_that_passed_filters))
        ]
    elif gen_num == 0:
        # Use the make_seed_list function to select the list to advance.
        # This list will be chosen strictly by
        ligs_to_advance = make_seed_list(
            params,
            ligs_that_passed_filters,
            gen_num,
            len(ligs_that_passed_filters),
            num_elite_from_prev_gen,
        )

    elif not has_dock_score:
        # Take the 1st num_elite_to_advance_from_previous_gen number of
        # molecules from ligands_which_passed_filters
        random.shuffle(ligs_that_passed_filters)
        ligs_to_advance = [
            ligs_that_passed_filters[x] for x in range(num_elite_from_prev_gen)
        ]
    else:
        # Use the make_seed_list function to select the list to advance. This
        # list will be chosen strictly by
        ligs_to_advance = make_seed_list(
            params, ligs_that_passed_filters, gen_num, 0, num_elite_from_prev_gen,
        )

    if gen_num == 0 or len(ligs_to_advance) >= num_elite_from_prev_gen:
        return ligs_to_advance
    return "Not enough ligands were chosen to advance to the next generation."


#############
# Saving Output files for generations and seeds
#############
def save_generation_smi(
    output_directory: str,
    generation_num: int,
    formatted_smile_list: List[PreDockedCompoundInfo],
    nomenclature_tag: Optional[str],
) -> Tuple[str, str]:
    """"
    This function saves a list of newly generated population of ligands as an
    .smi file. .smi file column 1 is the SMILES string and column 2 is its
    smile ID


    Inputs:
    :param dict output_directory: the directory of the run to save the
        generation
    :param int generation_num: the interger of the current generation
    :param list formatted_smile_list: list of the newly generated population
        of ligands
    :param str nomenclature_tag: The str describing the ligand list if None
        than don't add tag. It is the full list of all ligands for the
        generation
        If it says '_to_convert' its the list of ligands which will need to be
        converted to 3D -this may or may not have the ligands which pass
        through from the last gen.

    Returns:
    :returns: str output_file_name: name of the output file
    :returns: str new_gen_folder_path: the path to the folder containing all
        that will be in this generation
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


def save_ligand_list(
    output_directory: str,
    generation_num: int,
    chosen_ligands: List[PreDockedCompoundInfo],
    nomenclature_tag: str,
) -> None:
    """
    Save the list of ligands. nomenclature_tag is a string such as "Mutation"
    or "Crossover" or "Previous_Gen_choice" describing what this data is used
    for. If it says seeding it is the chosen mols from the previous generation
    being used to seed the next generation

    Inputs:
    :param dict output_directory: the directory of the run to save the
        generation
    :param int generation_num: The generation number
    :param list list_of_chosen_ligands: The formatted list of ligands to seed
        a generation
    :param str nomenclature_tag: The str describing the ligand list
        -ie seeding_mutations is the list that seeded the mutations while
            mutations would be the list of mutations generated from the
            seeding_mutations list
        -ie. mutation, crossover, Previous_Gen_choice
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
