"""
This script handles the docking and file conversion for docking.
"""
import __future__

import os
from typing import Any, Dict, List, Optional, Type, Union, cast


# from autogrow.docking.docking_class.docking_class_children \
#                           import VinaDocking, QuickVina2Docking

from autogrow.docking.docking_class.docking_file_conversion import *
from autogrow.docking.docking_class.get_child_class import get_all_subclasses
from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter
from autogrow.plugins.docking import DockingPluginManager
from autogrow.plugins.plugin_manager_base import get_plugin_manager
from autogrow.types import PostDockedCompound, PreDockedCompound


def pick_run_conversion_class_dict(
    conversion_choice: str,
) -> Type[ParentPDBQTConverter]:
    """
    This will retrieve all the names of every child class of the parent class
    ParentDocking

    Inputs:
    :param list conversion_choice: List with the User specified docking
        choices

    Returns:
    :returns: object child_dict[conversion_choice]: the class for running the
        chosen docking method
    """

    children = get_all_subclasses(ParentPDBQTConverter)

    child_dict = {child.__name__: child for child in children}
    return child_dict[conversion_choice]


def run_docking_common(
    params: Dict[str, Any],
    current_gen_int: int,
    current_generation_dir: str,
    smiles_file_new_gen: str,
    new_gen_predock_cmpds: List[PreDockedCompound],
) -> str:
    """
    This section runs the functions common to all Docking programs.

    IF ONE INCORPORATES A NEW DOCKING SOFTWARE, CONFIRM THAT ITS INPUT/OUTPUTS
    CONFORM TO THIS SECTION.

    ############## VERY IMPORTANT SECTION ########################

    Inputs:
    :param dict params: User variables which will govern how the programs runs
    :param int current_gen_int: the interger of the current generation indexed
        to zero
    :param str current_generation_dir: the current generation directory to
        find the subfolder with pdb files
    :param str smile_file_new_gen: the name of the file containing the
        molecules in the new population
    :param list new_gen_predock_cmpnds: a list of PreDockedCompound
        objects for the new generation

    Returns:
    :returns: str unweighted_ranked_smile_file: the name of the
        unweighted-ranked SMILES with their docking score
    """

    # Get directory string of PDB files for Ligands
    # current_generation_pdb_dir = f"{current_generation_dir}PDBs{os.sep}"

    # conversion_choice = params["conversion_choice"]

    # temp_vars = {
    #     key: params[key] for key in list(params.keys()) if key != "parallelizer"
    # }
    # file_conversion_cls = pick_run_conversion_class_dict(conversion_choice)

    # # TODO: new file_conversion_class_object automatically converts receptor. To
    # # convert ligand, must access object function. That's pretty awkward.
    # file_conversion_obj = file_conversion_cls(temp_vars, test_boot=False)

    docking_plugin_manager = cast(
        DockingPluginManager, get_plugin_manager("DockingPluginManager")
    )

    # dock_class = pick_docking_class_dict(dock_choice)
    # docking = dock_class(temp_vars, receptor, file_conversion_obj, test_boot=False)

    # # Find PDB's
    # pdbs_in_folder = docking_plugin_manager.find_pdb_ligands(current_generation_pdb_dir)
    # job_input_convert_lig = tuple(
    #     (docking_plugin_manager, pdb, file_conversion_obj) for pdb in pdbs_in_folder
    # )

    # print("Convert Ligand to PDBQT format Begun")
    # smiles_names_failed_to_convert = params["parallelizer"].run(
    #     job_input_convert_lig, lig_convert_multithread
    # )

    # print("Convert Ligand to PDBQT format Completed")
    # deleted_smiles_names_list_convert = [
    #     x for x in smiles_names_failed_to_convert if x is not None
    # ]
    # deleted_smiles_names_list_convert = list(set(deleted_smiles_names_list_convert))

    # if deleted_smiles_names_list_convert:
    #     print("THE FOLLOWING LIGANDS FAILED TO CONVERT:")
    #     print(deleted_smiles_names_list_convert)

    # # Docking the ligands which converted to PDBQT Find PDBQT's
    # pdbqts_in_folder = docking_plugin_manager.find_converted_ligands(
    #     current_generation_pdb_dir
    # )

    job_input_dock_lig = [
        (docking_plugin_manager, new_gen_predock_cmpd)
        for new_gen_predock_cmpd in new_gen_predock_cmpds
    ]
    print("Docking Begun")
    post_docked_compounds: List[PostDockedCompound] = params["parallelizer"].run(
        job_input_dock_lig, _run_dock_multithread_wrapper
    )

    # Remove those that failed to convert
    post_docked_compounds = [x for x in post_docked_compounds if x is not None]

    # Remove those not associated with a docked sdf file
    post_docked_compounds = [
        x for x in post_docked_compounds if x.docked_sdf_path is not None
    ]

    print("\nBegin Ranking and Saving results")
    unweighted_ranked_smile_file = docking_plugin_manager.rank_and_save_output_smi(
        current_generation_dir, current_gen_int, smiles_file_new_gen,
        post_docked_compounds
    )
    print("\nCompleted Ranking and Saving results\n")
    return unweighted_ranked_smile_file


def lig_convert_multithread(
    docking: DockingPluginManager,
    pdb: str,
    file_conversion_class_object: ParentPDBQTConverter,
) -> Union[str, None]:
    """
    Run the ligand conversion of a single molecule. If it failed
    failed_smiles_name will be a string of the SMILE which failed to convert
    If it converts failed_smiles_name will be a None.

    Inputs:
    :param object docking_object: the class for running the chosen docking
        method
    :param str pdb: the path to the pdb of a molecule

    Returns:
    :returns: list failed_smiles_name: if the molecule failed to convert to
        final format. (ie. pdbqt conversion fail)
    """

    return docking.run_ligand_handling_for_docking(pdb, file_conversion_class_object)


def _run_dock_multithread_wrapper(
    docking: DockingPluginManager, new_gen_predock_cmpd: PreDockedCompound,
) -> Optional[float]:
    """
    Run the docking of a single molecule.

    Inputs:
    :param object docking_object: the class for running the chosen docking
        method
    :param str pdb: the path to the pdb of a molecule

    Returns:
    :returns: list failed_smiles_names: any smiles which were deleted (ie.
        docking failed)
    """

    # print("Attempt to Dock complete: ", lig_pdb)
    return docking.run(predocked_cmpd=new_gen_predock_cmpd)
