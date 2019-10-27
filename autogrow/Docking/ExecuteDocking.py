import __future__

import glob
import os
import copy

import autogrow.Docking.Delete_failed_mol as Delete
from autogrow.Docking.Docking_Class.get_child_class import get_all_subclasses

from autogrow.Docking.Docking_Class.DockingClassChildren import *
from autogrow.Docking.Docking_Class.ParentDockClass import ParentDocking
# from autogrow.Docking.Docking_Class.DockingClassChildren import VinaDocking, QuickVina2Docking

from autogrow.Docking.Docking_Class.Docking_File_Conversion import *
from autogrow.Docking.Docking_Class.ParentPDBQTConverter import ParentPDBQTConverter
# from autogrow.Docking.Docking_Class.Docking_File_Conversion import Convert_with_MGLTOOLS, Convert_with_obabel


def pick_docking_class_dict(Dock_choice):
    """
    This will retrieve all the names of every child class of the parent class ParentDocking

    
    Input:
    :param list Dock_choice: List with the User specified docking choices
    
    Return:
    :returns: object child_dict[Dock_choice]: the class for running the chosen docking method
    """ 
    children = get_all_subclasses(ParentDocking)

    child_dict = {}
    for child in children:
        childName = child.__name__
        child_dict[childName] = child

    return child_dict[Dock_choice]
#

def pick_run_conversion_class_dict(Conversion_choice):
    """
    This will retrieve all the names of every child class of the parent class ParentDocking

    
    Input:
    :param list Conversion_choice: List with the User specified docking choices
    
    Return:
    :returns: object child_dict[Conversion_choice]: the class for running the chosen docking method
    """ 
    children = get_all_subclasses(ParentPDBQTConverter)

    child_dict = {}
    for child in children:
        childName = child.__name__
        child_dict[childName] = child
    
    return child_dict[Conversion_choice]
#

def run_docking_common(vars, current_gen_int, current_generation_dir, smile_file_new_gen):
    """
    This section runs the functions common to all Docking programs.

    IF ONE INCORPORATES A NEW DOCKING SOFTWARE, CONFIRM THAT ITS INPUT/OUTPUTS CONFORM TO THIS SECTION.
    ############## VERY IMPORTANT SECTION######################## 

    Input:
    :param dict vars: User variables which will govern how the programs runs
    :param int current_gen_int: the interger of the current generation indexed to zero
    :param str current_generation_dir: the current generation directory to find the subfolder with pdb files
    :param str smile_file_new_gen: the name of the file containing the molecules in the new population
    Return:
    :returns: str unweighted_ranked_smile_file: the name of the unweighted-ranked SMILES with their docking score
    """
    # Get directory string of PDB files for Ligands
    current_generation_PDB_dir = current_generation_dir + "PDBs" + os.sep

    Dock_choice = vars["Dock_choice"]
    Conversion_choice = vars["Conversion_choice"]
    receptor = vars["filename_of_receptor"]

    # Use a temp vars dict so you don't put mpi multiprocess info through itself...
    temp_vars = {}
    for key in list(vars.keys()):
        if key =="Parallelizer":
            continue
        temp_vars[key]= vars[key]

    file_conversion_class_object = pick_run_conversion_class_dict(Conversion_choice)
    file_conversion_class_object = file_conversion_class_object(temp_vars, receptor, test_boot=False)


    
    dock_class = pick_docking_class_dict(Dock_choice)
    dockingObject = dock_class(temp_vars, receptor, file_conversion_class_object,test_boot=False)

    if vars["docking_executable"] == None:
        docking_executable = dockingObject.get_docking_executable_file(vars)
        vars["docking_executable"] = docking_executable

    # Find PDB's
    pdbs_in_folder = dockingObject.find_pdb_ligands(current_generation_PDB_dir)
    job_input_convert_lig = tuple([tuple([dockingObject, pdb]) for pdb in pdbs_in_folder])

    print("####################")
    print("Convert Ligand to PDBQT format Begun")
    smiles_names_failed_to_convert = vars['Parallelizer'].run(job_input_convert_lig, lig_convert_multithread)
    
    print("Convert Ligand to PDBQT format Completed")
    deleted_smiles_names_list_convert = [x for x in smiles_names_failed_to_convert if x is not None]
    deleted_smiles_names_list_convert = list(set(deleted_smiles_names_list_convert))

    if len(deleted_smiles_names_list_convert) != 0:
        print("THE FOLLOWING LIGANDS WHICH FAILED TO CONVERT:")
        print(deleted_smiles_names_list_convert)
    print("####################")
    
    # Docking the ligands which converted to PDBQT
    # Find PDBQT's
    pdbqts_in_folder = dockingObject.find_converted_ligands(current_generation_PDB_dir)

    job_input_dock_lig = tuple([tuple([dockingObject, pdbqt]) for pdbqt in pdbqts_in_folder])
    print("####################")
    print("Docking Begun")
    smiles_names_failed_to_dock = vars['Parallelizer'].run(job_input_dock_lig, run_dock_multithread)

    print("")
    print("")
    print("")
    print("Docking Completed")    
    print("####################")

    deleted_smiles_names_list_dock = [x for x in smiles_names_failed_to_dock if x is not None]
    deleted_smiles_names_list_dock = list(set(deleted_smiles_names_list_dock))

    if len(deleted_smiles_names_list_dock) != 0:
        print("THE FOLLOWING LIGANDS WHICH FAILED TO DOCK:")
        print(deleted_smiles_names_list_dock)
        
    
    print("####################")
    deleted_smiles_names_list = deleted_smiles_names_list_convert + deleted_smiles_names_list_dock

    if len(deleted_smiles_names_list) != 0:
        print("")
        print("THE FOLLOWING LIGANDS WHERE DELETED FOR FAILURE TO CONVERT OR DOCK:")
        print(deleted_smiles_names_list)
        

    print("#################### ")
    print("")
    print("Begin Ranking and Saving results")
    unweighted_ranked_smile_file = dockingObject.rank_and_save_output_smi(vars, current_generation_dir, current_gen_int, smile_file_new_gen, deleted_smiles_names_list)
    print("")
    print("Completed Ranking and Saving results")
    print("")

    return unweighted_ranked_smile_file
#

def lig_convert_multithread(dockingObject, pdb):
    """
    Run the ligand conversion of a single molecule.
    If it failed failed_smiles_name will be a string of the SMILE which failed to convert
    If it converts failed_smiles_name will be a None. 

    Input:
    :param object dockingObject: the class for running the chosen docking method
    :param str pdb: the path to the pdb of a molecule

    Return:
    :returns: list failed_smiles_name: if the molecule failed to convert to final format. 
                            (ie. pdbqt conversion fail)
    """
    failed_smiles_name = dockingObject.run_ligand_handling_for_docking(pdb)
    return failed_smiles_name
#

def run_dock_multithread(dockingObject, pdb):
    """
    Run the docking of a single molecule.

    Input:
    :param object dockingObject: the class for running the chosen docking method
    :param str pdb: the path to the pdb of a molecule

    Return:
    :returns: list failed_smiles_names: any smiles which were deleted 
                            (ie. docking failed)
    """
    print("Attempt to Dock complete: ", pdb)
    failed_smiles_names = dockingObject.run_dock(pdb)
    return failed_smiles_names
#
