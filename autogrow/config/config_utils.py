
from typing import List
import os
from shutil import copyfile

def make_complete_children_dict(purpose_of_object):
    """
    This will retrieve all the names of every child class of the parent class
    This can be either filter, parent_pdbqt_converter, ParentDocking,
    or ParentScoring

    Inputs:
    :param str purpose_of_object: either filter, parent_pdbqt_converter,
        ParentDocking, or ParentScoring
    Returns:
    :returns: dict child_dict: Dictionary of all the class objects for either
        Filtering, docking, Dockingfile conversion or scoring
    """
    if purpose_of_object == "filter":
        import autogrow.operators.filter.filter_classes.filter_children_classes
        from autogrow.operators.filter.filter_classes.parent_filter_class import (
            ParentFilter as parent_object,
        )
        from autogrow.operators.filter.filter_classes.get_child_filter_class import (
            get_all_subclasses,
        )

    elif purpose_of_object == "parent_pdbqt_converter":
        import autogrow.docking.docking_class.docking_file_conversion
        from autogrow.docking.docking_class.parent_pdbqt_converter import (
            ParentPDBQTConverter as parent_object,
        )
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses
    elif purpose_of_object == "ParentDocking":
        import autogrow.docking.docking_class.docking_class_children
        from autogrow.docking.docking_class.parent_dock_class import (
            ParentDocking as parent_object,
        )
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses
    elif purpose_of_object == "ParentScoring":
        import autogrow.docking.scoring.scoring_classes.scoring_functions
        from autogrow.docking.scoring.scoring_classes.parent_scoring_class import (
            ParentScoring as parent_object,
        )
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses
    else:
        raise Exception("Invalid purpose_of_object")

    children = get_all_subclasses(parent_object)
    child_dict = {}
    for child in children:
        child_object = child()
        child_name = child_object.get_name()
        child_dict[child_name] = child_object

    return child_dict


def get_path_to_custom_script(
    custom_class, param_name: str, type_desc: str, path_prts: List[str]
) -> str:
    cname = custom_class[0]
    py_file_path = custom_class[1]

    if os.path.exists(py_file_path) is False:
        print(custom_class)
        # Check that the path to the original script exists.
        raise Exception(
            f"File can not be found for {param_name} \
            {py_file_path}\n This parameter must be a {type_desc}."
        )

    new_file = os.sep.join(
        [os.path.abspath(os.path.dirname(__file__)), ".."]
        + path_prts
        + [f"{os.path.basename(cname)}.py"]
    )

    if os.path.exists(new_file):
        # File has been copied to proper dir but is not being found by the code
        printout = f"A copy of the custom script {py_file_path} has been moved \
            to {new_file}\n"
        printout += "However, this script could not be imported."
        printout += f"Please check the file naming \
            corresponding to: {custom_class}\n\n"
        print(printout)
        raise Exception(printout)

    return new_file



def copy_new_custom_py_file(src_file: str, dest_file: str):
    print(
        "copying custom class file:\n"
        + f"\t Source: {src_file}\n\t Destination: {dest_file}\n\n"
        + f"Copying is done once, so if the script needs to be changed \
        please either remove or replace the script that is now in the \
        {os.path.dirname(dest_file)} folder.\n"
    )

    copyfile(src_file, dest_file)

    print(
        "\n########################################"
        + "#####################################\n"
        + "AutoGrow has incorporated the custom files.\n"
        + "AutoGrow needs to be restarted and should now be able to run custom scripts.\n"
        + "Please ensure you unit test this code properly before incorporating.\n"
        + "#####################################"
        + "########################################\n"
    )


