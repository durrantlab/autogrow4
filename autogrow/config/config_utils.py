from typing import Any, Dict, List
import os
from shutil import copyfile


def make_complete_children_dict(purpose_of_object: str) -> Dict[str, Any]:
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
    if purpose_of_object == "parent_pdbqt_converter":
        import autogrow.docking.docking_class.docking_file_conversion
        from autogrow.docking.docking_class.parent_pdbqt_converter import (
            ParentPDBQTConverter as parent_object,
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
