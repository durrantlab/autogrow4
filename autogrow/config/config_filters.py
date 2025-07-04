import sys
from typing import Any, Dict, List, Tuple

from autogrow.config.config_utils import (
    copy_new_custom_py_file,
    get_path_to_custom_script,
    make_complete_children_dict,
)
from autogrow.validation.validate_custom_classes import validate_custom_param_type


def setup_filters(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    This function handles selecting the user defined Ligand filters.

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: dict params: Dictionary of User variables with the
        chosen_ligand_filters added
    """
    if "No_Filters" in list(params.keys()):
        if params["No_Filters"] is True:
            chosen_ligand_filters = None
        else:
            chosen_ligand_filters, params = _picked_filters(params)
    else:
        chosen_ligand_filters, params = _picked_filters(params)
    params["chosen_ligand_filters"] = chosen_ligand_filters

    import autogrow.operators.filter.execute_filters as Filter

    # get child filter class object function dictionary
    params["filter_object_dict"] = Filter.make_run_class_dict(chosen_ligand_filters)

    return params


def _picked_filters(params: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
    """
    This will take the user params and return a list of the filters
    which a molecule must pass to move into the next generation.

    Inputs:
    :param dict params: Dictionary of User variables
    Returns:
    :returns: list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow
        the same nomenclature and enter that filter in the user
        params["alternative_filters"] as the name of that class and place
        that file in the same folder as the other filter classes.
    """
    filter_list = []
    vars_keys = list(params.keys())

    if "LipinskiStrictFilter" in vars_keys:
        if params["LipinskiStrictFilter"] is True:
            filter_list.append("LipinskiStrictFilter")
    else:
        params["LipinskiStrictFilter"] = False

    if "LipinskiLenientFilter" in vars_keys:
        if params["LipinskiLenientFilter"] is True:
            filter_list.append("LipinskiLenientFilter")
    else:
        params["LipinskiLenientFilter"] = False

    if "GhoseFilter" in vars_keys:
        if params["GhoseFilter"] is True:
            filter_list.append("GhoseFilter")
    else:
        params["GhoseFilter"] = False

    if "GhoseModifiedFilter" in vars_keys:
        if params["GhoseModifiedFilter"] is True:
            filter_list.append("GhoseModifiedFilter")
    else:
        params["GhoseModifiedFilter"] = False

    if "MozziconacciFilter" in vars_keys:
        if params["MozziconacciFilter"] is True:
            filter_list.append("MozziconacciFilter")
    else:
        params["MozziconacciFilter"] = False

    if "VandeWaterbeemdFilter" in vars_keys:
        if params["VandeWaterbeemdFilter"] is True:
            filter_list.append("VandeWaterbeemdFilter")
    else:
        params["VandeWaterbeemdFilter"] = False

    if "PAINSFilter" in vars_keys:
        if params["PAINSFilter"] is True:
            filter_list.append("PAINSFilter")
    else:
        params["PAINSFilter"] = False

    if "NIHFilter" in vars_keys:
        if params["NIHFilter"] is True:
            filter_list.append("NIHFilter")
    else:
        params["NIHFilter"] = False

    if "BRENKFilter" in vars_keys:
        if params["BRENKFilter"] is True:
            filter_list.append("BRENKFilter")
    else:
        params["BRENKFilter"] = False

    if "alternative_filter" in vars_keys:
        filter_list = _handle_alternative_filters(params, filter_list)
    else:
        params["alternative_filter"] = None

    # if there is no user specified ligand filters but they haven't set
    # filters to None ---> set filter to default of LipinskiLenientFilter.
    if len(filter_list) == 0:
        params["LipinskiLenientFilter"] = True
        filter_list.append("LipinskiLenientFilter")

    return filter_list, params


def _handle_alternative_filters(
    params: Dict[str, Any], filter_list: List[str]
) -> List[str]:
    """
    This will handle Custom Filters

    Inputs:
    :param dict params: Dictionary of User variables
    :param list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow the same
        nomenclature and enter that filter in the user params["alternative_filters"]
        as the name of that class and place that file in the same folder as the
        other filter classes.

    Returns:
    :returns: list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow the same
        nomenclature and enter that filter in the user params["alternative_filters"]
        as the name of that class and place that file in the same folder as the
        other filter classes.
    """
    if params["alternative_filter"] is None:
        return filter_list

    validate_custom_param_type(
        "alternative_filter",
        params,
        list,
        list,
        "list of lists [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]",
    )

    full_children_dict = make_complete_children_dict("filter")

    scripts_to_copy = []
    for custom_class in params["alternative_filter"]:
        if custom_class[0] not in full_children_dict.keys():
            new_file = get_path_to_custom_script(
                custom_class,
                "alternative_filter",
                "list of lists [[name_filter1, Path/to/name_filter1.py], [name_filter2, Path/to/name_filter2.py]]",
                ["operators", "filter", "filter_classes", "filter_children_classes"],
            )

            # Add to list of scripts to copy into the filter folder
            scripts_to_copy.append([custom_class[1], new_file])
        filter_list.append(custom_class[0])
    if scripts_to_copy:
        for filter_info in scripts_to_copy:
            copy_new_custom_py_file(filter_info[0], filter_info[1])

        # Technically Exit intentionally but maybe should be a raise Exception
        sys.exit(0)
    return filter_list
