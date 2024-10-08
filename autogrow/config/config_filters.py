import sys
from typing import Any, Dict, List, Tuple


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

    # if there is no user specified ligand filters but they haven't set
    # filters to None ---> set filter to default of LipinskiLenientFilter.
    if len(filter_list) == 0:
        params["LipinskiLenientFilter"] = True
        filter_list.append("LipinskiLenientFilter")

    return filter_list, params
