from typing import Any, Dict, List, Tuple


def setup_filters(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Handle selecting the user-defined Ligand filters.

    This function processes the user parameters to determine which ligand
    filters should be applied. If "No_Filters" is specified and set to True,
    no filters will be used. Otherwise, it calls _picked_filters to determine
    the list of filters to apply.

    Args:
        params (Dict[str, Any]): Dictionary of user variables.

    Returns:
        Dict[str, Any]: Updated dictionary of user variables with the
            chosen_ligand_filters added.
    """
    chosen_ligand_filters, params = _picked_filters(params)
    params["chosen_ligand_filters"] = chosen_ligand_filters

    return params


def _picked_filters(params: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
    """
    Determine the list of filters to apply based on user parameters.

    This function examines the user parameters and compiles a list of filters
    that will be used to check for drug-likeness in each generation. If no
    filters are specified, it defaults to using the LipinskiLenientFilter.

    Args:
        params (Dict[str, Any]): Dictionary of user variables.

    Returns:
        Tuple[List[str], Dict[str, Any]]: A tuple containing:
            - A list of filter names to be applied.
            - The updated params dictionary with filter flags set.
    """
    # TODO: This is not how filters should work anymore. Need to fix this.
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
