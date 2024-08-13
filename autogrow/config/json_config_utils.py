import json
import copy
import os
import sys


def convert_json_params_from_unicode(params_unicode):
    """
    Set the parameters that will control this ConfGenerator object.

    :param dict params_unicode: The parameters. A dictionary of {parameter name:
                value}.
    Returns:
    :returns: dict params: Dictionary of User variables
    """
    # Also, rdkit doesn't play nice with unicode, so convert to ascii

    # Because Python2 & Python3 use different string objects, we separate their
    # usecases here.
    params = {}
    if sys.version_info < (3,):
        for param in params_unicode:
            val = params_unicode[param]
            if isinstance(val, unicode):
                val = str(val).encode("utf8")
            key = param.encode("utf8")
            params[key] = val
    else:
        for param in params_unicode:
            val = params_unicode[param]
            key = param
            params[key] = val
    return params

def save_vars_as_json(vars):
    """
    This function saves the vars dictionary as a json file. This can be used
    later to track experiments and is necessary for several of the utility
    scripts.
    It saves all variables except the parallelizer class object.

    It saves the file to the output_directory + "vars.json"
        -If AutoGrow has been run multiple times for the same directory it
        will save the new vars file as append a number to the file name
        starting with 2. The util scripts will only look at the original "vars.json"
            ie) output_directory + "vars_2.json"

    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    """
    output_directory = vars["output_directory"]

    vars_file = output_directory + os.sep + "vars.json"
    if os.path.exists(vars_file):
        # vars.json already exists. lets make the next file.
        path_exists = True
        i = 2
        while path_exists:
            vars_file = f"{output_directory}{os.sep}vars_{i}.json"
            if os.path.exists(vars_file):
                i = i + 1
            else:
                path_exists = False

    temp_vars = {}
    for k in vars.keys():
        if "parallelizer" in k or k == "filter_object_dict":
            continue

        temp_vars[k] = copy.deepcopy(vars[k])

    with open(vars_file, "w") as fp:
        json.dump(temp_vars, fp, indent=4)

