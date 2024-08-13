import os

def validate_custom_param_type(
    param_name: str, params: dict, outer_type, inner_type, type_desc: str
):
    # Note that this is called from custom_classes.py

    if type(params[param_name]) != outer_type:
        # raise Exception(
        #     "If you want to add custom filters to the filter \
        #     child classes Must be a list of lists \
        #     [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]"
        # )
        print(params[param_name])
        raise Exception(f"Parameter {param_name} must be a {type_desc}.")
    if type(params[param_name][0]) != inner_type:
        print(params[param_name])
        raise Exception(f"Parameter {param_name} must be a {type_desc}.")
        # raise Exception(
        #     "If you want to add custom filters to the filter \
        #     child classes Must be a list of lists \
        #     [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]"
        # )

def validate_custom_params(params):
    # Check if Custom docking option if so there's a few things which need to
    # also be specified. if not lets flag the error.
    if params["dock_choice"] == "Custom":
        if params["docking_executable"] is None:
            raise ValueError(
                "TO USE Custom DOCKING OPTION, MUST SPECIFY THE \
                PATH TO THE docking_executable AND THE DOCKING_CLASS"
            )
        if not os.path.exists(params["docking_executable"]):
            raise ValueError(
                f"Custom docking_executable could not be found at:\
                {params['docking_executable']}"
            )
        if type(
            params["custom_docking_script"]
        ) != list or not os.path.exists(params["custom_docking_script"][1]):
            raise ValueError(
                "TO USE Custom DOCKING OPTION, MUST SPECIFY THE \
                PATH TO THE Custom DOCKING SCRIPT"
            )
    if params["conversion_choice"] == "Custom" and (
        type(params["custom_conversion_script"]) != list
        or not os.path.exists(params["custom_conversion_script"][1])
    ):

        raise ValueError(
            "TO USE Custom conversion_choice OPTION, \
            MUST SPECIFY THE PATH TO THE custom Conversion SCRIPT"
        )

    if params["scoring_choice"] == "Custom" and (
        type(params["custom_scoring_script"]) != list
        or not os.path.exists(params["custom_scoring_script"][1])
    ):

        raise ValueError(
            "TO USE custom scoring_choice OPTION, \
            MUST SPECIFY THE PATH TO THE Custom SCORING SCRIPT"
        )

    # Mutation Settings
    if params["rxn_library"] == "Custom":
        if (
            params["rxn_library_file"] == ""
            or params["function_group_library"] == ""
        ):
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
                 THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
            )
        if not os.path.exists(params["rxn_library_file"]):
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
                THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
            )

        if params["complementary_mol_directory"] == "":
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
                TO THE REACTION LIBRARY USING INPUT PARAMETER function_group_library"
            )
        if not os.path.isdir(params["complementary_mol_directory"]):
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
                TO THE REACTION LIBRARY USING INPUT PARAMETER complementary_mol_directory"
            )
    else:  # Using default settings
        if params["rxn_library_file"] != "":
            raise ValueError(
                "You have selected a Custom rxn_library_file group \
            library but not chosen to use the Custom option for rxn_library. \
            Please use either the provided rxn_library options or chose the Custom \
            option for rxn_library"
            )
        if params["function_group_library"] != "":
            raise ValueError(
                "You have selected a Custom function_group_library but \
            not chosen to use the Custom option for rxn_library. Please use \
            either the provided rxn_library options or chose the Custom option \
            for rxn_library"
            )
        if params["complementary_mol_directory"] != "":
            raise ValueError(
                "You have selected a Custom complementary_mol_directory\
            but not chosen to use the Custom option for rxn_library. \
            Please use either the provided rxn_library options or chose the Custom\
            option for rxn_library"
            )

