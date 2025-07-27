"""
Utility functions for reaction-based mutation plugins.
"""
import copy
import json
import os
from typing import Any, Dict, Optional

from autogrow.plugins.registry_base import plugin_managers
from autogrow.types import Compound
import autogrow.utils.mol_object_handling as MOH


def validate_rxn_library_path(params: dict):
    """
    Validate the rxn_library_path parameter.

    Args:
        params (dict): A dictionary of parameters to validate.

    Raises:
        ValueError: If rxn_library_path is not provided or is invalid.
    """
    if "rxn_library_path" not in params:
        raise ValueError("rxn_library_path must be provided.")

    if not os.path.exists(params["rxn_library_path"]):
        internal_lib = os.path.join(
            os.path.dirname(__file__),
            "reaction_libraries",
            params["rxn_library_path"],
        )
        if os.path.exists(internal_lib):
            params["rxn_library_path"] = internal_lib
        else:
            raise ValueError(
                "rxn_library_path is not a valid path. "
                "Please provide a valid path to the reaction library."
            )


def load_reaction_library(rxn_library_path: str) -> Dict[str, Any]:
    """
    Load and reformat the reaction library from a JSON file.

    Args:
        rxn_library_path (str): The path to the reaction library directory.

    Returns:
        Dict[str, Any]: The loaded and formatted reaction dictionary.
    """
    rxn_library_file = os.path.join(rxn_library_path, "rxn_library.json")
    try:
        with open(rxn_library_file, "r") as f:
            reaction_dict_raw = json.load(f)
    except Exception as e:
        raise Exception("Failed to load or parse rxn_library.json") from e
    return reformat_rxn_dict(reaction_dict_raw)


def reformat_rxn_dict(old_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively reformat a dictionary loaded from JSON to appropriate types.

    Args:
        old_dict (Dict[str, Any]): The dictionary loaded from JSON.

    Returns:
        Dict[str, Any]: The reformatted dictionary.
    """
    new_dict = {}
    for rxn_key, rxn_dic_old in old_dict.items():
        key_str = str(rxn_key)

        if isinstance(rxn_dic_old, dict):
            new_sub_dict = {}
            for key, item in rxn_dic_old.items():
                sub_key_str = str(key)
                if sub_key_str in [
                    "functional_groups",
                    "group_smarts",
                    "example_rxn_reactants",
                    "reverse_reaction_strings",
                ]:
                    if isinstance(item, list):
                        item = [str(i) for i in item]
                elif sub_key_str == "num_reactants":
                    item = int(item)
                else:
                    item = str(item)

                new_sub_dict[sub_key_str] = item
            new_dict[key_str] = new_sub_dict
        else:
            new_dict[key_str] = str(old_dict[rxn_key])
    return new_dict


def prepare_mol_for_reaction(smiles: str) -> Optional[Any]:
    """
    Prepare a molecule for reaction by sanitizing and protonating it.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        Optional[Any]: An RDKit molecule object if successful, else None.
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit
    try:
        mol = chemtoolkit.mol_from_smiles(smiles, sanitize=False)
    except Exception:
        return None

    mol = MOH.check_sanitization(mol)
    if mol is None:
        return None

    # Reactions generally work better with explicit hydrogens
    mol = MOH.try_reprotanation(mol)
    return mol


def validate_product(
    product_mol: Any, parent_info: Compound, plugin_managers_obj: Any
) -> Optional[str]:
    """
    Validate a reaction product.

    Args:
        product_mol (Any): The RDKit molecule object of the product.
        parent_info (Compound): The parent compound, needed for filters.
        plugin_managers_obj (Any): The plugin managers object.

    Returns:
        Optional[str]: The SMILES string of the validated product, or None.
    """
    product_mol = MOH.check_sanitization(product_mol)
    if product_mol is None:
        return None

    product_mol = MOH.handle_frag_check(product_mol)
    if product_mol is None:
        return None

    product_mol = MOH.check_for_unassigned_atom(product_mol)
    if product_mol is None:
        return None

    product_mol = MOH.try_reprotanation(product_mol)
    if product_mol is None:
        return None

    product_mol = MOH.try_deprotanation(product_mol)
    if product_mol is None:
        return None

    product_mol = MOH.check_sanitization(product_mol)
    if product_mol is None:
        return None

    chemtoolkit = plugin_managers_obj.ChemToolkit.toolkit
    product_smiles: str = chemtoolkit.mol_to_smiles(
        product_mol, isomeric_smiles=True
    )

    if product_smiles == parent_info.smiles:
        return None

    # Run through filters
    tmp_predock_cmpd = Compound(smiles=product_smiles, id="tmp")
    passed_filter = (
        len(plugin_managers_obj.SmilesFilter.run(predock_cmpds=[tmp_predock_cmpd])) > 0
    )

    if passed_filter and len(plugin_managers_obj.DeepFragFilter.plugins) > 0:
        tmp_predock_cmpd.parent_3D_mols = [parent_info.mol_3D]
        passed_filter = (
            len(
                plugin_managers_obj.DeepFragFilter.run(
                    input_params=plugin_managers_obj.DeepFragFilter.params,
                    compounds=[tmp_predock_cmpd],
                )
            )
            > 0
        )

    return product_smiles if passed_filter else None