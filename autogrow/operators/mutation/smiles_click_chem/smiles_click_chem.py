"""SMILECLICK Class"""

import __future__

import contextlib
import random
import os
import json
import copy
from typing import Any, Dict, List, Optional, Tuple, Union

from autogrow.plugins.plugin_manager_base import get_plugin_manager
from autogrow.types import PreDockedCompoundInfo

# Disable the unnecessary RDKit warnings
# rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


class SmilesClickChem(object):
    """    This class will take a molecule and Mutate it by reacting it.    """

    # def __init__(
    #     self,
    #     rxn_library_variables: List[str],
    #     list_of_already_made_smiles: List[PreDockedCompoundInfo],
    # ) -> None:
    #     """
    #     init for SmilesClickChem. This will set up all the reaction and
    #     functional dictionaries required to Mutate a molecular

    #     Inputs:
    #     :param list rxn_library_variables: a list of user variables which
    #         define the rxn_library, rxn_library_file,
    #         complementary_mol_directory, and function_group_library. ie.
    #         rxn_library_variables = [params['rxn_library'],
    #         params['rxn_library_file'],
    #         params['function_group_library'],params['complementary_mol_directory']]
    #     :param list list_of_already_made_smiles: a list of lists. Each
    #         sublist contains info about a smiles made in this generation via
    #         mutation ie.[['O=C([O-])',
    #         '(Gen_3_Mutant_37_747+ZINC51)Gen_4_Mutant_15_52']]
    #     """

    #     # Unpackage the rxn_library_variables

    #     # params['rxn_library'], params['rxn_library_file'], params['function_group_library']

    #     rxn_library = rxn_library_variables[0]
    #     rxn_library_file = rxn_library_variables[1]
    #     function_group_library = rxn_library_variables[2]
    #     complementary_mol_dir = rxn_library_variables[3]
    #     self.reaction_dict = self._load_rxn_lib(rxn_library, rxn_library_file)

    #     # Retrieve the dictionary containing
    #     # all the possible ClickChem Reactions
    #     self.list_of_reaction_names = list(self.reaction_dict.keys())

    #     self.functional_group_dict = self._load_functional_grps(
    #         rxn_library, function_group_library
    #     )
    #     self.complementary_mol_dict = self._load_complementary_mols(
    #         rxn_library, complementary_mol_dir
    #     )

    #     # List of already predicted smiles
    #     self.list_of_already_made_smiles = [
    #         x.smiles for x in list_of_already_made_smiles
    #     ]

    def add_mutant_smiles(
        self, list_of_already_made_smiles_infos: List[PreDockedCompoundInfo]
    ) -> None:
        """
        This updates the list of Smiles which have been made in this
        generation via mutation.

        Inputs:
        :param list list_of_already_made_smiles_infos: a list of lists. Each sublist
            contains info about a smiles made in this generation via mutation.
            ie. [['O=C([O-])',
            '(Gen_3_Mutant_37_747+ZINC51)Gen_4_Mutant_15_52']]
        """
        list_of_already_made_smiles = [
            x.smiles for x in list_of_already_made_smiles_infos
        ]
        self.list_of_already_made_smiles.extend(list_of_already_made_smiles)

    def _make_reactant_order_list(
        self,
        substructure_search_results: List[int],
        has_substructure_matches_count: int,
    ) -> List[int]:
        """
        make an ordered list of reactants which composed of 0 and 1. This list
        will be used (in later steps) to determine which reactant is the
        ligand and which requires a complementary molecule.

        Inputs:
        :param list substructure_search_results: list composed of 0 and 1. 1
            for if it has the substructure 0 for not
        :param int has_substructure_matches_count: how many substructure
            matches there are
        Returns:
        :returns: list reactant_order_list: an ordered list of reactants which
            composed of 0 and 1.
        """
        # TODO: Is this used anywhere?

        # create a list to be used to determine which reactants need
        # complementary mol and which will use the Ligand
        reactant_order_list = []

        # for mols w at least 1 substructure
        if has_substructure_matches_count == 1:
            reactant_order_list = substructure_search_results
        elif has_substructure_matches_count > 1:
            # if more than 1 reactant is found in the ligand than we need to
            # randomly pick 1 to be the molecule in the reaction and the
            # other(s) to be mols chosen from the complementary molecule
            # dictionary

            chosen_as_mol_num = random.randint(0, has_substructure_matches_count - 1)
            counter_of_matches = 0
            for substructure_search_result in substructure_search_results:
                if (
                    substructure_search_result == 1
                    and counter_of_matches == chosen_as_mol_num
                ):
                    reactant_order_list.append(1)
                    counter_of_matches = counter_of_matches + 1
                elif substructure_search_result == 1:
                    reactant_order_list.append(0)
                    counter_of_matches = counter_of_matches + 1
                else:
                    reactant_order_list.append(0)
        return reactant_order_list
