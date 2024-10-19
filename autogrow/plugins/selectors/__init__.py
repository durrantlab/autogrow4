from abc import abstractmethod
import random
from typing import Any, List, Optional, Tuple, cast
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompound, ScoreType
from autogrow.utils.logging import LogLevel, log_info


class SelectorBase(PluginBase):
    def run(self, **kwargs) -> Any:
        """Run the plugin with provided arguments."""
        usable_smiles: List[PreDockedCompound] = kwargs["usable_smiles"]
        score_type: ScoreType = kwargs["score_type"]
        num_to_choose: int = kwargs["num_to_choose"]
        favor_most_negative: bool = kwargs["favor_most_negative"]

        return self.run_selector(
            usable_smiles=usable_smiles,
            num_to_choose=num_to_choose,
            score_type=score_type,
            favor_most_negative=favor_most_negative,
        )

    @abstractmethod
    def run_selector(
        self,
        usable_smiles: List[PreDockedCompound],
        num_to_choose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompound]:
        pass

    @abstractmethod
    def finalize_composite_docking_diversity_list(
        self,
        docking_diversity_list: List[PreDockedCompound],
        usable_smiles: List[PreDockedCompound],
    ) -> List[PreDockedCompound]:
        pass

    def get_chosen_mol_full_data_list(
        self,
        chosen_mol_list: List[PreDockedCompound],
        usable_smiles: List[PreDockedCompound],
    ) -> List[PreDockedCompound]:
        """
        This function will take a list of chosen molecules and a list of all the
        SMILES which could have been chosen and all of the information about those
        SMILES (ie. ligand name, SMILES string, docking score, diversity score...)

        It will iterated through the list of chosen mols (chosen_mol_list), get
        all the information from the usable_smiles Then it appends the
        corresponding item in usable_smiles to a new list
        weighted_order_list

        --- an issue to be aware of is that there may be redundancies in both
            chosen_mol_list and usable_smiles this causes a many-to-many
            problem so if manipulating this section you need to solve for
            one-to-many
        ---for this reason if this gets altered it will raise an
            AssertionError if the one-to-many is violated.

        It then shuffles the order of the list which to prevent biasing by the
        order of the ligands.

        It will return that list of the chosen molecules in a randomly shuffled
        order.

        Inputs:
        :param list chosen_mol_list: a list of chosen molecules
        :param list usable_smiles: List of all the possibly chosen ligs
            and all the of the info about it (ie. ligand name, SMILES string, docking
            score, diversity score...) ["CCCC"  "zinc123"   1    -0.1  -0.1]

        Returns:
        :returns: list weighted_order_list: a list of all the SMILES with all of
            the associated information in a random order
        """
        sorted_list = sorted(
            usable_smiles, key=lambda x: x.get_previous_score(ScoreType.DOCKING)
        )
        weighted_order_list: List[PreDockedCompound] = []
        for smile in chosen_mol_list:
            for smile_pair in sorted_list:
                if smile == smile_pair.smiles:
                    weighted_order_list.append(smile_pair)
                    break

        if len(weighted_order_list) != len(chosen_mol_list):
            raise AssertionError(
                "weighted_order_list not the same length as the chosen_mol_list"
            )

        random.shuffle(weighted_order_list)

        return weighted_order_list


class SelectorPluginManager(PluginManagerBase):
    def execute(self, **kwargs) -> List[PreDockedCompound]:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        selectors = self.get_selected_plugins_from_params()

        if selectors is None or len(selectors) == 0:
            raise Exception(
                f"You must specify a selector! Choose from {str(self.plugins.keys())}"
            )
        if len(selectors) > 1:
            raise Exception(
                f"Only one selector can be selected at a time! You selected {selectors}"
            )

        # Get the selector plugin to use
        selector = cast(SelectorBase, self.plugins[selectors[0]])

        docking_fitness_smiles_list: List[PreDockedCompound] = []
        diversity_smile_list: List[PreDockedCompound] = []

        # Select the molecules based on the docking score
        num_usable_smiles = len(kwargs["usable_smiles"])
        num_to_choose = kwargs["num_seed_dock_fitness"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by docking score from {num_usable_smiles} available compounds"
            )
            with LogLevel():
                docking_fitness_smiles_list = selector.run(
                    **{
                        "usable_smiles": kwargs["usable_smiles"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DOCKING,
                        "favor_most_negative": kwargs["favor_most_negative"],
                    }
                )

        # Select the molecules based on the diversity score
        num_to_choose = kwargs["num_seed_diversity"]
        if num_to_choose > 0:
            log_info(
                f"{selector.name}: Selecting {num_to_choose} compounds by diversity score from {num_usable_smiles} available compounds"
            )

            with LogLevel():
                diversity_smile_list = selector.run(
                    **{
                        "usable_smiles": kwargs["usable_smiles"],
                        "num_to_choose": num_to_choose,
                        "score_type": ScoreType.DIVERSITY,
                        "favor_most_negative": kwargs["favor_most_negative"],
                    }
                )

        # Combine the two lists
        docking_diversity_list = list(docking_fitness_smiles_list)
        docking_diversity_list.extend(diversity_smile_list)

        # Finalize the list
        return selector.finalize_composite_docking_diversity_list(
            docking_diversity_list, kwargs["usable_smiles"]
        )
