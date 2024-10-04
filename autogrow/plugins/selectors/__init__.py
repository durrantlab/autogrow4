from abc import abstractmethod
from typing import Any, List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompoundInfo, ScoreType


class SelectorBase(PluginBase):
    def run(self, **kwargs) -> Any:
        """Run the plugin with provided arguments."""
        usable_smiles: List[PreDockedCompoundInfo] = kwargs["usable_smiles"]
        score_type: ScoreType = kwargs["score_type"]
        num_to_chose: int = kwargs["num_to_chose"]
        favor_most_negative: bool = kwargs["favor_most_negative"]

        self.run_selector(
            usable_smiles=usable_smiles,
            num_to_chose=num_to_chose,
            score_type=score_type,
            favor_most_negative=favor_most_negative,
        )

    def throw_exception_if_tourn_size_defined(self):
        if "tourn_size" in self.params:
            raise Exception(
                "You specified the tourn_size parameter, but you are not using the tournament selector."
            )

    @abstractmethod
    def run_selector(
        self,
        usable_smiles: List[PreDockedCompoundInfo],
        num_to_chose: int,
        score_type: ScoreType,
        favor_most_negative: bool = True,
    ) -> List[PreDockedCompoundInfo]:
        pass


class SelectorPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> List[PreDockedCompoundInfo]:
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

        selector = self.plugins[selectors[0]]
        return selector.run(**kwargs)
