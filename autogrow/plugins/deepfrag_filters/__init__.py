"""
Defines the base class for managing the execution of DeepFrag plugins.
"""

from typing import List, cast
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilterBase


class DeepFragFilterManager(PluginManagerBase):
    """
    Manages the execution of DeepFrag plugins.

    This class is responsible for selecting and executing DeepFrag plugins.
    """

    def execute(self, **kwargs) -> List[Compound]:
        """
        Execute the selected DeepFrag plugin with provided arguments.

        Args:
            **kwargs: A dictionary of arguments to pass to the plugin.

        Returns:
            List[Compound]: A list of Compound objects.

        Raises:
            Exception: If multiple DeepFrag filters are selected.
        """
        deepfrag = self.get_selected_plugins_from_params()

        if deepfrag is None or len(deepfrag) == 0:
            return kwargs["compounds"]
        if len(deepfrag) > 1:
            raise Exception(
                f"Only one DeepFrag filter can be selected at a time! You selected {deepfrag}"
            )

        # Get the selector plugin to use
        deepfrag = cast(DeepFragFilterBase, self.plugins[deepfrag[0]])
        return deepfrag.run(**kwargs)
