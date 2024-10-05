from abc import abstractmethod
import random
from typing import Any, Dict, List, Optional, Tuple, Union, cast
from autogrow.config.argparser import ArgumentVars
from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompoundInfo, ScoreType

class DockingBase(PluginBase):
    def run(self, **kwargs) -> Any:
        """Run the plugin with provided arguments."""
        lig_pdbqt_filename: str = kwargs["lig_pdbqt_filename"]
        file_conversion_class_object: ParentPDBQTConverter = kwargs["file_conversion_class_object"]

        return self.run_docking(
            lig_pdbqt_filename=lig_pdbqt_filename,
            file_conversion_class_object=file_conversion_class_object
        )

    @abstractmethod
    def run_docking(self, lig_pdbqt_filename: str, file_conversion_class_object: ParentPDBQTConverter) -> Optional[str]:
        """
        run_docking is needs to be implemented in each class.

        Inputs:
        :param str pdbqt_filename: a string for docking process raise exception
            if missing

        Returns:
        :returns: str smile_name: name of smiles if it failed to dock returns
            None if it docked properly
        """

        # raise NotImplementedError("run_dock() not implemented")
        pass

    @abstractmethod
    def rank_and_save_output_smi(
        self,
        params: Dict[str, Any],
        current_generation_dir: str,
        current_gen_int: int,
        smiles_file: str,
        deleted_smiles_names_list: List[str],
    ) -> str:
        """
        rank_and_save_output_smi is needs to be implemented in each class.
        raise exception if missing

        Given a folder with PDBQT's, rank all the SMILES based on docking
        score (High to low). Then format it into a .smi file. Then save the
        file.

        Inputs:
        :param dict params: params needs to be threaded here because it has the
            paralizer object which is needed within Scoring.run_scoring_common
        :param str current_generation_dir: path of directory of current
            generation
        :param int current_gen_int: the interger of the current generation
            indexed to zero
        :param str smiles_file: File path for the file with the ligands for the
            generation which will be a .smi file
        :param list deleted_smiles_names_list: list of SMILES which may have
            failed the conversion process

        Returns:
        :returns: str output_ranked_smile_file: the path of the output ranked
            .smi file
        """

        pass

class DockingPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> Any:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        dockings = self.get_selected_plugins_from_params()

        if dockings is None or len(dockings) == 0:
            raise Exception(
                f"You must specify a docking program! Choose from {str(self.plugins.keys())}"
            )
        if len(dockings) > 1:
            raise Exception(
                f"Only one docking program can be selected at a time! You selected {dockings}"
            )

        # Get the selector plugin to use
        docking = cast(DockingBase, self.plugins[dockings[0]])

        docking.run(**kwargs)
