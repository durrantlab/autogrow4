"""
This script holds the parent class for scoring/rescoring.
This is used as the basis for all scoring/rescoring classes.
"""
import __future__
from abc import ABC, abstractmethod
from typing import Union

from autogrow.types import PostDockedCompoundInfo


class ParentScoring(ABC):
    """
    This is a script containing all of the scoring functions.
    """

    def get_name(self) -> str:
        """
        Returns the current class name.

        Returns:
        :returns: str self.__class__.__name__: the current class name.
        """
        return self.__class__.__name__

    @abstractmethod
    def run_scoring(self, file_path: str) -> Union[PostDockedCompoundInfo, None]:
        """
        run_scoring is needs to be implemented in each class.

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list lig_info: a list containing all info from
            self.smiles_dict for a given ligand and the ligands short_id_name and
            the docking score from the best pose.
        """

        # raise NotImplementedError("run_scoring() not implemented")
        pass
