"""
This script holds the parent class for filtering.
This is used as the basis for all filter classes.
"""
import __future__
from abc import ABC, abstractmethod
import rdkit


class ParentFilter(ABC):
    """
    This is a script containing all of the filters for drug likeliness

    Filters for orally bio-available drugs:
        1) Lipinski

    Filters for for lead-likeness:
        1) GhoseFilter
        2) GhoseModifiedFilter
        3) MozziconacciFilter

    Filters for CNS/Blood Brain Barrier Permeable:
        1) VandeWaterbeemdFilter

    False-Positive/Metabolite substructure searches:
        1) PAINSFilter
        2) NIHFilter
        3) BRENKFilter
    """

    def get_name(self) -> str:
        """
        Returns the current class name.

        Returns:
        :returns: str self.__class__.__name__: the current class name.
        """

        return self.__class__.__name__

    @abstractmethod
    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        run_filter is needs to be implemented in each class.

        Inputs:
        :param rdkit.Chem.rdchem.Mol mol: a molecule to filter

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        # raise NotImplementedError("run_filter() not implemented")
        pass
