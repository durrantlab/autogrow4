from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import Dict, List, Optional, Tuple, Union

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompound
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy

from autogrow.plugins.plugin_base import PluginBase
import autogrow.utils.mol_object_handling as MOH


class SmilesFilterBase(PluginBase):
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

    def run(self, **kwargs) -> bool:
        """Run the plugin(s) with provided arguments."""
        return self.run_filter(kwargs["mol"])

    @abstractmethod
    def run_filter(self, mol: Chem.rdchem.Mol) -> bool:
        """
        run_filter is needs to be implemented in each class.

        Inputs:
        :param rdkit.Chem.rdchem.Mol mol: a molecule to filter

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class SmilesFilterPluginManager(PluginManagerBase):
    def execute(self, **kwargs) -> List:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        # Make sure it is a list of strings
        assert isinstance(kwargs["smiles"], str), "smiles must be a list of strings"
        assert isinstance(
            kwargs["smiles"][0], str
        ), "smiles must be a list of strings"

        # Run filter on a single smiles string.
        return self._run_filter_on_just_smiles(kwargs["smiles"])

    def _run_filter_on_just_smiles(self, smiles: str) -> List[str]:
        """
        This takes a smiles and the selected filter list (child_dict) and
        runs it through the selected filters.

        Inputs:
        :param str smile_string: A smiles. example: smiles_info
            ["CCCCCCC","zinc123"]
        :param dict child_dict: This dictionary contains all the names of the
            chosen filters as keys and the the filter objects as the items Or None if
            User specifies no filters

        Returns:
        :returns: str smile_string: smile_string if it passed the filter. returns
            False If the mol fails a filter.
        """

        passed_smiles: List[str] = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            # try sanitizing, which is necessary later
            mol = MOH.check_sanitization(mol)
            if mol is None:
                continue

            mol = MOH.try_deprotanation(mol)
            if mol is None:
                continue

            # run through the filters
            passed = self._run_all_selected_filters(mol)
            if passed:
                passed_smiles.append(smi)

        return passed_smiles

    def _run_all_selected_filters(self, mol: Chem.rdchem.Mol) -> bool:
        """
        Iterate through all of the filters specified by the user for a single
        molecule. returns True if the mol passes all the chosen filters. returns
        False if the mol fails any of the filters.

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be tested
            if it passes the filters
        :param dict child_dict: This dictionary contains all the names of the
            chosen filters as keys and the the filter objects as the items

        Returns:
        returns bol bol: True if the mol passes all the filters. False if the mol
            fails any filters.
        """

        filters_failed = 0
        mol = MOH.check_sanitization(mol)
        if mol is None:
            return False
        for plugin_name in self.plugins:
            mol_copy = copy.deepcopy(mol)
            plugin = self.plugins[plugin_name]
            filter_function = plugin.run
            if not filter_function(mol=mol_copy):
                filters_failed = filters_failed + 1
                print(f"Failed {plugin_name} filter: {Chem.MolToSmiles(mol_copy)}")

        if filters_failed == 0:
            return True

        # failed one or more filters
        return False
