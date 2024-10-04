from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import Dict, List, Optional, Tuple, Union
from autogrow.config.argparser import ArgumentVars, register_argparse_group

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PreDockedCompoundInfo
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
import copy

from autogrow.plugins.plugin_base import PluginBase
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


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
    def run(self, **kwargs) -> List:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        if "smiles" in kwargs:
            # Make sure it is a list of strings
            assert isinstance(kwargs["smiles"], str), "smiles must be a list of strings"
            assert isinstance(
                kwargs["smiles"][0], str
            ), "smiles must be a list of strings"

            # Run filter on a single smiles string.
            return self._run_filter_on_just_smiles(kwargs["smiles"])
        elif "predocked_compounds" in kwargs:
            # Run filter on a list of smiles strings.
            return self._run_filter_on_predocked_compounds(
                kwargs["predocked_compounds"]
            )
        else:
            raise ValueError(
                "No smiles strings or predocked compounds provided to filter!"
            )

    def _run_filter_on_predocked_compounds(
        self, list_of_new_ligands: List[PreDockedCompoundInfo]
    ) -> List[PreDockedCompoundInfo]:
        """
        This will run a filter of the Users chosing.

        This will take a list of lists of ligands to filter. list_of_new_ligands =
        [["CCC","Zinc123],["CCCC","Zinc1234]]

        Inputs:
        :param dict params: User variables which will govern how the programs runs
        :param list list_of_new_ligands: list of lists containing all the newly
            generated ligands and their names

        Returns:
        :returns: list ligands_which_passed_filter: a list of only the molecules
            which passed the filter. Excludes all molecules which failed.
        """

        # make a list of tuples for multi-processing Filter
        job_input = []
        for smiles_info in list_of_new_ligands:
            # TODO: I don't thiknk you need to wrap this in a tuple.
            job_input.append((smiles_info,))
        job_input = tuple(job_input)

        results = self.params["parallelizer"].run(job_input, self._run_filters_mol)

        # remove mols which fail the filter
        return [x for x in results if x is not None]

    def _run_filters_mol(
        self, smiles_info: PreDockedCompoundInfo
    ) -> Optional[PreDockedCompoundInfo]:
        """
        This takes a smiles_string and the selected filter list (child_dict) and
        runs it through the selected filters.

        Inputs:
        :param list smiles_info: A list with info about a ligand, the SMILES string
            is idx=0 and the name/ID is idx=1. example: smiles_info
            ["CCCCCCC","zinc123"]
        :param dict child_dict: This dictionary contains all the names of the
            chosen filters as keys and the the filter objects as the items Or None if
            User specifies no filters

        Returns:
        :returns: list smiles_info: list of the smiles_info if it passed the filter.
            returns None If the mol fails a filter.
        """

        smiles_string = smiles_info.smiles

        mol = Chem.MolFromSmiles(smiles_string, sanitize=False)
        # try sanitizing, which is necessary later
        mol = MOH.check_sanitization(mol)
        if mol is None:
            return None

        mol = MOH.try_deprotanation(mol)
        if mol is None:
            return None

        mol = MOH.check_sanitization(mol)
        if mol is None:
            return None

        # remove charge from mol objects. This affects some properties
        # such as: logP, Mol refractivity, and polar surface area
        # which can impact filters such as Ghose and VandeWaterbeemd
        # This is done because logP is traditionally applied to neutral molecules
        uncharger_obj = rdMolStandardize.Uncharger()
        mol = uncharger_obj.uncharge(mol)
        if mol is None:
            return None

        # run through the filters
        filter_result = self._run_all_selected_filters(mol)

        # see if passed. If it passed return the smiles_info
        return smiles_info if filter_result else None

    def _run_filter_on_just_smiles(self, smiles: str) -> List[str]:
        """
        This takes a smiles_string and the selected filter list (child_dict) and
        runs it through the selected filters.

        Inputs:
        :param str smile_string: A smiles_string. example: smiles_info
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
