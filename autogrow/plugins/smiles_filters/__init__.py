from abc import abstractmethod
from typing import Any, List, Optional, cast

from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound
from autogrow.plugins.registry_base import plugin_managers

from autogrow.plugins.plugin_base import PluginBase
from autogrow.utils.logging import log_warning
import autogrow.utils.mol_object_handling as MOH


class SmilesFilterBase(PluginBase):
    """
    This is a script containing all of the filters for drug likeliness.

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
        """
        Execute the filter plugin with the provided molecule.

        This method serves as a standardized entry point for all filter plugins,
        delegating the actual filtering logic to the plugin-specific
        `run_filter` implementation.

        Args:
            **kwargs: Keyword arguments containing filter parameters, must
                include: cmpd (Compound): The molecule to be
                filtered

        Returns:
            bool: True if the molecule passes the filter criteria, False
                otherwise
        """
        return self.run_filter(kwargs["cmpd"])

    @abstractmethod
    def run_filter(self, cmpd: Compound) -> bool:
        """
        run_filter is needs to be implemented in each class.

        Inputs:
        :param Compound cmpd: a molecule to filter

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass

    def cmpd_to_rdkit_mol(self, cmpd: Compound) -> Optional[Any]:
        """
        Convert a Compound object to an RDKit molecule object.

        Args:
            cmpd (Compound): The Compound object to
                convert.

        Returns:
            Optional[Chem.Mol]: The RDKit molecule object. None if fails.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        mol = chemtoolkit.mol_from_smiles(cmpd.smiles, sanitize=False)
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

        return mol


class SmilesFilterPluginManager(PluginManagerBase):
    """Manages and executes filter plugins in the autogrow framework."""

    def execute(self, **kwargs) -> List:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        # # Make sure it is a list of strings
        # assert isinstance(kwargs["smiles"], str), "smiles must be a list of strings"
        # assert isinstance(kwargs["smiles"][0], str), "smiles must be a list of strings"

        # Run filter on a single smiles string.
        passed_cmpds: List[Compound] = []
        for cmpd in kwargs["predock_cmpds"]:
            # run through the filters
            passed = self._run_all_selected_filters(cmpd)
            if passed:
                passed_cmpds.append(cmpd)

        return passed_cmpds

    def _run_all_selected_filters(self, cmpd: Compound) -> bool:
        """
        Determine if mol passes filters.

        Iterate through all of the filters specified by the user for a single
        molecule. returns True if the mol passes all the chosen filters. returns
        False if the mol fails any of the filters.

        Inputs:
        :param Compound cmpd: An rdkit mol object to be tested
            if it passes the filters

        Returns:
        returns bol bol: True if the mol passes all the filters. False if the mol
            fails any filters.
        """
        filters_failed = 0
        for plugin_name in self.plugins:
            # mol_copy = copy.deepcopy(mol)
            plugin = cast(SmilesFilterBase, self.plugins[plugin_name])
            filter_function = plugin.run
            if not filter_function(cmpd=cmpd):
                filters_failed = filters_failed + 1
                log_warning(f"Failed {plugin_name} filter: {cmpd.smiles}")

        if filters_failed == 0:
            return True

        # failed one or more filters
        return False
