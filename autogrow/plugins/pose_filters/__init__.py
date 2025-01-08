from abc import abstractmethod
from typing import List, cast
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound
from rdkit import Chem  # type: ignore
from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore
from autogrow.plugins.plugin_base import PluginBase


class PoseFilterBase(PluginBase):
    """
    This is a script containing all the filters for docked compounds

    1) Filters based on ProLIF
    """

    def run(self, **kwargs) -> bool:
        """
        Execute the filter plugin with the provided receptor and molecule.

        This method serves as a standardized entry point for all filter plugins,
        delegating the actual filtering logic to the plugin-specific
        `run_filter` implementation.

        Args:
            **kwargs: Keyword arguments containing filter parameters, must
                include: receptor: A rdkit object representing the receptor,
                docked_cmpd: A rdkit object representing the docked molecule to be filtered

        Returns:
            bool: True if the molecule passes the filter criteria, False
                otherwise
        """
        return self.run_filter(kwargs["receptor"], kwargs["docked_cmpd"])

    @abstractmethod
    def run_filter(self, receptor, docked_cmpd) -> bool:
        """
        run_filter is needs to be implemented in each class.

        Inputs:
        :param receptor: A rdkit object representing the receptor,
        :param docked_cmpd: A rdkit object representing the docked molecule to be filtered

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class PoseFilterPluginManager(PluginManagerBase):
    def execute(self, **kwargs) -> List:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin. It must contain the path to the
        receptor (receptor_path), and a list containing Compound objects that represent docked molecules.

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        # Run filter on a single smiles string.
        passed_cmpds: List[Compound] = []
        receptor = Chem.MolFromPDBFile(kwargs["receptor_path"], sanitize=True)
        for docked_cmpd in kwargs["docked_cmpds"]:
            # run through the filters
            passed = self._run_all_selected_filters(receptor, docked_cmpd)
            if passed:
                passed_cmpds.append(docked_cmpd)

        return passed_cmpds

    def _run_all_selected_filters(self, receptor, docked_cmpd: Compound) -> bool:
        """
        Iterate through all the filters specified by the user for a single
        molecule. returns True if the mol passes all the chosen filters. returns
        False if the mol fails any of the filters.

        Inputs:
        :param receptor: A rdkit object representing the receptor,
        :param docked_cmpd: A Compound object containing the path to the SDF file of the docked molecule

        Returns:
        returns bool: True if the mol passes all the filters. False if the mol
            fails any filters.
        """
        r = Chem.SDMolSupplier(docked_cmpd.sdf_path)
        for lig in r:
            docked_cmpd_mol = lig
            break
        r.reset()

        filters_failed = 0
        for plugin_name in self.plugins:
            # mol_copy = copy.deepcopy(mol)
            plugin = cast(PoseFilterBase, self.plugins[plugin_name])
            filter_function = plugin.run
            if not filter_function(receptor=receptor, docked_cmpd=docked_cmpd_mol):
                filters_failed = filters_failed + 1
                print(f"Failed {plugin_name} filter: {docked_cmpd.smiles}")

        if filters_failed == 0:
            return True

        # failed one or more filters
        return False
