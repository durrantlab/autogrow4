from abc import abstractmethod
from typing import List, Optional, cast, Type
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import Compound
from autogrow.plugins.plugin_base import PluginBase
from autogrow.utils.logging import log_warning
from autogrow.plugins.registry_base import plugin_managers


class PoseFilterBase(PluginBase):
    """
    This is a script containing all the filters for docked compounds

    1) Filters based on ProLIF
    """

    @property
    def plugin_type_name(self) -> str:
        """Return the user-friendly name of the plugin type."""
        return "Pose Filter"

    @property
    def plugin_description(self) -> str:
        """Return a brief description of the plugin type."""
        return "Filters docked poses based on interactions or other criteria"

    def run(self, **kwargs) -> bool:
        """
        Execute the filter plugin with the provided receptor and molecule.

        This method serves as a standardized entry point for all filter plugins, delegating the actual filtering logic
         to the plugin-specific `run_filter` implementation.

        Args:
        **kwargs:a dictionary of arguments to pass to the plugin. It must contain the path to the
        receptor (receptor_path), a list containing Compound objects that represent docked molecules, and a
        dictionary containing the input parameters specified at the command line

        Returns:
            bool: True if the molecule passes the filter criteria, False
                otherwise
        """
        return self.run_filter(**kwargs)

    @abstractmethod
    def run_filter(self, **kwargs) -> bool:
        """
        run_filter is needs to be implemented in each class.

        Args:
        **kwargs:a dictionary of arguments to pass to the plugin. It must contain the path to the
        receptor (receptor_path), a list containing Compound objects that represent docked molecules, and a
        dictionary containing the input parameters specified at the command line

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass


class PoseFilterPluginManager(PluginManagerBase):
    """Manages and executes interaction-based filter plugins in the autogrow framework."""
    @property
    def default_plugin(self) -> Optional[str]:
        """Return the name of the default plugin for this manager."""
        return None

    def __init__(self, plugin_base_class: Type[PluginBase]):
        self.setup_filter_logger_file = True
        super().__init__(plugin_base_class)

    def execute(self, **kwargs) -> List:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin. It must contain the path to the
        receptor (receptor_path), a list containing Compound objects that represent docked molecules, and a
        dictionary containing the input parameters specified at the command line

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        # Run filter on a single smiles string.
        passed_cmpds: List[Compound] = []
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        receptor = chemtoolkit.mol_from_pdb_file(kwargs["docking_plugin_manager_params"]["receptor_path"], remove_hs=False, sanitize=True)
        for docked_cmpd in kwargs["docked_cmpds"]:
            # run through the filters
            passed = self._run_all_selected_filters(receptor, docked_cmpd, kwargs["docking_plugin_manager_params"])
            if passed:
                passed_cmpds.append(docked_cmpd)

        return passed_cmpds

    def _run_all_selected_filters(self, receptor, docked_cmpd: Compound, docking_plugin_manager_params) -> bool:
        """
        Iterate through all the filters specified by the user for a single  molecule. returns True if the mol passes
        all the chosen filters. Returns False if the mol fails any of the filters.

        Inputs:
        :param receptor: A rdkit object representing the receptor,
        :param docked_cmpd: A Compound object containing the path to the SDF file of the docked molecule
        :param docking_plugin_manager_params: dictionary received as input at the command line

        Returns:
        returns bool: True if the mol passes all the filters. False if the mol
            fails any filters.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        mols = chemtoolkit.mols_from_sdf_file(docked_cmpd.sdf_path, remove_hs=False)
        docked_cmpd_mol = None
        if mols:
            docked_cmpd_mol = mols[0]

        if docked_cmpd_mol is None:
            return False

        filters_failed = 0
        for plugin_name in self.plugins:
            # mol_copy = copy.deepcopy(mol)

            plugin = cast(PoseFilterBase, self.plugins[plugin_name])
            filter_function = plugin.run
            if not filter_function(receptor=receptor, docked_cmpd=docked_cmpd_mol,
                                   docking_plugin_manager_params=docking_plugin_manager_params):
                filters_failed = filters_failed + 1
                log_warning(f"Failed {plugin_name} filter: {docked_cmpd.smiles} (id: {docked_cmpd.id})")
                self.filter_logger_file.info(f"Failed {plugin_name} filter: {docked_cmpd.smiles} (id: {docked_cmpd.id})")
        if filters_failed == 0:
            return True

        # failed one or more filters
        return False
