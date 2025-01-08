"""
NIH Filter for removing ligands with undesirable functional groups.

This module implements the NIH filter to eliminate ligands with undesirable
functional groups. It uses the RDKit predefined FilterCatalog, which is
maintained by RDKit.

If using the NIH filter, please cite: Jadhav A, et al. Quantitative Analyses of
Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for
Inhibitors of a Thiol Protease. J Med Chem 53 (2009) 37D51.
doi:10.1021/jm901070c.
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound


from typing import Any, List, Tuple
from autogrow.config.argparser import ArgumentVars


class NIHFilter(SmilesFilterBase):
    """
    Implements NIH filter for removing ligands with bad functional groups.

    This class uses the NIH screening filter to eliminate ligands with
    undesirable functional groups. It relies on the RDKit predefined
    FilterCatalog, which is maintained by RDKit.

    If using the NIH filter, please cite: Jadhav A, et al. Quantitative Analyses
    of Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for
    Inhibitors of a Thiol Protease. J Med Chem 53 (2009) 37D51.
    doi:10.1021/jm901070c.
    """

    def __init__(self) -> None:
        """Initialize the NIH filter by loading the required filters."""
        self._filters = None  # Don't load filters in __init__

    @property
    def filters(self):
        """Lazy load filters only when needed."""
        if self._filters is None:
            self._filters = self.get_filters()
        return self._filters

    def get_filters(self) -> Any:
        """
        Load and return the NIH filters.

        Returns:
            FilterCatalog.FilterCatalog: A set of RDKit NIH Filters.
        """
        # TODO: Seems like this won't work for OpenEye.
        from autogrow.plugins.plugin_manager_instances import plugin_managers

        return plugin_managers.ChemToolkit.toolkit.get_nih_filter()

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the NIH filter on a given molecule.

        This method applies the NIH filter to eliminate ligands with undesirable
        functional groups. It matches common false positive molecules to the
        current mol.

        Based on the PAINS filter implementation in RDKit described in
        http://rdkit.blogspot.com/2016/04/changes-in-201603-release-filtercatalog.html

        Args:
            cmpd (Compound): A Compound to be tested.

        Returns:
            bool: True if the molecule passes the filter (no matches found in
                the filter list), False otherwise.
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). if No matches are found to filter list this will
        # return a True as it Passed the filter.
        mol = self.cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        return self.filters.HasMatch(mol) is not True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the NIH
        Filter. It allows users to enable the filter via command-line options.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of ArgumentVars objects defining the
                command-line arguments.
        """
        return (
            "SMILES Filters",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="NIH filters against molecules with undersirable functional groups \
        using substructure a search.",
                )
            ],
        )
