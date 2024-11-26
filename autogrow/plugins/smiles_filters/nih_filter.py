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
import rdkit  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore


from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars


class NIHFilter(SmilesFilterBase):
    """
    Implements the NIH filter for removing ligands with undesirable functional
    groups.

    This class uses the NIH screening filter to eliminate ligands with
    undesirable functional groups. It relies on the RDKit predefined
    FilterCatalog, which is maintained by RDKit.

    If using the NIH filter, please cite: Jadhav A, et al. Quantitative Analyses
    of Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for
    Inhibitors of a Thiol Protease. J Med Chem 53 (2009) 37D51.
    doi:10.1021/jm901070c.
    """

    def __init__(self) -> None:
        """
        Initialize the NIH filter by loading the required filters.
        """
        self.filters = self.get_filters()

    def get_filters(self) -> FilterCatalog.FilterCatalog:
        """
        Load and return the NIH filters.

        Returns:
            FilterCatalog.FilterCatalog: A set of RDKit NIH Filters.
        """
        # Make a list of the NIH filter.
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        # This is our set of all the NIH filters
        return FilterCatalog.FilterCatalog(params)

    def run_filter(self, predock_cmpd: Compound) -> bool:
        """
        Run the NIH filter on a given molecule.

        This method applies the NIH filter to eliminate ligands with undesirable
        functional groups. It matches common false positive molecules to the
        current mol.

        Based on the PAINS filter implementation in RDKit described in
        http://rdkit.blogspot.com/2016/04/changes-in-201603-release-filtercatalog.html

        Args:
            predock_cmpd (PostDockedCompound): A PostDockedCompound to be tested.

        Returns:
            bool: True if the molecule passes the filter (no matches found in
                the filter list), False otherwise.
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). if No matches are found to filter list this will
        # return a True as it Passed the filter.
        mol = self.predock_cmpd_to_rdkit_mol(predock_cmpd)
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
