"""
PAINS Filter for removing Pan Assay Interference Compounds.

This module implements the PAINS (Pan Assay Interference Compounds) filter to
eliminate compounds that are likely to interfere with biological assays. It uses
substructure search to identify and filter out these compounds.

This implementation includes PAINS_A, PAINS_B, and PAINS_C filtering.

The filter relies on the RDKit predefined FilterCatalog, which is maintained by
RDKit.

If using the PAINS filter, please cite: Baell JB, Holloway GA. New Substructure
Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening
Libraries and for Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40.
doi:10.1021/jm901137j.
"""
import __future__
from typing import Any, List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound


class PAINSFilter(SmilesFilterBase):
    """
    Implements the PAINS filter for removing Pan Assay Interference Compounds.

    This class filters ligands using PAINS filters, which eliminate Pan Assay
    Interference Compounds using substructure search. It includes PAINS_A,
    PAINS_B, and PAINS_C filtering.

    The filter relies on the RDKit predefined FilterCatalog, which is maintained
    by RDKit.

    If using the PAINS filter, please cite: Baell JB, Holloway GA. New
    Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS)
    from Screening Libraries and for Their Exclusion in Bioassays. J Med Chem 53
    (2010) 2719D40. doi:10.1021/jm901137j.
    """

    def __init__(self) -> None:
        """Initialize the PAINS filter by loading the required filters."""
        self._filters_list = None  # Don't load filters in __init__

    @property
    def filters_list(self):
        """Lazy load filters only when needed."""
        if self._filters_list is None:
            self._filters_list = self.get_filters_list()
        return self._filters_list

    def get_filters_list(self) -> List[Any]:
        """
        Load and return the list of PAINS filters.

        Returns:
            List[FilterCatalog.FilterCatalog]: A list of RDKit PAINS Filters
                including PAINS_A, PAINS_B, PAINS_C, and general PAINS.
        """
        # Make a list of all the different PAINS Filters. PAINS should include
        # PAINS_A,PAINS_B, and PAINS_C, but because RDKit documentation
        # doesn't specify this explicitly we have included all 4 of the PAINS
        # FilterCatalogs for precaution.

        # TODO: This is going to have to be reworked for OpenEye
        from autogrow.plugins.registry_base import plugin_managers

        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        return [
            chemtoolkit.get_pains_a_filter(),
            chemtoolkit.get_pains_b_filter(),
            chemtoolkit.get_pains_c_filter(),
        ]

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the PAINS filter on a given molecule.

        This method applies the PAINS filter to eliminate Pan Assay Interference
        Compounds using substructure search. It checks the molecule against
        PAINS_A, PAINS_B, PAINS_C, and general PAINS filters.

        Based on the PAINS filter implementation in RDKit described in
        http://rdkit.blogspot.com/2016/04/changes-in-201603-release-filtercatalog.html

        Args:
            cmpd (Compound): A Compound to be tested.

        Returns:
            bool: True if the molecule passes all PAINS filters (no matches
                found in any filter list), False otherwise.
        """
        mol = self.cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        # This is our set of all the PAINS filters
        for filters in self.filters_list:

            # If the mol matches a mol in the filter list we return a False
            # (as it failed the filter)
            if filters.HasMatch(mol) is True:
                return False

        # if No matches are found to filter list this will return a True as it
        # Passed the filter.
        return True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the PAINS
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
                    help="PAINS filters against Pan Assay Interference Compounds using \
        substructure a search.",
                )
            ],
        )
