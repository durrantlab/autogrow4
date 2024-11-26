"""
Implements the BRENK filter for lead-likeness screening in AutoGrow.

This module provides the BRENKFilter class that filters ligands using the BRENK
screening filter for lead-likeness. It matches common false positive molecules
to the current molecule.

The BRENK filter relies on the RDKit predefined FilterCatalog. FilterCatalog is
maintained by RDKit.

If using the BRENK filter, please cite: Brenk R et al. Lessons Learnt from
Assembling Screening Libraries for Drug Discovery for Neglected Diseases.
ChemMedChem 3 (2008) 435-444. doi:10.1002/cmdc.200700139.
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound
import rdkit  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars


class BRENKFilter(SmilesFilterBase):
    """
    A filter that uses the BRENK screening filter for lead-likeness.

    This class implements a filter that matches common false positive molecules
    to the current molecule using the BRENK screening filter. It relies on the
    RDKit predefined FilterCatalog.

    Attributes:
        filters (FilterCatalog.FilterCatalog): A set of RDKit BRENK filters.
    """

    def __init__(self) -> None:
        """
        Initialize the BRENKFilter by loading the BRENK filters.
        """
        self.filters = self.get_filters()

    def get_filters(self) -> FilterCatalog.FilterCatalog:
        """
        Load the BRENK filters from RDKit.

        Returns:
            FilterCatalog.FilterCatalog: A set of RDKit BRENK filters.
        """
        # Make a list of the BRENK filter.
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        # This is our set of all the BRENK filters
        return FilterCatalog.FilterCatalog(params)

    def run_filter(self, predock_cmpd: Compound) -> bool:
        """
        Run the BRENK filter on a given molecule.

        This method matches common false positive molecules to the current
        molecule to filter for lead-likeness. It's based on the PAINS filter
        implementation in RDKit described in:
        http://rdkit.blogspot.com/2016/04/changes-in-201603-release-filtercatalog.html

        Args:
            predock_cmpd (PostDockedCompound): A PostDockedCompound to be tested.

        Returns:
            bool: True if the molecule passes the filter; False if it fails.

        Note:
            If the molecule matches any molecule in the filter list, it fails
            the filter (returns False). If no matches are found, it passes the
            filter (returns True).
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). If No matches are found to filter list this will
        # return a True as it Passed the filter.

        mol = self.predock_cmpd_to_rdkit_mol(predock_cmpd)
        if mol is None:
            return False

        return self.filters.HasMatch(mol) is not True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the BRENK filter.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("SMILES Filters")
                - A list with one ArgumentVars object defining the argument to
                  enable the BRENK filter
        """
        return (
            "SMILES Filters",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="BRENK filter for lead-likeliness, by matching common false positive \
        molecules to the current mol.",
                )
            ],
        )
