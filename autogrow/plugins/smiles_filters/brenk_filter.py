"""#BRENK filter
This will filter a ligand using the BRENK filter for lead-likeliness, by
matching common false positive molecules to the current mol..

This script relies on the RDKit predefined FilterCatalog. FilterCatalog is
maintained by RDKit.

If using the BRENK filter please cite: Brenk R et al. Lessons Learnt from
Assembling Screening Libraries for Drug Discovery for Neglected Diseases.
ChemMedChem 3 (2008) 435-444. doi:10.1002/cmdc.200700139.
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
import rdkit  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars


class BRENKFilter(SmilesFilterBase):
    """
    This will filter a ligand using a BRENK screening filter for
    lead-likeliness, by matching common false positive molecules to the
    current mol.

    This script relies on the RDKit predefined FilterCatalog. FilterCatalog is
        maintained by RDKit.

    If using the BRENK filter please cite: Brenk R et al. Lessons Learnt from
    Assembling Screening Libraries for Drug Discovery for Neglected Diseases.
    ChemMedChem 3 (2008) 435-444. doi:10.1002/cmdc.200700139.

    Inputs:
    :param class ParentFilter: a parent class to initialize off of.
    """

    def __init__(self) -> None:
        """
        This loads in the filters which will be used.
        """
        self.filters = self.get_filters()

    def get_filters(self) -> FilterCatalog.FilterCatalog:
        """
        This loads in the filters which will be used.

        Returns:
        :returns: rdkit.Chem.rdfiltercatalog.FilterCatalog filters: A set of
            RDKit Filters
        """
        # Make a list of the BRENK filter.
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        # This is our set of all the BRENK filters
        return FilterCatalog.FilterCatalog(params)

    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        Runs a BRENK filter by matching common false positive molecules to the
        current mol. Filters for for lead-likeliness.

        Based on the PAINS filter implementation in RDKit described in
        http://rdkit.blogspot.com/2016/04/changes-in-201603-release-filtercatalog.html

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be
            tested if it passes the filters

        Returns:
        :returns: bool bool: True if the mol passes the filter; False if it
            fails the filter
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). If No matches are found to filter list this will
        # return a True as it Passed the filter.

        return self.filters.HasMatch(mol) is not True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
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
