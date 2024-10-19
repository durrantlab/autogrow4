"""#No filter for SMILES
Useful for skipping the filter step.
"""
import __future__

from autogrow.plugins.smiles_filters import SmilesFilterBase
import rdkit  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars


class NoSmilesFilter(SmilesFilterBase):
    """
    This essentially skip the smiles filtering step.

    Inputs:
    :param class ParentFilter: a parent class to initialize off of.
    """

    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        # All compounds pass
        return True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
        return (
            "SMILES Filters",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Apply no SMILES filter. Useful to skip the SMILES-filter step. NOTE: You can also simply omit all SMILES filters from the user parameters to achieve the same effect.",
                )
            ],
        )
