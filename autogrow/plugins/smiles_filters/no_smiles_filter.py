"""
No SMILES Filter module.

This module provides a dummy filter that allows skipping the SMILES filtering
step in the pipeline.
"""
import __future__

# TODO: Needed?

from autogrow.plugins.smiles_filters import SmilesFilterBase
import rdkit  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars


class NoSmilesFilter(SmilesFilterBase):
    """
    A dummy filter class that skips the SMILES filtering step.

    This class essentially allows all compounds to pass through without any
    filtering, effectively skipping the SMILES filtering step in the pipeline.
    """

    def run_filter(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        Run the dummy filter on a given molecule.

        This method always returns True, allowing all compounds to pass.

        Args:
            mol (rdkit.Chem.rdchem.Mol): An RDKit mol object to be "filtered".

        Returns:
            bool: Always returns True, indicating that all compounds pass.
        """
        return True

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the No SMILES
        Filter. It allows users to explicitly choose to skip the SMILES
        filtering step via command-line options.

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
                    help="Apply no SMILES filter. Useful to skip the SMILES-filter step. NOTE: You can also simply omit all SMILES filters from the user parameters to achieve the same effect.",
                )
            ],
        )
