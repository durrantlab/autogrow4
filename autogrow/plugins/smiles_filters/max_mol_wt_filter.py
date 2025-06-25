"""
Implements a filter for maximum molecular weight.
"""
import __future__
from typing import List, Tuple

from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.registry_base import plugin_managers
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound


class MaxMolWtFilter(SmilesFilterBase):
    """
    A filter that screens compounds based on a maximum molecular weight.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the
        MaxMolWtFilter.

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
                    help="Filter compounds by a maximum molecular weight.",
                ),
                ArgumentVars(
                    name="max_mol_weight",
                    type=float,
                    default=550.0,
                    help=f"The maximum molecular weight allowed for a compound. Used with the {self.name} plugin.",
                ),
            ],
        )

    def validate(self, params: dict):
        """
        Validate the arguments provided for the MaxMolWtFilter.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.

        Raises:
            ValueError: If 'max_mol_weight' is not provided or is not a
                valid number.
        """
        if "max_mol_weight" not in params or params["max_mol_weight"] is None:
            raise ValueError(
                f"The 'max_mol_weight' parameter must be provided when using the {self.name} filter."
            )

        try:
            self.max_mw = float(params["max_mol_weight"])
        except (ValueError, TypeError):
            raise ValueError(
                f"Invalid value for max_mol_weight: {params['max_mol_weight']}. Must be a number."
            )

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the maximum molecular weight filter on a given molecule.

        Args:
            cmpd (Compound): A Compound object to be tested.

        Returns:
            bool: True if the molecule's molecular weight is less than or
                equal to the maximum, False otherwise.
        """
        mol = self.cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        exact_mwt = chemtoolkit.descriptors_exact_mol_wt(mol)

        return exact_mwt <= self.max_mw