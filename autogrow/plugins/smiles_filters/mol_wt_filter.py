"""
Implements a filter for molecular weight range.
"""
import __future__
from typing import List, Tuple

from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.registry_base import plugin_managers
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound

class MolWtFilter(SmilesFilterBase):
    """
    A filter that screens compounds based on a minimum and maximum molecular weight.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the
        MolWtFilter.

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
                    help="Enable filtering of compounds by a molecular weight range.",
                ),
                ArgumentVars(
                    name="min_mol_weight",
                    type=float,
                    default=0.0,
                    help="The minimum molecular weight allowed for a compound.",
                ),
                ArgumentVars(
                    name="max_mol_weight",
                    type=float,
                    default=550.0,
                    help="The maximum molecular weight allowed for a compound.",
                ),
            ],
        )

    def validate(self, params: dict):
        """
        Validate the arguments provided for the MolWtFilter.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.

        Raises:
            ValueError: If 'min_mol_weight' or 'max_mol_weight' are invalid or inconsistent.
        """
        try:
            self.min_mw = float(params["min_mol_weight"])
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid value for min_mol_weight: {params['min_mol_weight']}. Must be a number."
            ) from e

        try:
            self.max_mw = float(params["max_mol_weight"])
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid value for max_mol_weight: {params['max_mol_weight']}. Must be a number."
            ) from e

        if self.min_mw < 0:
            raise ValueError("min_mol_weight cannot be negative.")
        if self.max_mw < self.min_mw:
            raise ValueError("max_mol_weight cannot be less than min_mol_weight.")

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the molecular weight filter on a given molecule.

        Args:
            cmpd (Compound): A Compound object to be tested.

        Returns:
            bool: True if the molecule's molecular weight is within the
                specified range, False otherwise.
        """
        mol = self.cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        exact_mwt = chemtoolkit.descriptors_exact_mol_wt(mol)
        return self.min_mw <= exact_mwt <= self.max_mw