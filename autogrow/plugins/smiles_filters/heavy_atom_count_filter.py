"""
Implements a filter for heavy-atom count.
"""
import __future__
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.registry_base import plugin_managers
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound

class HeavyAtomCountFilter(SmilesFilterBase):
    """
    A filter that screens compounds based on a minimum and maximum heavy-atom count.
    """

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the
        HeavyAtomCountFilter.

        Returns:
            List[ArgumentVars]: A list of ArgumentVars objects defining the
            command-line arguments.
        """
        return [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Enable filtering of compounds by a heavy-atom count range. This filter removes molecules that have fewer than `min_heavy_atoms` or more than `max_heavy_atoms`.",
            ),
            ArgumentVars(
                name="min_heavy_atoms",
                type=int,
                default=0,
                help="The minimum number of heavy atoms allowed for a compound.",
            ),
            ArgumentVars(
                name="max_heavy_atoms",
                type=int,
                default=50,
                help="The maximum number of heavy atoms allowed for a compound.",
            ),
        ]

    def validate(self, params: dict):
        """
        Validate the arguments provided for the HeavyAtomCountFilter.

        Args:
            params (dict): A dictionary of parameters provided to the plugin.

        Raises:
            ValueError: If 'min_heavy_atoms' or 'max_heavy_atoms' are invalid or inconsistent.
        """
        try:
            self.min_heavy_atoms = int(params["min_heavy_atoms"])
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid value for min_heavy_atoms: {params['min_heavy_atoms']}. Must be an integer."
            ) from e

        try:
            self.max_heavy_atoms = int(params["max_heavy_atoms"])
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid value for max_heavy_atoms: {params['max_heavy_atoms']}. Must be an integer."
            ) from e

        if self.min_heavy_atoms < 0:
            raise ValueError("min_heavy_atoms cannot be negative.")
        if self.max_heavy_atoms < self.min_heavy_atoms:
            raise ValueError("max_heavy_atoms cannot be less than min_heavy_atoms.")

    def run_filter(self, cmpd: Compound) -> bool:
        """
        Run the heavy-atom count filter on a given molecule.

        Args:
            cmpd (Compound): A Compound object to be tested.

        Returns:
            bool: True if the molecule's heavy-atom count is within the
            specified range, False otherwise.
        """
        mol = self.cmpd_to_mol(cmpd)
        if mol is None:
            return False
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        heavy_atom_count = chemtoolkit.lipinski_heavy_atom_count(mol)
        return self.min_heavy_atoms <= heavy_atom_count <= self.max_heavy_atoms