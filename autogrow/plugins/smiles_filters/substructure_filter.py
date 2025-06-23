"""Plugin for filtering compounds based on a required substructure."""
from typing import List, Optional, Tuple

from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import Compound


class SubstructureFilter(SmilesFilterBase):
    """Filter that passes compounds containing a specified substructure."""

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments for the substructure filter.
        
        Returns:
            Tuple[str, List[ArgumentVars]]: Title and list of arguments.
        """
        return (
            "SMILES Filters",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Filter compounds based on a required substructure.",
                ),
                ArgumentVars(
                    name="--substructure_smiles",  # Note actually SMARTS string
                    help="SMARTS string of required substructure. Compounds will only "
                         "pass if they contain this substructure. Used with the " + self.name + " plugin.",
                    type=str,
                    default=None
                ),
                ArgumentVars(
                    name="--exclude_substructure",  # NOTE: Untested
                    action="store_true",
                    default=False,
                    help="Flag to exclude compounds containing the substructures specified via --substructure_smiles. "
                    "Used with the "
                    + self.name
                    + " plugin.",
                ),
            ]
        )

    def validate(self, params: dict):
        """Validate the substructure SMILES parameter.
        
        Args:
            params (dict): Parameter dictionary.
            
        Raises:
            ValueError: If substructure_smiles is not provided or invalid.
        """
        if "substructure_smiles" not in params or not params["substructure_smiles"]:
            raise ValueError("Must provide --substructure_smiles parameter")

        # Verify the SMILES can be converted to a valid molecule
        assert self.plugin_managers is not None, "Plugin managers not set"
        assert self.plugin_managers.ChemToolkit is not None, "ChemToolkit not set"
        chemtoolkit = self.plugin_managers.ChemToolkit.toolkit
        query_mol = chemtoolkit.mol_from_smarts(params["substructure_smiles"])
        if query_mol is None:
            raise ValueError(
                f"Invalid substructure SMILES: {params['substructure_smiles']}"
            )
        
        self.query_mol = query_mol
        self.exclude_substructure = params.get("exclude_substructure", False)

    def run_filter(self, cmpd: Compound) -> bool:
        """Check if compound contains the required substructure.
        
        Args:
            cmpd (Compound): Compound to check.
            
        Returns:
            bool: True if compound contains substructure, False otherwise.
        """
        # Convert SMILES to mol object for substructure matching
        mol = self.cmpd_to_rdkit_mol(cmpd)
        if mol is None:
            return False

        # Use the chemistry toolkit to do substructure matching
        assert self.plugin_managers is not None, "Plugin managers not set"
        assert self.plugin_managers.ChemToolkit is not None, "ChemToolkit not set"
        chemtoolkit = self.plugin_managers.ChemToolkit.toolkit
        matches = chemtoolkit.get_substruct_matches(mol, self.query_mol)

        return len(matches) == 0 if self.exclude_substructure else len(matches) > 0