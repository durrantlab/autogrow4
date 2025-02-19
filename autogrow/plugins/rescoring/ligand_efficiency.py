"""
Ligand efficiency plugin for re-scoring.
"""
import __future__

from autogrow.plugins.rescoring import RescoringBase
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
import rdkit.Chem.Descriptors as Descriptors  # type: ignore
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.types import Compound

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class LigandEfficiency(RescoringBase):
    """
    Ligand efficiency plugin for re-scoring.
    """

    def run_rescoring(self, docked_cmpds: List[Compound]) -> List[Compound]:
        """
        Calculate the ligand efficiency

        Args:
            docked_cmpds (List[Compound]): A list of Compound objects.

        Returns:
            List[Compound]: A list of Compound objects, each
            containing the re-scoring result.
        """
        for compound in docked_cmpds:
            compound.docking_score = compound.docking_score / rdkit.Chem.Lipinski.HeavyAtomCount(
                Chem.MolFromSmiles(compound.smiles))
        return docked_cmpds

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific for the
        ligand efficiency calculation.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of ArgumentVars objects defining the
                command-line arguments.
        """
        return (
            "Ligand efficiency",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Calculates the ligand efficiency.",
                )
            ],
        )

    def validate(self, params: dict):
        """
        Validate the provided arguments for Vina-like docking.

        Args:
            params (dict): A dictionary of plugin parameters to validate.
        """
        pass
