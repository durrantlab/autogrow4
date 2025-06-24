"""
Ligand efficiency plugin for re-scoring.
"""
import __future__

from autogrow.plugins.rescoring import RescoringBase
import rdkit  # type: ignore
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from autogrow.types import Compound
from autogrow.plugins.registry_base import plugin_managers

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class LigandEfficiency(RescoringBase):
    """
    Ligand efficiency plugin for re-scoring.
    """

    def run_rescoring(self, docked_cmpds: List[Compound]) -> List[float]:
        """
        Calculate the ligand efficiency and update the docking score. The
        original docking score is saved in the compound's history.

        Args:
            docked_cmpds (List[float]): A list of Compound objects.

        Returns:
            List[float]: A list of the new scores. RescoringBase adds them to
            the compounds.
        """
        rescores = []
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        for compound in docked_cmpds:
            if compound.docking_score is not None:
                mol = chemtoolkit.mol_from_smiles(compound.smiles)
                if mol is not None:
                    heavy_atom_count = chemtoolkit.lipinski_heavy_atom_count(mol)
                    if heavy_atom_count > 0:
                        original_score = compound.docking_score
                        new_score = original_score / heavy_atom_count
                        rescores.append(new_score)
        return rescores

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
                    help="Rescore ligands by ligand efficiency (docking score / number of heavy atoms). The original docking score is saved in the compound's history.",
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
