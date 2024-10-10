import __future__

from autogrow.plugins.smi_to_sdf import SmiToSdfBase
from autogrow.plugins.smiles_filters import SmilesFilterBase
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
import os


class ObabelSmiToSDF(SmiToSdfBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
        return (
            "SMILES-to-3D-SDF Convertor",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="",  # TODO: Add this.
                )
            ],
        )

    def validate(self, params: dict):
        """Validate the provided arguments."""
        if "obabel_path" not in params:
            raise Exception("You must provide the path to obabel via the --obabel_path parameter.")

    def run_smi_to_sdf_convertor(self, smi_file: str) -> None:
        """
        run_smi_to_sdf_convertor is needs to be implemented in each class.

        Inputs:
        :param str smi_file: The file path to the SMILES file. Note that it is
            a single with with many smi files in it.

        Returns:
        :returns: str: The file path to the 3D SDF file. This one file contains
            multiple molecules (all those that could be successfully converted).
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). If No matches are found to filter list this will
        # return a True as it Passed the filter.

        cmd = f"{self.params['obabel_path']} -ismi {smi_file} -osdf --gen3d --p 7.4 -m -e -O {smi_file}.sdf"

        # print(cmd)

        try:
            os.system(cmd)
        except Exception as e:
            raise Exception(f"Could not convert file with obabel: {smi_file}") from e
