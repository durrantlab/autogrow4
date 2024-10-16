import __future__

from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
import os


class FakeSmiTo3DSDF(SmiTo3DSdfBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
        return (
            "SMILES-to-3D-SDF Convertor",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Fake conversion to 3D SDF. Converts only to 2D for speed. This is for testing purposes only. Use together with FakeDocking.",
                )
            ],
        )

    def validate(self, params: dict):
        """Validate the provided arguments."""
        if "obabel_path" not in params:
            raise Exception(
                "You must provide the path to obabel via the --obabel_path parameter."
            )

    def run_smi_to_3d_sdf_convertor(
        self, predock_cmpds: List[PreDockedCompound], pwd: str
    ) -> List[PreDockedCompound]:
        """
        run_smi_to_sdf_convertor is needs to be implemented in each class.

        Inputs:
        :param str predock_cmpds: A list of PreDockedCompound objects. Each
            conains a SMILES string, a name, etc.
        :param str pwd: The path to the working directory.

        Returns:
        :returns: List[PreDockedCompound]: A list of PreDockedCompound,
            the same as the input, but with the sdf_3d_path field filled in.
        """

        for cmpd_idx, predock_cmpd in enumerate(predock_cmpds):
            out_file = f"{pwd}compound{cmpd_idx}.sdf"
            predock_cmpd.sdf_3d_path = out_file

            with open(out_file, "w") as f:
                f.write("fake 3D SDF file")

        return predock_cmpds
