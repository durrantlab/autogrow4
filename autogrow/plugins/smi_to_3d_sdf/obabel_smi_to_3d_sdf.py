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


class ObabelSmiTo3DSDF(SmiTo3DSdfBase):
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
            raise Exception(
                "You must provide the path to obabel via the --obabel_path parameter."
            )

    def run_smi_to_3d_sdf_convertor(
        self, predock_cmpd: PreDockedCompound, pwd: str, cmpd_idx: int
    ) -> PreDockedCompound:
        """
        run_smi_to_sdf_convertor is needs to be implemented in each class.

        Inputs:
        :param str predock_cmpd: A PreDockedCompound object. Conains a
            SMILES string, a name, etc.
        :param str pwd: The path to the working directory.
        :param int cmpd_idx: The index of the compound in the generation.

        Returns:
        :returns: PreDockedCompound: A PreDockedCompound,
            the same as the input, but with the sdf_3d_path field filled in.
        """
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). If No matches are found to filter list this will
        # return a True as it Passed the filter.

        # TODO: Can run on multiple processors

        base_file = f"{pwd}compound{cmpd_idx}"
        in_file = f"{base_file}.smi"
        out_file = f"{base_file}.sdf"
        obabel_path = self.params["obabel_path"]

        with open(in_file, "w") as f:
            f.write(predock_cmpd.smiles)

        convert_success = obabel_convert(
            in_file, out_file, obabel_path, prot_and_3d=True
        )

        if not convert_success:
            return predock_cmpd

        # cmd = f'{obabel_path} -:"{predock_cmpd.smiles}" -osdf --gen3d --p 7.4 -e -O "{out_file}"'

        # try:
        #     os.system(cmd)
        # except Exception as e:
        #     log_warning(
        #         f"Could not convert smiles to 3D SDF with obabel: {predock_cmpd.smiles}"
        #     )
        #     return predock_cmpd

        # if not os.path.exists(out_file):
        #     log_warning(
        #         f"Could not convert smiles to 3D SDF with obabel: {predock_cmpd.smiles}"
        #     )
        #     return predock_cmpd

        # with open(out_file, "r") as f:
        #     content = f.read().strip()
        #     if content == "":
        #         log_warning(
        #             f"Could not convert smiles to 3D SDF with obabel: {predock_cmpd.smiles}"
        #         )
        #         return predock_cmpd

        predock_cmpd.sdf_3d_path = out_file

        return predock_cmpd
