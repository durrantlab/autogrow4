import __future__

# from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.smi_to_3d_sdf import SmiTo3DSdfBase
from autogrow.plugins.smiles_filters import SmilesFilterBase
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import log_warning
from autogrow.utils.obabel import obabel_convert, obabel_convert_cmd
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
        # If the mol matches a mol in the filter list. we return a False (as it
        # failed the filter). If No matches are found to filter list this will
        # return a True as it Passed the filter.

        cmds = []
        out_files = []
        for cmpd_idx, predock_cmpd in enumerate(predock_cmpds):
            base_file = f"{pwd}compound{cmpd_idx}"
            in_file = f"{base_file}.smi"
            out_file = f"{base_file}.sdf"
            obabel_path = self.params["obabel_path"]

            with open(in_file, "w") as f:
                f.write(predock_cmpd.smiles)

            cmd = obabel_convert_cmd(
                in_file, out_file, obabel_path, extra_params="--gen3d --p 7.4"
            )

            cmds.append(cmd)
            out_files.append(out_file)

        # Get parallelizer plugin to use
        # TODO: Need to specify nprocs?
        assert self.plugin_managers is not None, "Plugin managers is None"
        assert (
            self.plugin_managers.ShellParallelizer is not None
        ), "Shell parallelizer is None"
        self.plugin_managers.ShellParallelizer.run(cmds=cmds)

        for cmpd_idx, predock_cmpd in enumerate(predock_cmpds):
            out_file = out_files[cmpd_idx]

            if not os.path.exists(out_file):
                log_warning(
                    f"Could not convert smiles to 3D SDF with obabel: {predock_cmpd.smiles}"
                )
                continue

            with open(out_file, "r") as f:
                content = f.read().strip()
                if content == "":
                    log_warning(
                        f"Could not convert smiles to 3D SDF with obabel: {predock_cmpd.smiles}"
                    )
                    continue

            predock_cmpd.sdf_3d_path = out_file

        return predock_cmpds
