"""
DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
fingerprints for receptor-parent pairs.
"""
import __future__

import rdkit
import logging
import numpy as np
from autogrow.config.argument_vars import ArgumentVars
from typing import List, Tuple
from autogrow.plugins.deepfrag_filters.deepfrag_filter_base import DeepFragFilterBase
import os
import wget
import sys

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

try:
    import torch
    import prody
    from io import StringIO
    from collagen.util import rand_rot
    from collagen.core.molecules.mol import Mol
    from apps.deepfrag.model import DeepFragModel
    from collagen.core.voxelization.voxelizer import VoxelParamsDefault

    numba_logger = logging.getLogger("numba")
    numba_logger.setLevel(logging.WARNING)
    prody.LOGGER._logger.disabled = True
except ImportError as e:
    print("DeepFrag environment (e.g., torch, prody) is not installed. DeepFrag filters will not be available. " + str(e) + "\n")


class DeepFragFilter(DeepFragFilterBase):
    """
    Abstract DeepFrag plugin to be extended depending on the type of fingerprints to be used.
    """

    url_by_in_house_model = {
        "all_best": "https://durrantlab.pitt.edu/apps/deepfrag2/models/all_best.ckpt",
        "gte_4_acid_best": "https://durrantlab.pitt.edu/apps/deepfrag2/models/gte_4_acid_best.ckpt",
        "gte_4_aliphatic_best": "https://durrantlab.pitt.edu/apps/deepfrag2/models/gte_4_aliphatic_best.ckpt",
        "gte_4_aromatic_best": "https://durrantlab.pitt.edu/apps/deepfrag2/models/gte_4_aromatic_best.ckpt",
        "gte_4_base_best": "https://durrantlab.pitt.edu/apps/deepfrag2/models/gte_4_base_best.ckpt",
        "gte_4_best": "https://durrantlab.pitt.edu/apps/deepfrag2/models/gte_4_best.ckpt",
        "lte_3_best": "https://durrantlab.pitt.edu/apps/deepfrag2/models/lte_3_best.ckpt",
    }

    recep = None
    model = None
    ckpt_filename = None
    cpu = False

    def validate(self, params: dict):
        """Validate the provided arguments."""
        super().validate(params)

        self.cpu = bool(params["DeepFragOnCPU"]) or not torch.cuda.is_available()

        if "DeepFragModel" not in params:
            raise Exception("The path of a DeepFrag model should be given as input using"
                            " the '--DeepFragModel' parameter.")

        df_model = params["DeepFragModel"]
        if df_model in self.url_by_in_house_model:
            df_model = self.download_deepfrag_ckpt(df_model + ".ckpt", self.url_by_in_house_model[df_model])
            params["DeepFragModel"] = df_model

        if not os.path.exists(df_model):
            raise Exception(f"DeepFrag model {df_model} is not an in-house model, or does not exist in "
                            f"the path specified.")

        self.ckpt_filename = params["DeepFragModel"]
        self.model = DeepFragModel.load_from_checkpoint(self.ckpt_filename)
        if not self.cpu:
            self.model = self.model.to(torch.device('cuda'))
        self.model.eval()

    def get_prediction_for_parent_receptor(self, parent_mol, receptor, branching_point):
        """
        Calculating DeepFrag fingerprints for the receptor-parent complex.

        Args:
            parent_mol: RDKit molecule representing the parent interacting with the receptor.
            receptor: .pdb file containing the receptor.
            branching_point: coordinates of the branching point.

        Returns:
           Numpy array containing the DeepFrag fingerprints.
        """
        voxel_params = VoxelParamsDefault.DeepFrag

        # Load the receptor
        if self.recep is None:
            with open(receptor, "r") as f:
                m = prody.parsePDBStream(StringIO(f.read()), model=1)
            prody_mol = m.select("all")
            self.recep = Mol.from_prody(prody_mol)

        center = np.array([branching_point.x, branching_point.y, branching_point.z])

        # Load the ligand
        lig = Mol.from_rdkit(parent_mol, strict=False)

        print(f"Using checkpoint {self.ckpt_filename}")

        # You're iterating through multiple checkpoints. This allows output
        # from multiple trained models to be averaged.
        fps = []
        for r in range(8):
            # Random rotations, unless debugging voxels
            rot = rand_rot()

            voxel = self.recep.voxelize(
                voxel_params, cpu=self.cpu, center=center, rot=rot
            )
            voxel = lig.voxelize(
                voxel_params, tensor=voxel, layer_offset=voxel_params.receptor_featurizer.size(), is_receptor=False,
                cpu=self.cpu, center=center, rot=rot
            )

            if not self.cpu:
                voxel = voxel.to(torch.device('cuda'))

            fps.append(self.model.forward(voxel))

        avg_over_ckpts_of_avgs = torch.mean(torch.stack(fps), dim=0)
        return avg_over_ckpts_of_avgs.cpu().detach().numpy()[0]

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        return (
            "DeepFragFilter",
            [
                ArgumentVars(
                    name=self.name,
                    type=float,
                    default=False,
                    help="An value representing the cosine similarity value to be "
                         "considered as cutoff in order to a molecule pass the filter or "
                         "not",
                ),
                ArgumentVars(
                    name="DeepFragModel",
                    type=str,
                    default=None,
                    help="path to the DeepFrag model that is .ckpt file",
                ),
                ArgumentVars(
                    name="DeepFragOnCPU",
                    action="store_true",
                    default=False,
                    help="Use CPU to run DeepFrag.",
                )
            ],
        )

    def download_deepfrag_ckpt(self, deepfrag_model_ckpt, deepfrag_model_url):
        """Download an in-house DeepFrag model checkpoint."""

        current_directory = os.getcwd() + os.sep + "deepfrag_models"
        if not os.path.exists(current_directory):
            os.makedirs(current_directory, exist_ok=True)

        deepfrag_model_path = current_directory + os.sep + deepfrag_model_ckpt
        if not os.path.exists(deepfrag_model_path):
            print(f"Starting download of the DeepFrag model ({deepfrag_model_url} --> {deepfrag_model_path}): {deepfrag_model_ckpt}")
            wget.download(
                deepfrag_model_url,
                deepfrag_model_path,
                self.bar_progress,
            )

        return deepfrag_model_path

    def bar_progress(self, current: float, total: float, width=80):
        """Progress bar for downloading Molbert model.

        Args:
            current (float): Current progress.
            total (float): Total progress.
            width (int, optional): Width of the progress bar. Defaults to 80.
        """
        progress_message = "Downloading DeepFrag model: %d%% [%d / %d] bytes" % (
            current / total * 100,
            current,
            total,
        )
        # Don't use print() as it will print in new line every time.
        sys.stdout.write("\r" + progress_message)
        sys.stdout.flush()
