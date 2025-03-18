"""
DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
fingerprints for receptor-parent pairs.
"""
import __future__

import torch
import rdkit
import rdkit.Chem as Chem
import logging
import numpy as np
from autogrow.types import Compound
from typing import List, Tuple
from autogrow.config.argument_vars import ArgumentVars
from apps.deepfrag.model import DeepFragModel
from collagen.core.voxelization.voxelizer import VoxelParamsDefault
from autogrow.plugins.deepfrag_filters.deepfrag_filter import DeepFragFilterBase

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

try:
    import prody
    from io import StringIO
    from collagen.util import rand_rot
    from collagen.core.molecules.mol import Mol

    numba_logger = logging.getLogger("numba")
    numba_logger.setLevel(logging.WARNING)
    prody.LOGGER._logger.disabled = True
except:
    raise Exception("DeepFrag environment is not installed to be used")


class DeepFragFilterRDKit(DeepFragFilterBase):
    """
    DeepFrag plugin calculating RDKit fingerprints for chemical fragments and DeepFrag
    fingerprints for receptor-parent pairs.
    """

    recep = None
    ckpt_filename = None

    def run(self, **kwargs) -> List[Compound]:
        if "DeepFragModel" not in kwargs["input_params"]:
            raise Exception("The path of a DeepFrag model should be given as input using the '--DeepFragModel' parameter.")

        self.ckpt_filename = kwargs["input_params"]["DeepFragModel"]
        return super().run(**kwargs)

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
        device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")

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
        model = DeepFragModel.load_from_checkpoint(self.ckpt_filename)
        model.eval()

        # You're iterating through multiple checkpoints. This allows output
        # from multiple trained models to be averaged.
        fps = []
        for r in range(8):
            # Random rotations, unless debugging voxels
            rot = rand_rot()

            voxel = self.recep.voxelize(
                voxel_params, cpu=device, center=center, rot=rot
            )
            voxel = lig.voxelize(
                voxel_params, tensor=voxel, layer_offset=voxel_params.receptor_featurizer.size(), is_receptor=False,
                cpu=device, center=center, rot=rot
            )

            fps.append(model.forward(voxel))

        avg_over_ckpts_of_avgs = torch.mean(torch.stack(fps), dim=0)
        return avg_over_ckpts_of_avgs.cpu().detach().numpy()[0]

    def get_fingerprints_for_fragment(self, fragment):
        """
        Calculating a RDKit binary fingerprints for a chemical fragment.

        Args:
            fragment: RDKit molecule representing a chemical fragment.

        Returns:
           Numpy array containing the binary fingerprints calculated with RDKit library.
        """
        try:
            fp = Chem.rdmolops.RDKFingerprint(fragment, maxPath=10, fpSize=2048)
            n_fp = list(map(int, list(fp.ToBitString())))
            return np.array(n_fp)
        except:
            return np.zeros(2048)

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
                )
            ],
        )
