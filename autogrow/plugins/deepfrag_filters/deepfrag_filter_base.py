"""
DeepFrag plugin.
"""
import __future__
from typing import List
from autogrow.types import Compound
from autogrow.plugins.plugin_base import PluginBase
from abc import abstractmethod
from autogrow.config.argument_vars import ArgumentVars
from scipy.spatial.distance import cosine
import numpy as np
from autogrow.utils.logging import LogLevel, log_info, log_warning
import os
import hashlib
from autogrow.plugins.registry_base import plugin_managers
from autogrow.plugins.deepfrag_filters.common_substructure import find_mcs_and_fragments

DEEPFRAG_DEBUG = False
if DEEPFRAG_DEBUG:
    import os
    import copy

class DeepFragFilterBase(PluginBase):
    """
    Abstract base class for DeepFrag plugins.

    This class defines the interface for DeepFrag plugins and provides a common
    run method that calls the abstract run_deepfrag_filter method.
    """

    apply_on_crossover = False
    fps_fragment_cache = {}
    filter_logger_file = None

    @property
    def plugin_type_name(self) -> str:
        """Return the user-friendly name of the plugin type."""
        return "DeepFrag Filter"

    @property
    def plugin_description(self) -> str:
        """Return a brief description of the plugin type."""
        return "Filters molecules based on DeepFrag predictions"

    def set_log_file(self, filter_logger_file):
        self.filter_logger_file = filter_logger_file

    def run(self, **kwargs) -> List[Compound]:
        """
        Run the DeepFrag plugin.

        Args:
            predocked_cmpds (List[Compound]): A list of Compound objects after filtering
            out molecules that are not suitable for the docking according to the DeepFrag
            filter.
            cutoff: the cosine similarity value to be considered as cutoff.

        Returns:
            List[Compound]: A list of Compound objects that passed the filter,
                with their `sort_score` attribute populated.
        """
        compounds = kwargs["compounds"]
        cutoff = kwargs["input_params"]["deepfrag_cosine_similarity_cutoff"]
        receptor = kwargs["input_params"]["receptor_path"]
        generation_num = kwargs["input_params"].get("generation_num", "NA")
        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        log_info(
            f"Applying DeepFrag filter with cutoff {cutoff} to {len(compounds)} compounds."
        )
        with LogLevel():
            final_compound_list: List[Compound] = []
            for compound in compounds:
                log_info(
                    # compound.id is temporary, so no need to mention it
                    # f"Processing compound {compound.id} with smiles string {compound.smiles}"
                    f"Processing compound with smiles string {compound.smiles}"
                )
                with LogLevel():
                    passed_filter = False
                    similarity_str = "None"
                    child_mol = chemtoolkit.mol_from_smiles(compound.smiles)

                    if len(compound.parent_3D_mols) == 1:
                        # This is a standard mutation
                        parent_mol = compound.parent_3D_mols[0]
                        mcs_mol, frag_mol, fragments, mcs_smarts, connection_points_3d = find_mcs_and_fragments(
                            parent_mol,
                            child_mol,
                            compound_id=compound.id,
                        )

                        if fragments is not None and len(fragments) > 1:
                            log_warning(
                                f"Compound {compound.id} ({compound.smiles}) has {len(fragments)} R-groups attached to the MCS."
                            )

                        similarity = self.__compute_cosine_similarity(
                            receptor, parent_mol, fragments
                        )

                        passed_filter = similarity is not None and similarity >= cutoff
                        if passed_filter:
                            compound.sort_score = similarity
                            final_compound_list.append(compound)

                        similarity_str = (
                            f"{similarity:.3f}" if similarity is not None else "None"
                        )
                        if DEEPFRAG_DEBUG and passed_filter:
                            self._save_debug_image(
                                parent_mol=parent_mol,
                                child_mol=child_mol,
                                mcs_smarts=mcs_smarts,
                                mcs_mol_with_coords=mcs_mol,
                                fragments_mol=frag_mol,
                                connection_points_3d=connection_points_3d,
                                similarity=similarity,
                                passed_filter=passed_filter,
                                compound_id=compound.id,
                                generation_num=generation_num,
                            )
                    elif self.apply_on_crossover and len(compound.parent_3D_mols) == 2:
                        # This is a standard crossover
                        parent_mol_1 = compound.parent_3D_mols[0]
                        mcs_mol_1, frag_mol_1, fragments_1, mcs_smarts_1, connection_points_3d_1 = find_mcs_and_fragments(
                            parent_mol_1,
                            child_mol,
                            compound_id=f"{compound.id}_parent1",
                        )
                        if fragments_1 is not None and len(fragments_1) > 1:
                            log_warning(
                                f"Compound {compound.id} ({compound.smiles}) with parent 1 has {len(fragments_1)} R-groups attached to the MCS."
                            )
                        similarity_0 = self.__compute_cosine_similarity(
                            receptor, parent_mol_1, fragments_1
                        )

                        parent_mol_2 = compound.parent_3D_mols[1]
                        mcs_mol_2, frag_mol_2, fragments_2, mcs_smarts_2, connection_points_3d_2 = find_mcs_and_fragments(
                            parent_mol_2,
                            child_mol,
                            compound_id=f"{compound.id}_parent2",
                        )
                        if fragments_2 is not None and len(fragments_2) > 1:
                            log_warning(
                                f"Compound {compound.id} ({compound.smiles}) with parent 2 has {len(fragments_2)} R-groups attached to the MCS."
                            )
                        similarity_1 = self.__compute_cosine_similarity(
                            receptor, parent_mol_2, fragments_2
                        )

                        sim_0_passes = similarity_0 is not None and similarity_0 >= cutoff
                        sim_1_passes = similarity_1 is not None and similarity_1 >= cutoff
                        passed_filter = sim_0_passes or sim_1_passes

                        if passed_filter:
                            max_similarity = -1.0
                            if sim_0_passes and sim_1_passes:
                                max_similarity = max(similarity_0, similarity_1)
                            elif sim_0_passes:
                                max_similarity = similarity_0
                            else:
                                max_similarity = similarity_1
                            compound.sort_score = max_similarity
                            final_compound_list.append(compound)
                        similarity_str = (
                            f"{similarity_0:.3f} and {similarity_1:.3f}"
                            if similarity_0 is not None and similarity_1 is not None
                            else "None"
                        )
                        if DEEPFRAG_DEBUG and passed_filter:
                            self._save_debug_image(
                                parent_mol=parent_mol_1,
                                child_mol=child_mol,
                                mcs_smarts=mcs_smarts_1,
                                mcs_mol_with_coords=mcs_mol_1,
                                fragments_mol=frag_mol_1,
                                connection_points_3d=connection_points_3d_1,
                                similarity=similarity_0,
                                passed_filter=passed_filter,
                                compound_id=f"{compound.id}_parent1",
                                generation_num=generation_num,
                            )
                            self._save_debug_image(
                                parent_mol=parent_mol_2,
                                child_mol=child_mol,
                                mcs_smarts=mcs_smarts_2,
                                mcs_mol_with_coords=mcs_mol_2,
                                fragments_mol=frag_mol_2,
                                connection_points_3d=connection_points_3d_2,
                                similarity=similarity_1,
                                passed_filter=passed_filter,
                                compound_id=f"{compound.id}_parent2",
                                generation_num=generation_num,
                            )

                    compound.mol_3D = None
                    compound.parent_3D_mols = None

                    # Log the outcome
                    if passed_filter:
                        log_info(
                            f"Candidate molecule {compound.id} with smiles string {compound.smiles} passed the similarity criterion using DeepFrag: {similarity_str}"
                        )
                    else:
                        mesg = f"Candidate molecule {compound.id} with smiles string {compound.smiles} did not fulfill with the similarity criterion using DeepFrag: {similarity_str}"
                        log_warning(mesg)
                        if self.filter_logger_file:
                            self.filter_logger_file.info(mesg)
            return final_compound_list

    @abstractmethod
    def get_prediction_for_parent_receptor(self, parent_mol, receptor, branching_point):
        pass

    @abstractmethod
    def get_fingerprints_for_fragment(self, fragment):
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        self.apply_on_crossover = bool(params["deepfrag_filter_for_crossover"])

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the
        DeepFrag filter.
        It allows users to enable the filter via command-line options.

        Returns:
            List[ArgumentVars]: A list of ArgumentVars objects defining the
            command-line arguments.
        """
        return [
            ArgumentVars(
                name="deepfrag_cosine_similarity_cutoff",
                type=float,
                default=0.35,
                help="The minimum cosine similarity score required for a molecule to pass any enabled DeepFrag filter. Default: 0.35",
            ),
            ArgumentVars(
                name="deepfrag_filter_for_crossover",
                action="store_true",
                default=False,
                help="Apply a DeepFrag filter on the new compounds obtained by crossover.",
            ),
        ]

    def _save_debug_image(self, parent_mol, child_mol, mcs_smarts, mcs_mol_with_coords, fragments_mol, connection_points_3d, similarity, passed_filter, compound_id, generation_num):
        """
        Saves a debug image showing the parent, all fragments, and the chosen fragment.

        Args:
            parent_mol (Any): The RDKit molecule object of the parent.
            child_mol (Any): The RDKit molecule object of the child.
            mcs_smarts (str): The SMARTS string of the MCS.
            mcs_mol_with_coords (Any): The RDKit molecule of the MCS with parent coordinates.
            fragments_mol (Any): The RDKit molecule containing the fragment(s).
            connection_points_3d (dict): A dictionary of connection points.
            similarity (float): The calculated cosine similarity.
            passed_filter (bool): Whether the compound passed the filter.
            compound_id (str): The ID of the compound.
            generation_num (int or str): The generation number.
        """
        os.makedirs("./deepfrag_debug", exist_ok=True)
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        # Create a unique hash for the image based on its contents
        parent_smiles = chemtoolkit.mol_to_smiles(parent_mol)
        child_smiles = chemtoolkit.mol_to_smiles(child_mol)
        fragment_smiles = chemtoolkit.mol_to_smiles(fragments_mol)
        branching_points_str = str(sorted(connection_points_3d.keys()))

        unique_str = f"{parent_smiles}|{child_smiles}|{fragment_smiles}|{branching_points_str}|{compound_id}"
        file_hash = hashlib.md5(unique_str.encode()).hexdigest()

        filename_prefix = "_PASSED_" if passed_filter else ""
        filename = f"./deepfrag_debug/{filename_prefix}gen{generation_num}_{file_hash}.png"

        if os.path.exists(filename):
            log_info(f"Debug image already exists, skipping: {filename}")
            return
        # Create a 2D representation of the parent molecule for better visualization.
        parent_mol_2d = copy.deepcopy(parent_mol)
        chemtoolkit.compute_2d_coords(parent_mol_2d)
        branching_points_indices = list(connection_points_3d.keys())
        mols_to_draw = [parent_mol_2d, child_mol, chemtoolkit.mol_from_smarts(mcs_smarts), mcs_mol_with_coords, fragments_mol]
        similarity_text = f"Similarity: {similarity:.3f}" if similarity is not None else "Similarity: None"
        legends = ["3D Parent, 2D Vis", "Child", f"MCS ({mcs_smarts})", "MCS with Parent Coords", f"Fragment(s)\n{similarity_text}"]
        highlight_list = [branching_points_indices] + [[] for _ in range(len(mols_to_draw) - 1)]
        img = chemtoolkit.mols_to_grid_image(
            mols_to_draw,
            mols_per_row=5,
            sub_img_size=(300, 300),
            highlight_atom_lists=highlight_list,
            legends=legends,
        )
        img.save(filename)
        log_info(f"Saved DeepFrag MCS debug image to {filename}")

    def __compute_cosine_similarity(self, receptor, parent_mol, fragments):
        """
        Compute the cosine similarity between the fingerprints of the receptor-parent complex and the fragment.

        Args:
            receptor: .pdb file containing the receptor.
            parent_mol: RDKit molecule representing the parent interacting with the receptor.
            fragments: list of dictionaries, where each dictionary contains information of a chemical fragment.

        Returns:
            The cosine similarity value.
        """
        similarity = 0.0
        if not fragments:
            return 0.0

        for fragment_info in fragments:
            fragment_smiles = fragment_info['fragment_smiles']
            branching_point = fragment_info['coordinates']

            # No cache because this calculation depends on the receptor and a specific branching point,
            # and it is strange that a branching point can be used twice in the same or different molecules
            fps_receptor_parent_np = self.get_prediction_for_parent_receptor(parent_mol, receptor, branching_point)
            fps_receptor_parent = fps_receptor_parent_np.tolist()

            fps_fragment_from_cache = self.fps_fragment_cache.get(fragment_smiles)
            if fps_fragment_from_cache is None:
                chemtoolkit = plugin_managers.ChemToolkit.toolkit
                fragment_mol = chemtoolkit.mol_from_smiles(fragment_smiles)
                if fragment_mol is None:
                    log_warning(f"Could not create molecule from fragment SMILES: {fragment_smiles}. Treating as zero vector.")
                    fps_fragment_np = np.zeros(2048)
                else:
                    fps_fragment_np = self.get_fingerprints_for_fragment(fragment_mol)

                fps_fragment = fps_fragment_np.tolist()
                self.fps_fragment_cache[fragment_smiles] = fps_fragment
            else:
                fps_fragment_np = np.array(fps_fragment_from_cache)
                fps_fragment = fps_fragment_from_cache

            # Check for zero vectors to avoid nan from cosine distance
            if np.all(fps_receptor_parent_np == 0) or np.all(fps_fragment_np == 0):
                current_similarity = 0.0
            else:
                cosine_distance = cosine(fps_receptor_parent, fps_fragment)
                if np.isnan(cosine_distance):
                    log_warning(f"Cosine similarity resulted in NaN for fragment {fragment_smiles}. Treating as 0 similarity.")
                    current_similarity = 0.0
                else:
                    current_similarity = 1.0 - cosine_distance

            similarity += current_similarity

        similarity = (similarity / len(fragments)) if len(fragments) > 0 else 1.0
        return similarity
