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

DEEPFRAG_DEBUG = True
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

        with LogLevel():
            log_info(
                f"Applying DeepFrag filter with cutoff {cutoff} to {len(compounds)} compounds."
            )
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
                        mcs_mol, frag_mol, fragments, mcs_smarts, connection_points_3d = self.__find_mcs_and_fragments(
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
                        mcs_mol_1, frag_mol_1, fragments_1, mcs_smarts_1, connection_points_3d_1 = self.__find_mcs_and_fragments(
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
                        mcs_mol_2, frag_mol_2, fragments_2, mcs_smarts_2, connection_points_3d_2 = self.__find_mcs_and_fragments(
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
                            f"Docked molecule {compound.id} with smiles string {compound.smiles} passed the similarity criterion using DeepFrag: {similarity_str}"
                        )
                    else:
                        mesg = f"Docked molecule {compound.id} with smiles string {compound.smiles} did not fulfill with the similarity criterion using DeepFrag: {similarity_str}"
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

    def _find_bridge_atoms(self, mol):
        """
        Identifies bridge atoms in a molecule.

        A bridge atom is an atom whose removal would increase the number of
        connected components (fragments) in the molecule.

        Args:
            mol (Any): An RDKit molecule object.

        Returns:
            List[int]: A list of indices of the bridge atoms.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        bridges = []
        num_frags_initial = len(chemtoolkit.get_mol_frags(mol))
        for atom in chemtoolkit.get_atoms(mol):
            emol = chemtoolkit.get_editable_mol(mol)
            chemtoolkit.remove_atom_from_editable_mol(emol, chemtoolkit.get_idx(atom))
            frag_mol = chemtoolkit.get_noneditable_mol(emol)
            num_frags_after_removal = len(chemtoolkit.get_mol_frags(frag_mol))
            if num_frags_after_removal > num_frags_initial:
                bridges.append(chemtoolkit.get_idx(atom))
        return bridges

    def _find_single_fragment_mcs(self, mol1, mol2):
        """
        Finds the largest non-fragmenting Maximum Common Substructure (MCS).

        This function identifies the largest common substructure between two
        molecules that, when removed, does not break either molecule into
        multiple chemically meaningful pieces. A "meaningful" piece is defined
        as a fragment containing at least one heavy (non-hydrogen) atom. This is
        crucial for identifying a single, coherent fragment for analysis by
        DeepFrag.

        The process can be computationally intensive and uses a two-phase approach
        for efficiency:
        1.  **Initial MCS & Pruning**: It first finds the absolute largest MCS. If
            this MCS causes fragmentation, it performs a "smart pruning" step. It
            identifies all "bridge atoms" within the MCS (atoms that connect
            different parts of the parent molecules) and removes them to create a
            smaller, more stable MCS candidate.
        2.  **Iterative Search**: This new, smaller candidate becomes the starting
            point for a more exhaustive iterative search, which removes one atom
            at a time until a valid, non-fragmenting MCS is found. This two-step
            process significantly reduces the search space and speeds up the
            calculation for complex molecules.

        Args:
            mol1 (Chem.Mol): The first RDKit molecule.
            mol2 (Chem.Mol): The second RDKit molecule.

        Returns:
            str: The SMARTS string of the optimal MCS, or None if no suitable MCS is found.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        if mol1 is None or mol2 is None:
            raise ValueError("Invalid RDKit molecule provided.")

        # Helper function to get fragments by removing MCS atoms
        def _get_fragments_for_mol(mol, mcs_mol):
            match = chemtoolkit.get_substruct_match(mol, mcs_mol)
            if not match:
                return [chemtoolkit.mol_to_smiles(mol)]  # Should not happen with valid MCS
            emol = chemtoolkit.get_editable_mol(mol)
            atoms_to_remove = sorted(list(match), reverse=True)
            for atom_idx in atoms_to_remove:
                chemtoolkit.remove_atom_from_editable_mol(emol, atom_idx)
            frag_mol = chemtoolkit.get_noneditable_mol(emol)
            try:
                # Important: Sanitize=False to handle radical fragments, then sanitize later if needed.
                frag_smiles = chemtoolkit.mol_to_smiles(frag_mol, isomeric_smiles=True)
                # Split into fragments and remove empty strings from result
                raw_frags = [s for s in frag_smiles.split('.') if s]
                meaningful_frags = []
                for frag_smi in raw_frags:
                    frag_mol_obj = chemtoolkit.mol_from_smiles(frag_smi, sanitize=False)
                    if frag_mol_obj is not None and chemtoolkit.lipinski_heavy_atom_count(frag_mol_obj) > 0:
                        meaningful_frags.append(frag_smi)
                return meaningful_frags
            except:
                return []

        def _get_atoms_sorted_by_connectivity(mol):
            """
            Returns atoms sorted by connectivity (degree), with terminal atoms first.
            
            Args:
                mol: RDKit molecule object
                
            Returns:
                list: Atoms sorted by degree (ascending), then by atom index for consistency
            """
            atoms_with_degree = []
            for atom in chemtoolkit.get_atoms(mol):
                degree = chemtoolkit.get_atom_degree(atom)
                atoms_with_degree.append((degree, chemtoolkit.get_idx(atom), atom))
            # Sort by degree (ascending, so terminal atoms first), then by index for consistency
            atoms_with_degree.sort(key=lambda x: (x[0], x[1]))
            
            return [atom for _, _, atom in atoms_with_degree]

        # 1. Find the initial, absolute MCS as a starting point.
        initial_mcs = chemtoolkit.find_mcs([mol1, mol2], timeout=30)
        if initial_mcs.numAtoms == 0:
            return None

        # Our list of candidates to check, starting with the biggest.
        # We store tuples of (mol_object, smarts_string).
        mcs_mol = chemtoolkit.mol_from_smarts(initial_mcs.smartsString)
        candidates = [(mcs_mol, initial_mcs.smartsString)]

        # Keep track of SMARTS we've already processed to avoid redundant work.
        seen_smarts = {initial_mcs.smartsString}

        # Check the initial MCS first
        initial_frags1 = _get_fragments_for_mol(mol1, mcs_mol)
        initial_frags2 = _get_fragments_for_mol(mol2, mcs_mol)
        if len(initial_frags1) <= 1 and len(initial_frags2) <= 1:
            return initial_mcs.smartsString

        # Initial MCS is fragmenting, so we start the optimization
        log_warning(
            "Initial MCS causes fragmentation. Starting optimized search for a "
            "non-fragmenting MCS. This may be slow."
        )

        # OPTIMIZATION: Prune by removing bridge atoms from the MCS
        bridge_atoms_in_mol1 = self._find_bridge_atoms(mol1)
        bridge_atoms_in_mol2 = self._find_bridge_atoms(mol2)

        parent1_match = chemtoolkit.get_substruct_match(mol1, mcs_mol)
        parent2_match = chemtoolkit.get_substruct_match(mol2, mcs_mol)

        problem_mcs_indices = set()
        if parent1_match:
            for i, parent_idx in enumerate(parent1_match):
                if parent_idx in bridge_atoms_in_mol1:
                    problem_mcs_indices.add(i)

        if parent2_match:
            for i, parent_idx in enumerate(parent2_match):
                if parent_idx in bridge_atoms_in_mol2:
                    problem_mcs_indices.add(i)

        if problem_mcs_indices:
            emol = chemtoolkit.get_editable_mol(mcs_mol)
            for idx in sorted(list(problem_mcs_indices), reverse=True):
                chemtoolkit.remove_atom_from_editable_mol(emol, idx)

            pruned_mcs_mol_base = chemtoolkit.get_noneditable_mol(emol)
            sub_frags = chemtoolkit.get_mol_frags(pruned_mcs_mol_base, as_mols=True)
            if sub_frags:
                largest_sub_frag = max(sub_frags, key=lambda m: chemtoolkit.get_num_atoms(m))
                try:
                    pruned_smarts = chemtoolkit.mol_to_smarts(largest_sub_frag)
                    if pruned_smarts not in seen_smarts:
                        candidates.append((largest_sub_frag, pruned_smarts))
                        seen_smarts.add(pruned_smarts)
                except:
                    pass  # Ignore if SMARTS generation fails

        # Now, proceed with the iterative search, which will start with better candidates
        while candidates:
            # 2. Always check the largest candidate first.
            candidates.sort(key=lambda x: chemtoolkit.get_num_atoms(x[0]), reverse=True)
            current_mcs_mol, current_mcs_smarts = candidates.pop(0)

            # 3. Check if this candidate is "valid".
            frags1 = _get_fragments_for_mol(mol1, current_mcs_mol)
            frags2 = _get_fragments_for_mol(mol2, current_mcs_mol)

            if len(frags1) <= 1 and len(frags2) <= 1:
                # SUCCESS: We found the largest MCS that satisfies the condition.
                return current_mcs_smarts

            # 4. If not valid, generate smaller candidates by removing atoms in connectivity order.
            if chemtoolkit.get_num_atoms(current_mcs_mol) > 1:
                # Get atoms sorted by connectivity (terminal atoms first)
                sorted_atoms = _get_atoms_sorted_by_connectivity(current_mcs_mol)
                
                for atom in sorted_atoms:
                    # Create a copy to modify
                    emol = chemtoolkit.get_editable_mol(current_mcs_mol)
                    chemtoolkit.remove_atom_from_editable_mol(emol, chemtoolkit.get_idx(atom))
                    # Removing an atom can disconnect the MCS itself. We take the largest resulting piece.
                    sub_frags = chemtoolkit.get_mol_frags(chemtoolkit.get_noneditable_mol(emol), as_mols=True)
                    if not sub_frags:
                        continue
                    largest_sub_frag = max(sub_frags, key=lambda m: chemtoolkit.get_num_atoms(m))
                    try:
                        new_smarts = chemtoolkit.mol_to_smarts(largest_sub_frag)
                        if new_smarts not in seen_smarts:
                            # Add new, smaller candidate to our list for checking.
                            candidates.append((largest_sub_frag, new_smarts))
                            seen_smarts.add(new_smarts)
                    except:
                        continue
        
        # If the loop finishes, no suitable MCS was found.
        return None

    # Create a new MCS molecule with 3D coordinates from parent
    def __create_mcs_molecule(self, parent, child):
        """
        Create a new molecule representing just the MCS with 3D coordinates from parent.

        This function orchestrates the creation of a 3D-aware MCS molecule. It
        first calls `_find_single_fragment_mcs` to determine the appropriate
        MCS SMARTS pattern, then constructs a new RDKit molecule with 3D
        coordinates inherited from the parent molecule.

        Returns:
            - mcs_mol: The MCS molecule with 3D coordinates
            - parent_to_mcs_map: Mapping from parent atom indices to MCS atom indices
            - child_to_mcs_map: Mapping from child atom indices to MCS atom indices
            - mcs_smarts: The SMARTS string of the MCS.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        try:
            parent = chemtoolkit.remove_hs(parent)
            child = chemtoolkit.remove_hs(child)
        except Exception as e:
            log_warning(
                f"Failed to remove hydrogens from parent or child molecule, skipping MCS. Error: {e}"
            )
            return None, {}, {}, ""
        # Find the Maximum Common Substructure that does not fragment molecules
        mcs_smarts = self._find_single_fragment_mcs(parent, child)

        if mcs_smarts is None:
            log_warning(f"Could not find single-fragment MCS for {chemtoolkit.mol_to_smiles(parent)} and {chemtoolkit.mol_to_smiles(child)}")
            return None, {}, {}, ""

        # Create a molecule from the MCS SMARTS pattern
        mcs_mol = chemtoolkit.mol_from_smarts(mcs_smarts)
        log_info(f"MCS found: {mcs_smarts}")
        log_info(f"Number of atoms in MCS: {chemtoolkit.get_num_atoms(mcs_mol)}")
        log_info(f"Number of bonds in MCS: {len(mcs_mol.GetBonds())}")
        # Match the MCS in both molecules
        parent_match = chemtoolkit.get_substruct_match(parent, mcs_mol)
        child_match = chemtoolkit.get_substruct_match(child, mcs_mol)
        if not parent_match or not child_match:
            log_warning("No match found in one or both molecules.")
            return None, {}, {}, ""

        # Create a new editable molecule for the MCS
        mcs_editable = chemtoolkit.create_empty_editable_mol()
        # Add atoms to the new molecule
        atom_mapping = {}  # Maps MCS query atom idx to the new molecule atom idx
        for i, atom in enumerate(chemtoolkit.get_atoms(mcs_mol)):
            # Get the corresponding atom from parent based on the match
            parent_atom_idx = parent_match[i]
            parent_atom = chemtoolkit.get_atom_with_idx(parent, parent_atom_idx)
            # Create a new atom with the same properties
            new_atom = chemtoolkit.copy_atom(parent_atom)
            # Add the atom to the new molecule
            new_idx = chemtoolkit.add_atom_to_mol(mcs_editable, new_atom)
            atom_mapping[i] = new_idx

        # Add bonds to the new molecule
        for bond in mcs_mol.GetBonds():
            begin_atom_idx = atom_mapping[bond.GetBeginAtomIdx()]
            end_atom_idx = atom_mapping[bond.GetEndAtomIdx()]
            bond_type = chemtoolkit.get_single_bond_type()  # Default to single bond
            # Get the corresponding bond from parent
            parent_begin = parent_match[bond.GetBeginAtomIdx()]
            parent_end = parent_match[bond.GetEndAtomIdx()]
            parent_bond = chemtoolkit.get_bond_between_atoms(parent, parent_begin, parent_end)
            if parent_bond:
                bond_type = chemtoolkit.get_bond_type(parent_bond)
            chemtoolkit.add_bond_to_mol(mcs_editable, begin_atom_idx, end_atom_idx, bond_type)
        # Create the final MCS molecule
        new_mcs_mol = chemtoolkit.get_noneditable_mol(mcs_editable)
        # Add a conformer to the new molecule to store 3D coordinates
        conf = chemtoolkit.create_conformer(chemtoolkit.get_num_atoms(new_mcs_mol))
        # Copy coordinates from parent to the new MCS molecule
        parent_conformer = chemtoolkit.get_conformer(parent)
        for i, idx in enumerate(parent_match):
            new_idx = atom_mapping[i]
            pos = chemtoolkit.get_atom_position(parent_conformer, idx)
            chemtoolkit.set_atom_position_in_conformer(conf, new_idx, pos)
        chemtoolkit.add_conformer_to_mol(new_mcs_mol, conf)
        # Try to sanitize the MCS molecule
        try:
            chemtoolkit.sanitize_mol(new_mcs_mol)
        except:
            log_warning("MCS molecule sanitization failed")

        # Add MCS to parent atom mapping
        mcs_to_parent_map = {
            atom_mapping[i]: parent_match[i] for i in range(len(parent_match))
        }
        # Add MCS to child atom mapping
        mcs_to_child_map = {
            atom_mapping[i]: child_match[i] for i in range(len(child_match))
        }
        return new_mcs_mol, mcs_to_parent_map, mcs_to_child_map, mcs_smarts

    # Function to find MCS and remove it from the second molecule
    def __find_mcs_and_fragments(self, parent, child, compound_id=None):
        """
        Find the MCS, identify fragments, and determine 3D connection points.

        This is the main orchestrator for preparing a molecule for DeepFrag
        analysis. It performs three key steps:
        1. Calls `__create_mcs_molecule` to obtain a 3D MCS molecule and atom
           mappings between the parent, child, and MCS.
        2. Identifies the fragments of the child molecule that are not part of
           the MCS.
        3. Determines the 3D coordinates of the "connection points" on the
           parent molecule where these fragments would attach.

        Args:
            parent (Any): The parent RDKit molecule.
            child (Any): The child RDKit molecule.
            compound_id (Optional[str], optional): The ID of the compound, used
                for debugging. Defaults to None.

        Returns:
            tuple: A tuple containing the MCS molecule, the fragment molecule, a
                list of fragment information dictionaries, the MCS SMARTS string,
                and a dictionary of 3D connection points.
        """
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        # Create an explicit MCS molecule with 3D coordinates from parent
        mcs_mol, mcs_to_parent_map, mcs_to_child_map, mcs_smarts = self.__create_mcs_molecule(
            parent, child
        )
        if mcs_mol is None:
            log_warning("Failed to create MCS molecule.")
            return None, None, [], None, {}

        # Create reverse mappings
        child_to_mcs_map = {v: k for k, v in mcs_to_child_map.items()}

        # Create an RWMol for editing
        rwmol = chemtoolkit.get_editable_mol(child)
        # Find atoms connected to the MCS but not part of it (meaning, attachment points)
        attachment_bonds = []
        connection_points_3d = {}  # Store connection info
        assert chemtoolkit.get_num_conformers(parent) == 1
        parent_conf = chemtoolkit.get_conformer(parent)
        # Find attachment points using the MCS molecule and mappings
        for mcs_atom_idx, child_atom_idx in mcs_to_child_map.items():
            child_atom = chemtoolkit.get_atom_with_idx(child, child_atom_idx)
            for neighbor in chemtoolkit.get_neighbors(child_atom):
                neighbor_idx = chemtoolkit.get_idx(neighbor)
                if neighbor_idx not in child_to_mcs_map:
                    # This is a bond to an attachment point
                    # Find the corresponding atom in the parent using the MCS mapping
                    parent_atom_idx = mcs_to_parent_map.get(mcs_atom_idx)

                    if parent_atom_idx is not None:
                        # Store only the necessary connection details
                        connection_points_3d[parent_atom_idx] = {
                            "mcs_atom_idx": mcs_atom_idx,
                            "neighbor_symbol": chemtoolkit.get_atom_with_idx(
                                child, neighbor_idx
                            ).GetSymbol(),
                            "coordinates": chemtoolkit.get_atom_position(parent_conf, parent_atom_idx),
                            # Only needed for matching later
                        }

                    attachment_bonds.append((child_atom_idx, neighbor_idx))

        # Create dummy atoms at attachment points
        dummy_atoms = []
        for child_atom_idx, neighbor_idx in attachment_bonds:
            # Add a dummy atom (R group)
            dummy_idx = chemtoolkit.add_atom_to_mol(rwmol, chemtoolkit.create_atom(0))
            # Add a bond from the dummy atom to the neighbor atom
            chemtoolkit.add_bond_to_mol(rwmol, dummy_idx, neighbor_idx, chemtoolkit.get_single_bond_type())
            dummy_atoms.append(
                (dummy_idx, child_atom_idx)
            )  # Store the dummy atom and its corresponding MCS atom
            # Find the parent atom index that corresponds to this child atom index
            for parent_idx, connection_info in connection_points_3d.items():
                if connection_info.get("mcs_atom_idx") == child_to_mcs_map.get(
                    child_atom_idx
                ):
                    connection_info["dummy_idx"] = dummy_idx
                    atom = chemtoolkit.get_atom_with_idx(rwmol, dummy_idx)
                    chemtoolkit.set_atom_property(
                        atom, "atom_connecting_mcs", str(connection_info.get("mcs_atom_idx"))
                    )
        # Convert the child_to_mcs_map keys to a set for faster lookups
        mcs_atoms_in_child = set(child_to_mcs_map.keys())

        # Now delete atoms in the MCS (in reverse order to avoid index shifting issues)
        atoms_to_delete = sorted(list(mcs_atoms_in_child), reverse=True)
        for atom_idx in atoms_to_delete:
            chemtoolkit.remove_atom_from_editable_mol(rwmol, atom_idx)
        # Convert to a molecule
        frag_mol = chemtoolkit.get_noneditable_mol(rwmol)
        # Get fragments (in case there are disconnected fragments)
        # First get the atom mapping for each fragment
        atom_frags = chemtoolkit.get_mol_frags(frag_mol)
        # Then get the actual fragment molecules
        frags = [
            chemtoolkit.get_mol_frags(frag_mol, as_mols=True, sanitize_frags=False)[i]
            for i in range(len(atom_frags))
        ]
        # For each fragment, identify which dummy atom it contains
        fragment_info = []
        for i, frag in enumerate(frags):
            frag_smiles = chemtoolkit.mol_to_smiles(frag)
            # We'll track which connections this fragment has
            fragment_connections = []

            # Look for dummy atoms in this fragment
            for atom in chemtoolkit.get_atoms(frag):
                if atom.GetSymbol() == "*":
                    # For each dummy atom, identify its connections
                    # dummy_neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]

                    # Try to match this fragment to a connection point
                    matched_connection = False

                    for cp_info in connection_points_3d.values():
                        # Generic matching based on atom indices and connectivity
                        # Use the dummy atom's neighbor to match with the appropriate connection point
                        if str(cp_info["mcs_atom_idx"]) in chemtoolkit.get_atom_property(
                            atom, "atom_connecting_mcs"
                        ):
                            connection_info = {
                                "fragment_smiles": frag_smiles,
                                "mcs_atom_idx": cp_info["mcs_atom_idx"],
                                "coordinates": cp_info["coordinates"],
                            }

                            fragment_connections.append(connection_info)
                            matched_connection = True
                            break  # Take only the first match

                    if not matched_connection:
                        log_warning(
                            f"Could not match fragment {frag_smiles} to a connection point"
                        )
            # Add all connection points for this fragment
            fragment_info.extend(fragment_connections)

        # Sanitize to fix any valence issues
        try:
            chemtoolkit.sanitize_mol(frag_mol)
        except:
            log_warning("Sanitization failed, the fragment may have valence issues")
            # Try to get the molecule anyway
            chemtoolkit.rdmolops_get_sssr(frag_mol)
        # Return without connection_points_3d
        return mcs_mol, frag_mol, fragment_info, mcs_smarts, connection_points_3d

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
