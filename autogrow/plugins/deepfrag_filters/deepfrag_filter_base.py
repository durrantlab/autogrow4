"""
DeepFrag plugin.
"""
import __future__
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import rdFMCS
from typing import List, Tuple
from autogrow.types import Compound
from autogrow.plugins.plugin_base import PluginBase
from abc import abstractmethod
from autogrow.config.argument_vars import ArgumentVars
from scipy.spatial.distance import cosine
from autogrow.utils.logging import LogLevel, log_debug, log_info, log_warning
import random
import os
import hashlib

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

DEEPFRAG_DEBUG = True
if DEEPFRAG_DEBUG:
    from rdkit.Chem import Draw
    import os


class DeepFragFilterBase(PluginBase):
    """
    Abstract base class for DeepFrag plugins.

    This class defines the interface for DeepFrag plugins and provides a common
    run method that calls the abstract run_deepfrag_filter method.
    """

    apply_on_crossover = False
    fps_fragment_cache = {}
    filter_logger_file = None

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
        cutoff = kwargs["input_params"][self.name]
        receptor = kwargs["input_params"]["receptor_path"]
        generation_num = kwargs["input_params"].get("generation_num", "NA")

        with LogLevel():
            log_info(
                f"Applying DeepFrag filter with cutoff {cutoff} to {len(compounds)} compounds."
            )
            final_compound_list: List[Compound] = []
            for compound in compounds:
                log_info(
                    f"Processing compound {compound.id} with smiles string {compound.smiles}"
                )
                with LogLevel():
                    passed_filter = False
                    similarity_str = "None"
                    child_mol = Chem.MolFromSmiles(compound.smiles)
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
        self.apply_on_crossover = bool(params["DeepFragFilterForCrossover"])

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific for the
        DeepFrag filter.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of ArgumentVars objects defining the
                command-line arguments.
        """
        return (
            "DeepFragFilter",
            [
                ArgumentVars(
                    name=self.name,
                    type=float,
                    default=False,
                    help="The minimum cosine similarity score required for a molecule to pass the filter. \
                        This filter is applied to the products of multi-reactant fragment additions. See \
                        also --DeepFragFilterForCrossover.",
                ),

                ArgumentVars(
                    name="DeepFragFilterForCrossover",
                    action="store_true",
                    default=False,
                    help="Apply a DeepFrag filter on the new compounds obtained by crossover.",
                )
            ],
        )

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
        
        # Create a unique hash for the image based on its contents
        parent_smiles = Chem.MolToSmiles(parent_mol)
        child_smiles = Chem.MolToSmiles(child_mol)
        fragment_smiles = Chem.MolToSmiles(fragments_mol)
        branching_points_str = str(sorted(connection_points_3d.keys()))

        unique_str = f"{parent_smiles}|{child_smiles}|{fragment_smiles}|{branching_points_str}|{compound_id}"
        file_hash = hashlib.md5(unique_str.encode()).hexdigest()

        filename_prefix = "_PASSED_" if passed_filter else ""
        filename = f"./deepfrag_debug/{filename_prefix}gen{generation_num}_{file_hash}.png"

        if os.path.exists(filename):
            log_info(f"Debug image already exists, skipping: {filename}")
            return

        branching_points_indices = list(connection_points_3d.keys())
        
        mols_to_draw = [parent_mol, child_mol, Chem.MolFromSmarts(mcs_smarts), mcs_mol_with_coords, fragments_mol]
        similarity_text = f"Similarity: {similarity:.3f}" if similarity is not None else "Similarity: None"
        legends = ["Parent", "Child", f"MCS ({mcs_smarts})", "MCS with Parent Coords", f"Fragment(s)\n{similarity_text}"]
        
        highlight_list = [branching_points_indices] + [[] for _ in range(len(mols_to_draw) - 1)]
        
        img = Draw.MolsToGridImage(
            mols_to_draw,
            legends=legends,
            molsPerRow=5,
            subImgSize=(300, 300),
            highlightAtomLists=highlight_list
        )
        img.save(filename)
        log_info(f"Saved DeepFrag MCS debug image to {filename}")

    def _find_single_fragment_mcs(self, mol1, mol2):
        """
        Finds the largest common substructure (MCS) between two molecules that does
        not result in more than one fragment for either of the input molecules
        when the MCS is removed.

        Args:
            mol1 (Chem.Mol): The first RDKit molecule.
            mol2 (Chem.Mol): The second RDKit molecule.

        Returns:
            str: The SMARTS string of the optimal MCS, or None if no suitable MCS is found.
        """
        if mol1 is None or mol2 is None:
            raise ValueError("Invalid RDKit molecule provided.")

        # Helper function to get fragments by removing MCS atoms
        def _get_fragments_for_mol(mol, mcs_mol):
            match = mol.GetSubstructMatch(mcs_mol)
            if not match:
                return [Chem.MolToSmiles(mol)] # Should not happen with valid MCS

            emol = Chem.RWMol(mol)
            atoms_to_remove = sorted(list(match), reverse=True)
            for atom_idx in atoms_to_remove:
                emol.RemoveAtom(atom_idx)
            
            frag_mol = emol.GetMol()
            
            try:
                # Important: Sanitize=False to handle radical fragments, then sanitize later if needed.
                frag_smiles = Chem.MolToSmiles(frag_mol, isomericSmiles=True)
                # Split into fragments and remove empty strings from result
                frags = [s for s in frag_smiles.split('.') if s]
                return frags
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
            for atom in mol.GetAtoms():
                degree = atom.GetDegree()
                atoms_with_degree.append((degree, atom.GetIdx(), atom))
            
            # Sort by degree (ascending, so terminal atoms first), then by index for consistency
            atoms_with_degree.sort(key=lambda x: (x[0], x[1]))
            
            return [atom for _, _, atom in atoms_with_degree]

        # 1. Find the initial, absolute MCS as a starting point.
        initial_mcs = rdFMCS.FindMCS([mol1, mol2], timeout=30)
        if initial_mcs.numAtoms == 0:
            return None

        # Our list of candidates to check, starting with the biggest.
        # We store tuples of (mol_object, smarts_string).
        mcs_mol = Chem.MolFromSmarts(initial_mcs.smartsString)
        candidates = [(mcs_mol, initial_mcs.smartsString)]
        
        # Keep track of SMARTS we've already processed to avoid redundant work.
        seen_smarts = {initial_mcs.smartsString}

        while candidates:
            # 2. Always check the largest candidate first.
            candidates.sort(key=lambda x: x[0].GetNumAtoms(), reverse=True)
            current_mcs_mol, current_mcs_smarts = candidates.pop(0)

            # 3. Check if this candidate is "valid".
            frags1 = _get_fragments_for_mol(mol1, current_mcs_mol)
            frags2 = _get_fragments_for_mol(mol2, current_mcs_mol)

            if len(frags1) <= 1 and len(frags2) <= 1:
                # SUCCESS: We found the largest MCS that satisfies the condition.
                return current_mcs_smarts

            # 4. If not valid, generate smaller candidates by removing atoms in connectivity order.
            if current_mcs_mol.GetNumAtoms() > 1:
                # Get atoms sorted by connectivity (terminal atoms first)
                sorted_atoms = _get_atoms_sorted_by_connectivity(current_mcs_mol)
                
                for atom in sorted_atoms:
                    # Create a copy to modify
                    emol = Chem.RWMol(current_mcs_mol)
                    emol.RemoveAtom(atom.GetIdx())
                    
                    # Removing an atom can disconnect the MCS itself. We take the largest resulting piece.
                    sub_frags = Chem.GetMolFrags(emol.GetMol(), asMols=True)
                    if not sub_frags:
                        continue

                    largest_sub_frag = max(sub_frags, key=lambda m: m.GetNumAtoms())
                    
                    try:
                        new_smarts = Chem.MolToSmarts(largest_sub_frag)
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
        Create a new molecule representing just the MCS with 3D coordinates from parent

        Returns:
            - mcs_mol: The MCS molecule with 3D coordinates
            - parent_to_mcs_map: Mapping from parent atom indices to MCS atom indices
            - child_to_mcs_map: Mapping from child atom indices to MCS atom indices
            - mcs_smarts: The SMARTS string of the MCS.
        """
        parent = Chem.RemoveHs(parent)
        child = Chem.RemoveHs(child)
        
        # Find the Maximum Common Substructure that does not fragment molecules
        mcs_smarts = self._find_single_fragment_mcs(parent, child)

        if mcs_smarts is None:
            log_warning(f"Could not find single-fragment MCS for {Chem.MolToSmiles(parent)} and {Chem.MolToSmiles(child)}")
            return None, {}, {}, ""

        # Create a molecule from the MCS SMARTS pattern
        mcs_mol = Chem.MolFromSmarts(mcs_smarts)
        log_info(f"MCS found: {mcs_smarts}")
        log_info(f"Number of atoms in MCS: {mcs_mol.GetNumAtoms()}")
        log_info(f"Number of bonds in MCS: {mcs_mol.GetNumBonds()}")
        
        # Match the MCS in both molecules
        parent_match = parent.GetSubstructMatch(mcs_mol)
        child_match = child.GetSubstructMatch(mcs_mol)

        if not parent_match or not child_match:
            log_warning("No match found in one or both molecules.")
            return None, {}, {}, ""

        # Create a new editable molecule for the MCS
        mcs_editable = Chem.EditableMol(Chem.Mol())

        # Add atoms to the new molecule
        atom_mapping = {}  # Maps MCS query atom idx to the new molecule atom idx
        for i, atom in enumerate(mcs_mol.GetAtoms()):
            # Get the corresponding atom from parent based on the match
            parent_atom_idx = parent_match[i]
            parent_atom = parent.GetAtomWithIdx(parent_atom_idx)

            # Create a new atom with the same properties
            new_atom = Chem.Atom(parent_atom.GetAtomicNum())
            new_atom.SetFormalCharge(parent_atom.GetFormalCharge())
            new_atom.SetChiralTag(parent_atom.GetChiralTag())
            new_atom.SetHybridization(parent_atom.GetHybridization())
            new_atom.SetNumExplicitHs(parent_atom.GetNumExplicitHs())
            new_atom.SetNoImplicit(parent_atom.GetNoImplicit())
            new_atom.SetIsAromatic(parent_atom.GetIsAromatic())

            # Add the atom to the new molecule
            new_idx = mcs_editable.AddAtom(new_atom)
            atom_mapping[i] = new_idx

        # Add bonds to the new molecule
        for bond in mcs_mol.GetBonds():
            begin_atom = atom_mapping[bond.GetBeginAtomIdx()]
            end_atom = atom_mapping[bond.GetEndAtomIdx()]
            bond_type = Chem.BondType.SINGLE  # Default to single bond

            # Get the corresponding bond from parent
            parent_begin = parent_match[bond.GetBeginAtomIdx()]
            parent_end = parent_match[bond.GetEndAtomIdx()]
            parent_bond = parent.GetBondBetweenAtoms(parent_begin, parent_end)

            if parent_bond:
                bond_type = parent_bond.GetBondType()

            mcs_editable.AddBond(begin_atom, end_atom, bond_type)

        # Create the final MCS molecule
        new_mcs_mol = mcs_editable.GetMol()

        # Add a conformer to the new molecule to store 3D coordinates
        conf = Chem.Conformer(new_mcs_mol.GetNumAtoms())

        # Copy coordinates from parent to the new MCS molecule
        for i, idx in enumerate(parent_match):
            new_idx = atom_mapping[i]
            pos = parent.GetConformer().GetAtomPosition(idx)
            conf.SetAtomPosition(new_idx, pos)

        new_mcs_mol.AddConformer(conf)

        # Try to sanitize the MCS molecule
        try:
            Chem.SanitizeMol(new_mcs_mol)
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
        rwmol = Chem.RWMol(child)

        # Find atoms connected to the MCS but not part of it (meaning, attachment points)
        attachment_bonds = []
        connection_points_3d = {}  # Store connection info
        assert parent.GetNumConformers() == 1
        parent_conf = parent.GetConformer()

        # Find attachment points using the MCS molecule and mappings
        for mcs_atom_idx, child_atom_idx in mcs_to_child_map.items():
            child_atom = child.GetAtomWithIdx(child_atom_idx)

            for neighbor in child_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in child_to_mcs_map:
                    # This is a bond to an attachment point
                    # Find the corresponding atom in the parent using the MCS mapping
                    parent_atom_idx = mcs_to_parent_map.get(mcs_atom_idx)

                    if parent_atom_idx is not None:
                        # Store only the necessary connection details
                        connection_points_3d[parent_atom_idx] = {
                            "mcs_atom_idx": mcs_atom_idx,
                            "neighbor_symbol": child.GetAtomWithIdx(
                                neighbor_idx
                            ).GetSymbol(),
                            "coordinates": parent_conf.GetAtomPosition(parent_atom_idx),
                            # Only needed for matching later
                        }

                    attachment_bonds.append((child_atom_idx, neighbor_idx))

        # Create dummy atoms at attachment points
        dummy_atoms = []
        for child_atom_idx, neighbor_idx in attachment_bonds:
            # Add a dummy atom (R group)
            dummy_idx = rwmol.AddAtom(Chem.Atom("*"))
            # Add a bond from the dummy atom to the neighbor atom
            rwmol.AddBond(dummy_idx, neighbor_idx, Chem.BondType.SINGLE)
            dummy_atoms.append(
                (dummy_idx, child_atom_idx)
            )  # Store the dummy atom and its corresponding MCS atom
            # Find the parent atom index that corresponds to this child atom index
            for parent_idx, connection_info in connection_points_3d.items():
                if connection_info.get("mcs_atom_idx") == child_to_mcs_map.get(
                    child_atom_idx
                ):
                    connection_info["dummy_idx"] = dummy_idx
                    atom = rwmol.GetAtomWithIdx(dummy_idx)
                    atom.SetProp(
                        "atom_connecting_mcs", str(connection_info.get("mcs_atom_idx"))
                    )
        # Convert the child_to_mcs_map keys to a set for faster lookups
        mcs_atoms_in_child = set(child_to_mcs_map.keys())

        # Now delete atoms in the MCS (in reverse order to avoid index shifting issues)
        atoms_to_delete = sorted(list(mcs_atoms_in_child), reverse=True)
        for atom_idx in atoms_to_delete:
            rwmol.RemoveAtom(atom_idx)

        # Convert to a molecule
        frag_mol = rwmol.GetMol()

        # Get fragments (in case there are disconnected fragments)
        # First get the atom mapping for each fragment
        atom_frags = Chem.GetMolFrags(frag_mol)

        # Then get the actual fragment molecules
        frags = [
            Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)[i]
            for i in range(len(atom_frags))
        ]
        # For each fragment, identify which dummy atom it contains
        fragment_info = []
        for i, frag in enumerate(frags):
            frag_smiles = Chem.MolToSmiles(frag)

            # We'll track which connections this fragment has
            fragment_connections = []

            # Look for dummy atoms in this fragment
            for atom in frag.GetAtoms():
                if atom.GetSymbol() == "*":
                    # For each dummy atom, identify its connections
                    # dummy_neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]

                    # Try to match this fragment to a connection point
                    matched_connection = False

                    for cp_info in connection_points_3d.values():
                        # Generic matching based on atom indices and connectivity
                        # Use the dummy atom's neighbor to match with the appropriate connection point
                        if str(cp_info["mcs_atom_idx"]) in atom.GetProp(
                            "atom_connecting_mcs"
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
            Chem.SanitizeMol(frag_mol)
        except:
            log_warning("Sanitization failed, the fragment may have valence issues")
            # Try to get the molecule anyway
            Chem.GetSSSR(frag_mol)

        # Return without connection_points_3d
        return mcs_mol, frag_mol, fragment_info, mcs_smarts, connection_points_3d

    def __compute_cosine_similarity(self, receptor, parent_mol, fragments):
        similarity = 0
        for fragment_info in fragments:
            fragment_smiles = fragment_info['fragment_smiles']
            branching_point = fragment_info['coordinates']

            # No cache because this calculation depends on the receptor and a specific branching point,
            # and it is strange that a branching point can be used twice in the same or different molecules
            fps_receptor_parent = self.get_prediction_for_parent_receptor(parent_mol, receptor, branching_point).tolist()

            fps_fragment = self.fps_fragment_cache.get(fragment_smiles)
            if fps_fragment is None:
                fps_fragment = self.get_fingerprints_for_fragment(Chem.MolFromSmiles(fragment_smiles)).tolist()
                self.fps_fragment_cache[fragment_smiles] = fps_fragment

            similarity = similarity + (1 - cosine(fps_receptor_parent, fps_fragment))

        similarity = (similarity / len(fragments)) if len(fragments) > 0 else 1
        return similarity
