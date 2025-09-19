"""
Functions for calculating a pruned maximum common substructure (MCS).
"""
import __future__
import time
from typing import Any, Dict, List, Optional, Tuple

SEARCH_TIMEOUT_SECONDS = 0.125  # Could be evaluating many, many candidates, so limit total time

if __name__ != "__main__":
    import autogrow.utils.mol_object_handling as MOH
    from autogrow.plugins.registry_base import plugin_managers
    from autogrow.utils.logging import log_info, log_warning

def _find_bridge_atoms(mol: Any) -> List[int]:
    """Identifies bridge atoms in a molecule.

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


def _get_fragments_for_mol(mol: Any, mcs_mol: Any) -> List[str]:
    """Get fragments by removing MCS atoms.

    Args:
        mol (Any): The molecule to fragment.
        mcs_mol (Any): The MCS to remove from the molecule.

    Returns:
        List[str]: A list of SMILES strings for the fragments.
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit
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


def _get_atoms_sorted_by_connectivity(mol: Any) -> List[Any]:
    """Returns atoms sorted by connectivity (degree), with terminal atoms first.

    Args:
        mol: RDKit molecule object

    Returns:
        list: Atoms sorted by degree (ascending), then by atom index for consistency
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit
    atoms_with_degree = []
    for atom in chemtoolkit.get_atoms(mol):
        degree = chemtoolkit.get_atom_degree(atom)
        atoms_with_degree.append((degree, chemtoolkit.get_idx(atom), atom))
    # Sort by degree (ascending, so terminal atoms first), then by index for consistency
    atoms_with_degree.sort(key=lambda x: (x[0], x[1]))
    return [atom for _, _, atom in atoms_with_degree]


def _is_valid_mcs(parent: Any, child: Any, mcs_mol: Any) -> bool:
    """
    Checks if an MCS candidate is valid by our two main rules:
    1. It must not fragment the parent or child molecule.
    2. It must not create a fragment connected by a bond that has
       inappropriately changed order (e.g., single to double).
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    # Rule 1: Check for fragmentation
    if len(_get_fragments_for_mol(parent, mcs_mol)) > 1:
        return False
    if len(_get_fragments_for_mol(child, mcs_mol)) > 1:
        return False

    # Rule 2: Check for invalid bond order changes
    parent_match = chemtoolkit.get_substruct_match(parent, mcs_mol)
    child_match = chemtoolkit.get_substruct_match(child, mcs_mol)
    if not parent_match or not child_match:
        return False # Should not happen if fragmentation check passed

    parent_to_mcs_map = {p_idx: mcs_idx for mcs_idx, p_idx in enumerate(parent_match)}
    child_to_mcs_map = {c_idx: mcs_idx for mcs_idx, c_idx in enumerate(child_match)}

    # Find connection bonds in the child and check their order
    for child_idx, mcs_idx in child_to_mcs_map.items():
        child_atom = chemtoolkit.get_atom_with_idx(child, child_idx)
        for neighbor in chemtoolkit.get_neighbors(child_atom):
            neighbor_idx = chemtoolkit.get_idx(neighbor)

            # Check if neighbor is a fragment atom (not in MCS)
            if neighbor_idx not in child_to_mcs_map:
                bond = chemtoolkit.get_bond_between_atoms(child, child_idx, neighbor_idx)
                # If the connection is a multiple bond, investigate further
                if chemtoolkit.get_bond_type(bond) != chemtoolkit.get_single_bond_type():
                    parent_idx = parent_match[mcs_idx]
                    parent_atom = chemtoolkit.get_atom_with_idx(parent, parent_idx)
                    
                    is_problem = True
                    for p_neighbor in chemtoolkit.get_neighbors(parent_atom):
                        p_neighbor_idx = chemtoolkit.get_idx(p_neighbor)
                        # Check if this neighbor in parent is also a fragment
                        if p_neighbor_idx not in parent_to_mcs_map:
                            p_bond = chemtoolkit.get_bond_between_atoms(parent, parent_idx, p_neighbor_idx)
                            # If the parent also had a bond of the same high order to a
                            # fragment, then it's not a problematic transformation.
                            if chemtoolkit.get_bond_type(p_bond) == chemtoolkit.get_bond_type(bond):
                                is_problem = False
                                break
                    if is_problem:
                        return False # This MCS is invalid

    return True


def _find_single_fragment_mcs(mol1: Any, mol2: Any) -> Optional[str]:
    """
    Finds the largest valid Maximum Common Substructure (MCS).

    This function identifies the largest common substructure between two
    molecules that is "valid" according to the `_is_valid_mcs` function.
    This includes checks for non-fragmentation and for valid bond orders
    at the connection points to ensure a chemically meaningful core.

    The process uses a two-phase approach for efficiency:
    1.  **Initial MCS & Pruning**: It first finds the absolute largest MCS. If
        this MCS is invalid, it performs a "smart pruning" step by removing
        any "bridge atoms" to create a smaller, more stable MCS candidate.
    2.  **Iterative Search**: This new, smaller candidate becomes the starting
        point for an exhaustive iterative search, which removes one atom
        at a time until a valid MCS is found.

    Args:
        mol1 (Any): The first RDKit molecule.
        mol2 (Any): The second RDKit molecule.

    Returns:
        Optional[str]: The SMARTS string of the optimal MCS, or None if no
            suitable MCS is found.
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit
    if mol1 is None or mol2 is None:
        raise ValueError("Invalid RDKit molecule provided.")

    # 1. Find the initial, absolute MCS as a starting point.
    initial_mcs = chemtoolkit.find_mcs([mol1, mol2], timeout=30)
    if initial_mcs.numAtoms == 0:
        return None

    # The list of candidates to check, starting with the biggest.
    # We store tuples of (mol_object, smarts_string).
    mcs_mol = chemtoolkit.mol_from_smarts(initial_mcs.smartsString)
    # Caching dictionary for the validity of a given MCS SMARTS. This avoids
    # re-computation.
    validity_cache: Dict[str, Optional[bool]] = {}
    # Check if the absolute largest MCS is valid from the start.
    is_valid = _is_valid_mcs(mol1, mol2, mcs_mol)
    validity_cache[initial_mcs.smartsString] = is_valid
    if is_valid:
        return initial_mcs.smartsString

    # Initial MCS is not valid, so we start the optimization search.
    log_warning(
        "Initial MCS is invalid (causes fragmentation or bond order change). "
        "Starting optimized search for a valid MCS. This may be slow."
    )
    candidates = [(mcs_mol, initial_mcs.smartsString)]
    # Prune by removing bridge atoms from the initial MCS
    bridge_atoms_in_mol1 = _find_bridge_atoms(mol1)
    bridge_atoms_in_mol2 = _find_bridge_atoms(mol2)

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
        sub_frags = chemtoolkit.get_mol_frags(pruned_mcs_mol_base, as_mols=True, sanitize_frags=False)
        if sub_frags:
            largest_sub_frag = max(sub_frags, key=lambda m: chemtoolkit.get_num_atoms(m))
            try:
                pruned_smarts = chemtoolkit.mol_to_smarts(largest_sub_frag)
                if pruned_smarts not in validity_cache:
                    candidates.append((largest_sub_frag, pruned_smarts))
                    validity_cache[pruned_smarts] = None  # Mark as seen, validity not yet known
            except:
                pass  # Ignore if SMARTS generation fails

    # Now, proceed with the iterative search, which will start with better candidates
    search_start_time = time.time()
    while candidates:
        if time.time() - search_start_time > SEARCH_TIMEOUT_SECONDS:
            log_warning(f"MCS search timed out after {SEARCH_TIMEOUT_SECONDS} seconds.")
            # DEBUG: Save two molecules to file for later analysis
            with open("mcs_timeout_debug.smi", "a") as f:
                f.write(chemtoolkit.mol_to_smiles(mol1) + "\t" + chemtoolkit.mol_to_smiles(mol2) + "\n")
            return None
        candidates.sort(key=lambda x: chemtoolkit.get_num_atoms(x[0]), reverse=True)
        current_mcs_mol, current_mcs_smarts = candidates.pop(0)
        # Check cache first. If not present (or None), compute and store.
        is_valid = validity_cache.get(current_mcs_smarts)
        if is_valid is None:
            is_valid = _is_valid_mcs(mol1, mol2, current_mcs_mol)
            validity_cache[current_mcs_smarts] = is_valid
        # Check if this candidate is valid using our comprehensive function.
        if is_valid:
            return current_mcs_smarts

        # If not valid, generate smaller candidates.
        if chemtoolkit.get_num_atoms(current_mcs_mol) > 1:
            # Get atoms sorted by connectivity (terminal atoms first)
            sorted_atoms = _get_atoms_sorted_by_connectivity(current_mcs_mol)
            for atom in sorted_atoms:
                # Create a copy to modify
                emol = chemtoolkit.get_editable_mol(current_mcs_mol)
                chemtoolkit.remove_atom_from_editable_mol(emol, chemtoolkit.get_idx(atom))
                # Removing an atom can disconnect the MCS itself. We take the largest resulting piece.
                sub_frags = chemtoolkit.get_mol_frags(chemtoolkit.get_noneditable_mol(emol), as_mols=True, sanitize_frags=False)
                if not sub_frags:
                    continue
                largest_sub_frag = max(sub_frags, key=lambda m: chemtoolkit.get_num_atoms(m))
                try:
                    new_smarts = chemtoolkit.mol_to_smarts(largest_sub_frag)
                    if new_smarts not in validity_cache:
                        # Add new, smaller candidate to our list for checking.
                        candidates.append((largest_sub_frag, new_smarts))
                        validity_cache[new_smarts] = None  # Mark as seen
                except:
                    continue
    # If the loop finishes, no suitable MCS was found.
    return None


def create_mcs_molecule(parent: Any, child: Any) -> Tuple[Optional[Any], Dict, Dict, str]:
    """
    Create a new molecule representing just the MCS with 3D coordinates from parent.

    This function orchestrates the creation of a 3D-aware MCS molecule. It
    first calls `_find_single_fragment_mcs` to determine the appropriate
    MCS SMARTS pattern, then constructs a new RDKit molecule with 3D
    coordinates inherited from the parent molecule.

    Args:
        parent (Any): The parent RDKit molecule.
        child (Any): The child RDKit molecule.

    Returns:
        Tuple[Optional[Any], Dict, Dict, str]: A tuple containing:
            - mcs_mol: The MCS molecule with 3D coordinates.
            - parent_to_mcs_map: Mapping from parent atom indices to MCS atom indices.
            - child_to_mcs_map: Mapping from child atom indices to MCS atom indices.
            - mcs_smarts: The SMARTS string of the MCS.
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit
    try:
        parent_no_h = chemtoolkit.remove_hs(parent, sanitize=False)
        parent = MOH.check_sanitization(parent_no_h)
        child_no_h = chemtoolkit.remove_hs(child, sanitize=False)
        child = MOH.check_sanitization(child_no_h)
        if parent is None or child is None:
            raise ValueError("Parent or child sanitization failed after removing hydrogens.")
    except Exception as e:
        log_warning(
            f"Failed to remove hydrogens from parent or child molecule, skipping MCS. Error: {e}"
        )
        return None, {}, {}, ""

    # Find the Maximum Common Substructure that does not fragment molecules
    mcs_smarts = _find_single_fragment_mcs(parent, child)
    if mcs_smarts is None:
        log_warning(f"Could not find single-fragment MCS for {chemtoolkit.mol_to_smiles(parent)} and {chemtoolkit.mol_to_smiles(child)}")
        return None, {}, {}, ""

    # Create a molecule from the MCS SMARTS pattern
    mcs_mol = chemtoolkit.mol_from_smarts(mcs_smarts)
    log_info(f"MCS found: {mcs_smarts}")
    # log_info(f"Number of atoms in MCS: {chemtoolkit.get_num_atoms(mcs_mol)}")
    # log_info(f"Number of bonds in MCS: {len(mcs_mol.GetBonds())}")

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


def find_mcs_and_fragments(parent: Any, child: Any, compound_id: Optional[str] = None) -> Tuple[Optional[Any], Optional[Any], List, Optional[str], Dict]:
    """
    Find the MCS, identify fragments, and determine 3D connection points.

    This is the main orchestrator for preparing a molecule for DeepFrag
    analysis. It performs three key steps:
    1. Calls `create_mcs_molecule` to obtain a 3D MCS molecule and atom
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
        Tuple[Optional[Any], Optional[Any], List, Optional[str], Dict]:
            A tuple containing:
                - The MCS molecule.
                - The fragment molecule.
                - A list of fragment information dictionaries.
                - The MCS SMARTS string.
                - A dictionary of 3D connection points for meaningful fragments.
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    # Create an explicit MCS molecule with 3D coordinates from parent
    mcs_mol, mcs_to_parent_map, mcs_to_child_map, mcs_smarts = create_mcs_molecule(
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

    # Get all fragments, which may include isolated hydrogens from the removed MCS
    all_frags_as_mols = chemtoolkit.get_mol_frags(frag_mol, as_mols=True, sanitize_frags=False)

    # Filter out fragments that are not chemically meaningful (i.e., contain no heavy atoms).
    # This prevents isolated hydrogens from being treated as fragments.
    frags = []
    for f in all_frags_as_mols:
        if any(chemtoolkit.get_atomic_num(atom) > 1 for atom in chemtoolkit.get_atoms(f)):
            frags.append(f)

    # Reconstruct the main frag_mol from only the meaningful fragments. This ensures
    # the returned molecule represents only the non-hydrogen parts of the fragment(s).
    if not frags:
        frag_mol = None
    elif len(frags) == 1:
        frag_mol = frags[0]
    else:
        combined_smiles = ".".join([chemtoolkit.mol_to_smiles(f) for f in frags])
        frag_mol = chemtoolkit.mol_from_smiles(combined_smiles, sanitize=False)

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
    
    # Filter the original connection_points_3d dictionary to only include points
    # that correspond to the valid fragments we just identified.
    final_connection_points = {}
    valid_mcs_indices = {info['mcs_atom_idx'] for info in fragment_info}
    for parent_idx, cp_info in connection_points_3d.items():
        if cp_info['mcs_atom_idx'] in valid_mcs_indices:
            final_connection_points[parent_idx] = cp_info

    # Sanitize to fix any valence issues
    if frag_mol:
        try:
            chemtoolkit.sanitize_mol(frag_mol)
        except:
            log_warning("Sanitization failed, the fragment may have valence issues")
            # Try to get the molecule anyway
            chemtoolkit.rdmolops_get_sssr(frag_mol)

    # Return the fully consistent and filtered results
    return mcs_mol, frag_mol, fragment_info, mcs_smarts, final_connection_points

if __name__ == "__main__":
    import sys
    import os

    # Add project root to path to allow for autogrow imports
    # The current file is autogrow/plugins/deepfrag_filters/common_substructure.py
    # Project root is 3 levels up.
    script_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir, os.pardir))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    # Now we can import from autogrow
    from autogrow.plugins.registry_base import plugin_managers
    import autogrow.utils.mol_object_handling as MOH
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def log_info(msg: str):
        print(f"INFO: {msg}")
    def log_warning(msg: str):
        print(f"WARNING: {msg}")

    # Setup a minimal params dict to initialize the toolkit
    params = {"RDKitToolkit": True, "obabel_path": "", "rxn_library_path": "", "procs_per_node": 0}

    # Initialize the plugin managers, which will set up the ChemToolkit
    plugin_managers.setup_plugin_managers(params)
    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    def test_find_mcs_and_fragments(parent_child_smiles: Tuple[str, str]):
        # Test1: Simple benzene and toluene
        mol1 = chemtoolkit.mol_from_smiles(parent_child_smiles[0])
        mol2 = chemtoolkit.mol_from_smiles(parent_child_smiles[1])

        # Add hydrogens and generate 3D conformers for testing
        mol1 = chemtoolkit.add_hs(mol1)
        AllChem.EmbedMolecule(mol1)

        mol2 = chemtoolkit.add_hs(mol2)
        AllChem.EmbedMolecule(mol2)

        return find_mcs_and_fragments(mol1, mol2)

    def test_create_mcs_molecule(parent_child_smiles: Tuple[str, str]):
        mol1 = chemtoolkit.mol_from_smiles(parent_child_smiles[0])
        mol2 = chemtoolkit.mol_from_smiles(parent_child_smiles[1])

        # Add hydrogens and generate 3D conformers for testing
        mol1 = chemtoolkit.add_hs(mol1)
        AllChem.EmbedMolecule(mol1)

        mol2 = chemtoolkit.add_hs(mol2)
        AllChem.EmbedMolecule(mol2)

        return create_mcs_molecule(mol1, mol2)

    # Test1: Simple benzene and toluene
    parent_child_smiles = ("c1ccccc1", "c1ccccc1C")
    mcs_mol, frag_mol, frag_info_dict, mcs_smarts, con_pts_dict = test_find_mcs_and_fragments(parent_child_smiles)
    assert Chem.MolToSmiles(mcs_mol) == "c1ccccc1"
    assert Chem.MolToSmiles(frag_mol) == "*C([H])([H])[H]"
    assert len(frag_info_dict) == 1
    assert frag_info_dict[0]["fragment_smiles"] == "*C([H])([H])[H]"
    assert mcs_smarts == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    assert len(con_pts_dict) == 1

    mcs_mol, parent_to_mcs_map, child_to_mcs_map, mcs_smarts = test_create_mcs_molecule(parent_child_smiles)
    assert Chem.MolToSmiles(mcs_mol) == "c1ccccc1"
    assert len(parent_to_mcs_map) == len(child_to_mcs_map) == 6
    assert mcs_smarts == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"

    # Test2: MCS would produce two fragments
    parent_child_smiles = ("c1ccccc1C", "c1ccccc1C(C)C")
    mcs_mol, frag_mol, frag_info_dict, mcs_smarts, con_pts_dict = test_find_mcs_and_fragments(parent_child_smiles)
    assert Chem.MolToSmiles(mcs_mol) == "c1ccccc1"
    assert Chem.MolToSmiles(frag_mol) == "*C([H])(C([H])([H])[H])C([H])([H])[H]"
    assert len(frag_info_dict) == 1
    assert frag_info_dict[0]["fragment_smiles"] == "*C([H])(C([H])([H])[H])C([H])([H])[H]"
    assert mcs_smarts == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    assert len(con_pts_dict) == 1

    mcs_mol, parent_to_mcs_map, child_to_mcs_map, mcs_smarts = test_create_mcs_molecule(parent_child_smiles)
    assert Chem.MolToSmiles(mcs_mol) == "c1ccccc1"
    assert len(parent_to_mcs_map) == len(child_to_mcs_map) == 6
    assert mcs_smarts == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1" 

    # NOTE: This scheme does allow for fragments with multiple connection
    # points. I believe the code handles this by averaging the deepfrag
    # predictions at each point. Not a great solution IMO, but should be a
    # fairly rare edge case. And technically allows DeepFrag to be applied to
    # crossover compounds, I think.
    # parent_child_smiles = ("IC1CCCCC1", "IC1C2C(CCC1)CCCC2")
    # mcs_mol, frag_mol, frag_info_dict, mcs_smarts, con_pts_dict = test_find_mcs_and_fragments(parent_child_smiles)

    parent_child_smiles = ("C1CCCCC1C", "C1CCCCC1C=C")
    mcs_mol, frag_mol, frag_info_dict, mcs_smarts, con_pts_dict = test_find_mcs_and_fragments(parent_child_smiles)
    assert Chem.MolToSmiles(mcs_mol) == "C1CCCCC1"
    assert Chem.MolToSmiles(frag_mol) == "*C([H])=C([H])[H]"
    assert len(frag_info_dict) == 1
    assert frag_info_dict[0]["fragment_smiles"] == "*C([H])=C([H])[H]"
    assert mcs_smarts == "[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1"
    assert len(con_pts_dict) == 1

    mcs_mol, parent_to_mcs_map, child_to_mcs_map, mcs_smarts = test_create_mcs_molecule(parent_child_smiles)
    assert Chem.MolToSmiles(mcs_mol) == "C1CCCCC1"  # This one fails. Chem.MolToSmiles(mcs_mol) => 'CC1CCCCC1' Why?
    assert len(parent_to_mcs_map) == len(child_to_mcs_map) == 6
    assert mcs_smarts == "[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1"

    # Try a case that will take a long time to find a valid MCS
    # print("\n\nHI\n\n")
    parent_child_smiles = ("CC(=O)N(CC(C)F)C(=O)NC(C)C(=O)Oc1n[nH]c(=O)c2ccccc12", "CC(=O)N(CC(C)OC(=O)NC(=S)NCC(=O)OCCC=CO)C(=O)NC(C)C(=O)Oc1n[nH]c(=O)c2ccccc12")
    mcs_mol, frag_mol, frag_info_dict, mcs_smarts, con_pts_dict = test_find_mcs_and_fragments(parent_child_smiles)


    print("\nAll tests passed.")

