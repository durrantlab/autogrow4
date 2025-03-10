from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS, Descriptors, rdmolops
import numpy as np
from collections import defaultdict

# Get the PDB atom names from a molecule
def get_pdb_atom_names(molecule):
    """Extract PDB atom names from a molecule's PDB block"""
    pdb_block = Chem.MolToPDBBlock(molecule)
    pdb_lines = pdb_block.split('\n')
    atom_names = {}
    
    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                atom_idx = int(line[6:11].strip()) - 1  # PDB is 1-indexed
                atom_name = line[12:16].strip()
                atom_names[atom_idx] = atom_name
            except (ValueError, IndexError):
                continue
    
    return atom_names

# Create a new MCS molecule with 3D coordinates from parent
def create_mcs_molecule(parent, child):
    """
    Create a new molecule representing just the MCS with 3D coordinates from parent
    
    Returns:
    - mcs_mol: The MCS molecule with 3D coordinates
    - parent_to_mcs_map: Mapping from parent atom indices to MCS atom indices
    - child_to_mcs_map: Mapping from child atom indices to MCS atom indices
    """
    # Find the Maximum Common Substructure
    mcs = rdFMCS.FindMCS([parent, child], 
                         completeRingsOnly=True,
                         ringMatchesRingOnly=True,
                         matchValences=True)
    
    # Create a molecule from the MCS SMARTS pattern
    mcs_smarts = mcs.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    
    print(f"MCS found: {mcs.smartsString}")
    print(f"Number of atoms in MCS: {mcs.numAtoms}")
    print(f"Number of bonds in MCS: {mcs.numBonds}")
    
    # Match the MCS in both molecules
    parent_match = parent.GetSubstructMatch(mcs_mol)
    child_match = child.GetSubstructMatch(mcs_mol)
    
    if not parent_match or not child_match:
        print("No match found in one or both molecules.")
        return None, {}, {}
    
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
        print("Warning: MCS molecule sanitization failed")
    
    # Add MCS to parent atom mapping
    mcs_to_parent_map = {atom_mapping[i]: parent_match[i] for i in range(len(parent_match))}
    
    # Add MCS to child atom mapping
    mcs_to_child_map = {atom_mapping[i]: child_match[i] for i in range(len(child_match))}
    
    return new_mcs_mol, mcs_to_parent_map, mcs_to_child_map

# Function to find MCS and remove it from the second molecule
def find_mcs_and_fragments(parent, child):
    # Create an explicit MCS molecule with 3D coordinates from parent
    mcs_mol, mcs_to_parent_map, mcs_to_child_map = create_mcs_molecule(parent, child)
    
    if mcs_mol is None:
        print("Failed to create MCS molecule.")
        return None, None, []
    
    # Create reverse mappings
    child_to_mcs_map = {v: k for k, v in mcs_to_child_map.items()}
    
    # Create an RWMol for editing
    rwmol = Chem.RWMol(child)
    
    # Find atoms connected to the MCS but not part of it (meaning, attachment points)
    attachment_bonds = []
    connection_points_3d = {}  # Store connection info
    
    # Find attachment points using the MCS molecule and mappings
    for mcs_atom_idx, child_atom_idx in mcs_to_child_map.items():
        child_atom = child.GetAtomWithIdx(child_atom_idx)
        
        for neighbor in child_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in child_to_mcs_map:
                # This is a bond to an attachment point
                # Find the corresponding atom in the parent using the MCS mapping
                parent_atom_idx = mcs_to_parent_map.get(mcs_atom_idx)
                
                if parent_atom_idx is not None and parent.GetNumConformers() > 0:
                    # Store only the necessary connection details
                    connection_points_3d[parent_atom_idx] = {
                        'mcs_atom_idx': mcs_atom_idx,
                        'neighbor_symbol': child.GetAtomWithIdx(neighbor_idx).GetSymbol() # Only needed for matching later
                    }
                
                attachment_bonds.append((child_atom_idx, neighbor_idx))
    
    # Create dummy atoms at attachment points
    dummy_atoms = []
    for atom_idx, neighbor_idx in attachment_bonds:
        # Add a dummy atom (R group)
        dummy_idx = rwmol.AddAtom(Chem.Atom('*'))
        # Add a bond from the dummy atom to the neighbor atom
        rwmol.AddBond(dummy_idx, neighbor_idx, Chem.BondType.SINGLE)
        dummy_atoms.append((dummy_idx, atom_idx))  # Store the dummy atom and its corresponding MCS atom
        
        # Find the parent atom index that corresponds to this child atom index
        for parent_idx, connection_info in connection_points_3d.items():
            if connection_info.get('mcs_atom_idx') == mcs_to_child_map.get(atom_idx):
                connection_info['dummy_idx'] = dummy_idx
    
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
    frags = [Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)[i] for i in range(len(atom_frags))]
    
    # For each fragment, identify which dummy atom it contains
    fragment_info = []
    for i, frag in enumerate(frags):
        frag_smiles = Chem.MolToSmiles(frag)
        
        # We'll track which connections this fragment has
        fragment_connections = []
        
        # Look for dummy atoms in this fragment
        for atom in frag.GetAtoms():
            if atom.GetSymbol() == '*':
                # For each dummy atom, identify its connections
                dummy_neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
                
                # Try to match this fragment to a connection point
                matched_connection = False
                
                for cp_info in connection_points_3d.values():
                    # Generic matching based on atom indices and connectivity
                    # Use the dummy atom's neighbor to match with the appropriate connection point
                    if cp_info['neighbor_symbol'] in dummy_neighbors:
                        
                        connection_info = {
                            'fragment_smiles': frag_smiles,
                            'mcs_atom_idx': cp_info['mcs_atom_idx']
                        }
                        
                        fragment_connections.append(connection_info)
                        matched_connection = True
                        break  # Take only the first match
                
                if not matched_connection:
                    print(f"Warning: Could not match fragment {frag_smiles} to a connection point")
        
        # Add all connection points for this fragment
        fragment_info.extend(fragment_connections)
    
    # Sanitize to fix any valence issues
    try:
        Chem.SanitizeMol(frag_mol)
    except:
        print("Warning: Sanitization failed, the fragment may have valence issues")
        # Try to get the molecule anyway
        Chem.GetSSSR(frag_mol)
    
    # Return without connection_points_3d
    return mcs_mol, frag_mol, fragment_info

def find_equivalent_atoms(mcs_mol, child_mol):
    """
    Creates a dictionary mapping each atom in MCS to its equivalent atoms in child molecule.
    
    Parameters:
    mcs_mol (rdkit.Chem.rdchem.Mol): rdkit.Chem.rdchem.Mol object of the MCS
    child_mol (rdkit.Chem.rdchem.Mol): rdkit.Chem.rdchem.Mol object of the child molecule
    
    Returns:
    dict: Contains:
        - is_substructure (bool): True if MCS is a substructure of child
        - mcs_to_child_map (dict): Dictionary where keys are atom indices in MCS
          and values are lists of equivalent atom indices in child molecule
        - num_matches (int): Number of matches found
    """
    
    # Check if the molecules were parsed correctly
    if mcs_mol is None or child_mol is None:
        return {
            "error": "Error parsing SMILES strings",
            "is_substructure": False,
            "num_matches": 0,
            "mcs_to_child_map": {}
        }
    
    # Get all possible substructure matches (including symmetric alternatives)
    matches = child_mol.GetSubstructMatches(mcs_mol, uniquify=False)
    
    # Check if any matches were found
    is_substructure = len(matches) > 0
    
    # Create dictionary mapping each MCS atom to all its equivalent atoms in child
    mcs_to_child_map = defaultdict(set)
    
    for match in matches:
        for mcs_idx, child_idx in enumerate(match):
            mcs_to_child_map[mcs_idx].add(child_idx)
    
    # Convert sets to sorted lists for cleaner output
    mcs_to_child_map = {k: sorted(list(v)) for k, v in mcs_to_child_map.items()}
    
    return {
        "is_substructure": is_substructure,
        "num_matches": len(matches),
        "mcs_to_child_map": mcs_to_child_map
    }

def main():
    # Create molecules directly from SMILES
    parent_smiles = "c1ccc(CC(=O)O)cc1"  # Phenylacetic acid
    child_smiles = "O=C(OCCO)Cc1cc(F)ccc1"
    
    # Create the parent molecule with 3D coordinates
    parent_3d_mol = Chem.MolFromSmiles(parent_smiles)
    parent_3d_mol = Chem.AddHs(parent_3d_mol)  # Add hydrogens for better 3D structure
    AllChem.EmbedMolecule(parent_3d_mol, randomSeed=42)  # Generate 3D coordinates
    AllChem.UFFOptimizeMolecule(parent_3d_mol)  # Optimize the 3D structure
    parent_3d_mol = Chem.RemoveHs(parent_3d_mol)  # Remove hydrogens for MCS matching
    
    # Create the child molecule without 3D coordinates
    child_1d_mol = Chem.MolFromSmiles(child_smiles)
    
    print("Molecules created successfully!")
    print(f"3D molecule (parent): {Chem.MolToSmiles(parent_3d_mol)}")
    print(f"1D molecule (child): {Chem.MolToSmiles(child_1d_mol)}")

    # Find MCS and fragments
    mcs_mol, fragments, fragment_info = find_mcs_and_fragments(parent_3d_mol, child_1d_mol)

    assert mcs_mol, "MCS molecule not found"
    assert fragments, "Fragments not found"

    # Get the SMILES of the MCS molecule
    mcs_smiles = Chem.MolToSmiles(mcs_mol)
    print(f"\nMCS molecule: {mcs_smiles}")
        
    # Get the SMILES with the R groups (asterisks)
    frag_smiles = Chem.MolToSmiles(fragments)
    print(f"\nFragments after removing MCS: {frag_smiles}")

    # Get all equivalent atoms
    equivalent_atoms = find_equivalent_atoms(mcs_mol, parent_3d_mol)
    
    # Print individual fragments with their connection point coordinates
    print("\n==== FRAGMENT SUMMARY ====\n")
    for info in fragment_info:
        mcs_atom_idx = info['mcs_atom_idx']
        equivalent_parent_atoms_idxs = equivalent_atoms["mcs_to_child_map"][mcs_atom_idx]

        for parent_atom_idx in equivalent_parent_atoms_idxs:
            coords = parent_3d_mol.GetConformer().GetAtomPosition(parent_atom_idx)

            print(f"Fragment: {info['fragment_smiles']}")
            print(f"  Coordinates: ({coords.x:.4f}, {coords.y:.4f}, {coords.z:.4f})")
            print("")

    # Save the MCS to a PDB for easy visualization
    mcs_pdb = Chem.MolToPDBBlock(mcs_mol)
    with open("mcs.pdb", "w") as f:
        f.write(mcs_pdb)

if __name__ == "__main__":
    main()