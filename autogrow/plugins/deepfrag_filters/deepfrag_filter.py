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
from autogrow.utils.logging import log_debug

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class DeepFragFilterBase(PluginBase):
    """
    Abstract base class for DeepFrag plugins.

    This class defines the interface for DeepFrag plugins and provides a common
    run method that calls the abstract run_deepfrag_filter method.
    """

    def run(self, **kwargs) -> List[Compound]:
        """
        Run the DeepFrag plugin.

        Args:
            predocked_cmpds (List[Compound]): A list of Compound objects after filtering
            out molecules that are not suitable for the docking according to the DeepFrag
            filter.
            cutoff: the cosine similarity value to be considered as cutoff.

        Returns:
            List[Compound]: A list of Compound objects.
        """
        compounds = kwargs["compounds"]
        cutoff = kwargs["input_params"][self.name]
        receptor = kwargs["input_params"]["receptor_path"]

        final_compound_list = []
        for compound in compounds:
            if len(compound.parent_3D_mols) == 1:
                mcs_mol, _, fragments = self.__find_mcs_and_fragments(compound.parent_3D_mols[0], Chem.MolFromSmiles(compound.smiles))
                similarity = self.__compute_cosine_similarity(receptor, mcs_mol, fragments)
                passed_filter = similarity >= cutoff
            elif len(compound.parent_3D_mols) == 2:
                mcs_mol, _, fragments = self.__find_mcs_and_fragments(compound.parent_3D_mols[0], Chem.MolFromSmiles(compound.smiles))
                similarity_0 = self.__compute_cosine_similarity(receptor, mcs_mol, fragments)
                mcs_mol, _, fragments = self.__find_mcs_and_fragments(compound.parent_3D_mols[1], Chem.MolFromSmiles(compound.smiles))
                similarity_1 = self.__compute_cosine_similarity(receptor, mcs_mol, fragments)
                passed_filter = similarity_0 >= cutoff or similarity_1 >= cutoff

            compound.mol_3D = None
            compound.parent_3D_mols = None

            if passed_filter:
                final_compound_list.append(compound)
            else:
                log_debug(
                    f"Docked molecule {compound.id} with smiles string {compound.smiles} did not fulfill with the similarity criterion using DeepFrag"
                )

        return final_compound_list

    @abstractmethod
    def get_prediction_for_parent_receptor(self, parent_mol, receptor, branching_point):
        pass

    @abstractmethod
    def get_fingerprints_for_fragment(self, fragment):
        pass

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass

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
                    help="An value representing the cosine similarity value to be "
                         "considered as cutoff in order to a molecule pass the filter or "
                         "not",
                )
            ],
        )

    # Create a new MCS molecule with 3D coordinates from parent
    def __create_mcs_molecule(self, parent, child):
        """
        Create a new molecule representing just the MCS with 3D coordinates from parent

        Returns:
        - mcs_mol: The MCS molecule with 3D coordinates
        - parent_to_mcs_map: Mapping from parent atom indices to MCS atom indices
        - child_to_mcs_map: Mapping from child atom indices to MCS atom indices
        """
        parent = Chem.RemoveHs(parent)
        child = Chem.RemoveHs(child)

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
    def __find_mcs_and_fragments(self, parent, child):

        # Create an explicit MCS molecule with 3D coordinates from parent
        mcs_mol, mcs_to_parent_map, mcs_to_child_map = self.__create_mcs_molecule(parent, child)

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
                            'mcs_atom_idx': mcs_atom_idx,
                            'neighbor_symbol': child.GetAtomWithIdx(neighbor_idx).GetSymbol(),
                            'coordinates': parent_conf.GetAtomPosition(parent_atom_idx)
                            # Only needed for matching later
                        }

                    attachment_bonds.append((child_atom_idx, neighbor_idx))

        for atom in rwmol.GetAtoms():
            atom_idx = atom.GetIdx()
            if atom_idx in child_to_mcs_map:
                atom.SetProp("is_mcs", "yes")
            else:
                atom.SetProp("is_mcs", "no")

        # Create dummy atoms at attachment points
        dummy_atoms = []
        for child_atom_idx, neighbor_idx in attachment_bonds:
            # Add a dummy atom (R group)
            dummy_idx = rwmol.AddAtom(Chem.Atom('*'))
            # Add a bond from the dummy atom to the neighbor atom
            rwmol.AddBond(dummy_idx, neighbor_idx, Chem.BondType.SINGLE)
            dummy_atoms.append((dummy_idx, child_atom_idx))  # Store the dummy atom and its corresponding MCS atom

            # Find the parent atom index that corresponds to this child atom index
            for parent_idx, connection_info in connection_points_3d.items():
                if connection_info.get('mcs_atom_idx') == child_to_mcs_map.get(child_atom_idx):
                    connection_info['dummy_idx'] = dummy_idx
                    atom = rwmol.GetAtomWithIdx(dummy_idx)
                    atom.SetProp('atom_connecting_mcs', str(connection_info.get('mcs_atom_idx')))

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
        frags = [Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)[i] for i in range(len(atom_frags))]

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
                    # dummy_neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]

                    # Try to match this fragment to a connection point
                    matched_connection = False

                    for cp_info in connection_points_3d.values():
                        # Generic matching based on atom indices and connectivity
                        # Use the dummy atom's neighbor to match with the appropriate connection point
                        if str(cp_info['mcs_atom_idx']) in atom.GetProp('atom_connecting_mcs'):
                            connection_info = {
                                'fragment_smiles': frag_smiles,
                                'mcs_atom_idx': cp_info['mcs_atom_idx'],
                                'coordinates': cp_info['coordinates']
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

    def __compute_cosine_similarity(self, receptor, parent_mol, fragments):
        similarity = 0
        for fragment_info in fragments:
            fragment_smiles = fragment_info['fragment_smiles']
            branching_point = fragment_info['coordinates']

            fps_receptor_parent = self.get_prediction_for_parent_receptor(parent_mol, receptor, branching_point).tolist()
            fps_fragment = self.get_fingerprints_for_fragment(Chem.MolFromSmiles(fragment_smiles)).tolist()

            similarity = similarity + (1 - cosine(fps_receptor_parent, fps_fragment))

        similarity = (similarity / len(fragments)) if len(fragments) > 0 else 1
        return similarity
