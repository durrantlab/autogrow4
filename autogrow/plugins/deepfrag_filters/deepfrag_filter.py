"""
DeepFrag plugin.
"""
import __future__

import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import rdFMCS
from rdkit.Geometry import Point3D
from typing import List, Tuple
from autogrow.types import Compound
from autogrow.plugins.plugin_base import PluginBase
import copy
from abc import abstractmethod
from autogrow.config.argument_vars import ArgumentVars
from scipy.spatial.distance import cosine

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
            passed_filter = False
            if len(compound.parent_3D_mols) == 1:
                parent_mol, branching_points, fragment_mols = self.__get_substructure_and_branching_points(compound.parent_3D_mols[0], compound)
                passed_filter = self.__compute_cosine_similarity(receptor, parent_mol, branching_points, fragment_mols) >= cutoff
            elif len(compound.parent_3D_mols) == 2:
                parent_mol, branching_points_0, fragment_mols_0 = self.__get_substructure_and_branching_points(compound.parent_3D_mols[0], compound)
                parent_mol, branching_points_1, fragment_mols_1 = self.__get_substructure_and_branching_points(compound.parent_3D_mols[1], compound)
                similarity_0 = self.__compute_cosine_similarity(receptor, parent_mol, branching_points_0, fragment_mols_0)
                similarity_1 = self.__compute_cosine_similarity(receptor, parent_mol, branching_points_1, fragment_mols_1)
                passed_filter = similarity_0 >= cutoff or similarity_1 >= cutoff

            compound.mol_3D = None
            compound.parent_3D_mols = None

            if passed_filter:
                final_compound_list.append(compound)

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

    def __compute_cosine_similarity(self, receptor, parent_mol, branching_points, fragment_mols):
        if parent_mol is None or receptor is None or branching_points is None or fragment_mols is None:
            return 0

        assert len(branching_points) == len(fragment_mols)
        similarity = 0
        for idx in range(len(branching_points)):
            fragment_mol = fragment_mols[idx]
            branching_point = branching_points[idx]

            fps_receptor_parent = self.get_prediction_for_parent_receptor(parent_mol, receptor, branching_point).tolist()
            fps_fragment = self.get_fingerprints_for_fragment(fragment_mol).tolist()

            similarity = similarity + (1 - cosine(fps_receptor_parent, fps_fragment))

        similarity = similarity / len(branching_points)
        return similarity

    def __get_substructure_and_branching_points(self, parent, predocked_cmpd):
        parent = Chem.RemoveHs(parent)
        parent_atoms_amount = parent.GetNumAtoms()

        predocked_cmpd.read_3D_structure_sdf()
        child = predocked_cmpd.mol_3D
        try:
            child = Chem.RemoveHs(child)
        except:
            # this exception is due to the fact that the .sdf file wrote after
            # docking cannot be read successfully.
            return None, None, None
        child_atoms_amount = child.GetNumAtoms()

        mcs_res = rdFMCS.FindMCS([parent, child])

        child_mcs_atom_indices = child.GetSubstructMatch(Chem.MolFromSmarts(mcs_res.smartsString), useChirality=False, useQueryQueryMatches=False)
        child_mcs_mol, _ = self.__get_substructure_with_coords(child, copy.deepcopy(child), child_mcs_atom_indices)
        branching_points = self.__get_branching_points(child_mcs_mol)

        if len(branching_points) == 0:
            assert child_atoms_amount <= parent_atoms_amount
            return None, None, None
        else:
            fragments, branching_points = self.__get_fragments_for_branching_points(copy.deepcopy(child_mcs_mol), branching_points, child_mcs_atom_indices)
            parent_mol, point_coords = self.__get_substructure_with_coords(child, Chem.MolFromSmarts(mcs_res.smartsString), child_mcs_atom_indices, True, branching_points)
            return parent_mol, point_coords, fragments

    def __get_fragments_for_branching_points(self, mol, branching_points, mcs_atom_idxs):
        branching_points = list(branching_points)
        branching_points.sort()
        bonds_to_remove = []
        atoms_by_branching_points = {}

        # Search the bonds where an atom belongs to the MCS and the another one does not
        for atom in mol.GetAtoms():
            assert atom.HasProp("is_mcs")
            if atom.GetProp("is_mcs") == "yes":
                # Get all the neighbor atoms for an atom belonging to the MCS
                for neighbor_atom in atom.GetNeighbors():
                    # Check that the neighbor atom does not belong to the MCS.
                    # if true, then remove the bond
                    if neighbor_atom.GetProp("is_mcs") == "no":
                        bond_to_remove = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor_atom.GetIdx())
                        bonds_to_remove.append(bond_to_remove.GetIdx())
                        if atom.GetIdx() not in atoms_by_branching_points:
                            atoms_by_branching_points[atom.GetIdx()] = set()
                        atoms_by_branching_points[atom.GetIdx()].add(neighbor_atom.GetIdx())

        # Convert fragments to separate molecules
        mol = Chem.FragmentOnBonds(mol, bonds_to_remove)

        # Convert tuple of indexes representing the atoms belonging to the mcs
        mcs_atom_idxs = set(mcs_atom_idxs)
        # Get the list of indexes of the atoms belonging to the fragments obtained after removing the bonds
        list_fragment_idxs = list(Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False))
        # Get the list of RDKit molecules corresponding to the fragments obtained after removing the bonds
        list_fragment_mols = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))

        # Search and remove the fragment matching with the parent molecule
        for pos, fragment_idxs in enumerate(list_fragment_idxs):
            fragment_atom_set = set(fragment_idxs)
            dif_fragment_parent = fragment_atom_set.intersection(mcs_atom_idxs)
            if len(dif_fragment_parent) != 0:
                del list_fragment_idxs[pos]
                del list_fragment_mols[pos]
                break

        # Ensure that the list of fragments to be returned is in the same order of the branching points
        if len(list_fragment_mols) > 0:
            c_list_fragment_mols = list_fragment_mols
            c_branching_points = branching_points
            list_fragment_mols = []
            branching_points = []
            for point in c_branching_points:
                bonded_atoms = atoms_by_branching_points[point]
                for pos, atoms_for_frag in enumerate(list_fragment_idxs):
                    if set(bonded_atoms).intersection(set(atoms_for_frag)):
                        list_fragment_mols.append(c_list_fragment_mols[pos])
                        branching_points.append(point)
                        del list_fragment_idxs[pos]
                        break

        assert len(branching_points) == len(list_fragment_mols)
        return list_fragment_mols, branching_points

    def __get_branching_points(self, mol):
        # Set to save all the possible branching points
        branching_point_idxs = set()

        for atom in mol.GetAtoms():
            assert atom.HasProp("is_mcs")
            # Check if the atom does not belong to the MCS.
            if atom.GetProp("is_mcs") == "no":
                # Get all neighbors
                for neighbor_atom in atom.GetNeighbors():
                    # If neighbor atom belong to the MCS, then it is a branching point
                    if neighbor_atom.GetProp("is_mcs") == "yes":
                        # Add the neighbor atom to the list to be returned
                        branching_point_idxs.add(neighbor_atom.GetIdx())

        return branching_point_idxs

    def __get_substructure_with_coords(self, mol, copy_mol, mcs_atom_indices, do_enumerate=False, branching_points=None):

        # Coordinates of the branching point
        branching_point_coord = []

        # Get the conformer from mol
        conf = mol.GetConformer()

        # Create new mol
        new_mol = Chem.RWMol(copy_mol)

        # Create conformer for new mol
        new_conf = Chem.Conformer(new_mol.GetNumAtoms())
        for atom in new_mol.GetAtoms():
            if atom.GetIdx() in mcs_atom_indices:
                atom.SetProp("is_mcs", "yes")
            else:
                atom.SetProp("is_mcs", "no")
            new_conf.SetAtomPosition(atom.GetIdx(), Point3D(0, 0, 0))

        # Set the coordinates
        if not do_enumerate:
            for atom_idx in mcs_atom_indices:
                new_conf.SetAtomPosition(atom_idx, conf.GetAtomPosition(atom_idx))
        else:
            for idx, atom_idx in enumerate(mcs_atom_indices):
                new_conf.SetAtomPosition(idx, conf.GetAtomPosition(atom_idx))
                if atom_idx in branching_points:
                    branching_point_coord.append(conf.GetAtomPosition(atom_idx))

        # Add new conf
        new_mol.RemoveConformer(0)
        new_mol.AddConformer(new_conf)

        # Convert to mol
        new_mol = new_mol.GetMol()

        return new_mol, branching_point_coord
