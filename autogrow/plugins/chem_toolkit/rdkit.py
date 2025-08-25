"""RDKit implementation of chemistry toolkit plugin."""
from typing import Any, Dict, List, Tuple, Optional
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.chem_toolkit import ChemToolkitBase

import rdkit
import rdkit.Chem as Chem
from rdkit import DataStructs
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.Chem import rdFMCS
import rdkit.Chem.MolSurf as MolSurf  # type: ignore
from rdkit.Chem import FilterCatalog  # type: ignore
from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
import rdkit.Chem.Lipinski as Lipinski  # type: ignore
import rdkit.Chem.Crippen as Crippen  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
from rdkit.Chem.BRICS import BRICSDecompose  # type: ignore
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
rdkit.RDLogger.DisableLog("rdApp.*")


class RDKitToolkit(ChemToolkitBase):
    """RDKit implementation of chemistry toolkit."""

    def add_arguments(self) -> List[ArgumentVars]:
        """
        Add command-line arguments specific to the RDkit toolkit.

        Returns:
            List[ArgumentVars]: A list with one ArgumentVars object defining the
            argument to enable the rdkit chemistry toolkit.
        """
        return [
            ArgumentVars(
                name=self.name,
                action="store_true",
                default=False,
                help="Use the RDKit library as the backend for all cheminformatics operations, such as molecule manipulation and fingerprint generation.",
            )
        ]

    def mol_from_smiles(self, smiles: str, sanitize: bool = True) -> Any:
        """Create a molecule from SMILES."""
        return Chem.MolFromSmiles(smiles, sanitize=sanitize)

    def mol_to_smiles(
        self, mol: Any, canonical: bool = True, isomeric_smiles: bool = True
    ) -> str:
        """Convert molecule to SMILES."""
        return Chem.MolToSmiles(
            mol, canonical=canonical, isomericSmiles=isomeric_smiles
        )

    def mol_from_smarts(self, smarts: str) -> Any:
        """Create a molecule from SMARTS."""
        return Chem.MolFromSmarts(smarts)

    def reaction_from_smarts(self, smarts: str) -> Any:
        """Create a reaction from SMARTS."""
        return AllChem.ReactionFromSmarts(smarts)

    def sanitize_mol(self, mol: Any, catch_errors: bool = False) -> Tuple[Any, Any]:
        """Sanitize a molecule."""
        resp = Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
            catchErrors=catch_errors,
        )
        return mol, resp

    def add_hs(self, mol: Any) -> Any:
        """Add hydrogens to a molecule."""
        return Chem.AddHs(mol)

    def remove_hs(self, mol: Any, sanitize: bool = True) -> Any:
        """Remove hydrogens from a molecule."""
        return Chem.RemoveHs(mol, sanitize=sanitize)

    def find_mcs(
        self,
        mols: List[Any],
        match_valences: bool = False,
        ring_matches_ring_only: bool = False,
        complete_rings_only: bool = False,
        timeout: int = 3600,
    ) -> Any:
        """Find maximum common substructure."""
        return rdFMCS.FindMCS(
            mols,
            matchValences=match_valences,
            ringMatchesRingOnly=ring_matches_ring_only,
            completeRingsOnly=complete_rings_only,
            timeout=timeout,
        )

    def get_mol_frags(
        self, mol: Any, as_mols: bool = False, sanitize_frags: bool = True
    ) -> List[Any]:
        """Get molecular fragments."""
        return Chem.GetMolFrags(mol, asMols=as_mols, sanitizeFrags=sanitize_frags)

    def get_morgan_fingerprint(
        self, mol: Any, radius: int, use_features: bool = False
    ) -> Any:
        """Generate Morgan fingerprint."""
        return GetMorganFingerprint(mol, radius, useFeatures=use_features)

    def dice_similarity(self, fp1: Any, fp2: Any) -> float:
        """Calculate Dice similarity."""
        return DataStructs.DiceSimilarity(fp1, fp2)

    def get_atoms(self, mol: Any) -> List[Any]:
        """Get atoms from molecule."""
        return mol.GetAtoms()

    def get_atomic_num(self, atom: Any) -> int:
        """Get atomic number."""
        return atom.GetAtomicNum()

    def get_formal_charge(self, atom: Any) -> int:
        """Get formal charge."""
        return atom.GetFormalCharge()

    def set_formal_charge(self, atom: Any, charge: int) -> None:
        """Set formal charge."""
        atom.SetFormalCharge(charge)

    def get_bonds(self, atom: Any) -> List[Any]:
        """Get bonds for atom."""
        return atom.GetBonds()

    def get_bond_type_as_double(self, bond: Any) -> float:
        """Get bond type as double value."""
        return bond.GetBondTypeAsDouble()

    def remove_atoms(self, mol: Any, atoms_to_remove: List[int]) -> Any:
        """Remove atoms from molecule."""
        atoms_to_remove.sort(reverse=True)
        em1 = Chem.EditableMol(mol)
        for atom in atoms_to_remove:
            em1.RemoveAtom(atom)

        return em1.GetMol()

    def descriptors_exact_mol_wt(self, mol: Any) -> float:
        """Get exact molecular weight."""
        return Chem.Descriptors.ExactMolWt(mol)

    def molsurf_tpsa(self, mol: Any) -> float:
        """Get topological polar surface area."""
        return MolSurf.TPSA(mol)

    def get_pains_a_filter(self) -> Any:
        """Get PAINS-A filter."""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        return FilterCatalog.FilterCatalog(params)

    def get_pains_b_filter(self) -> Any:
        """Get PAINS-B filter."""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        return FilterCatalog.FilterCatalog(params)

    def get_pains_c_filter(self) -> Any:
        """Get PAINS-C filter."""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        return FilterCatalog.FilterCatalog(params)

    def get_pains_filter(self) -> Any:
        """Get general PAINS filter."""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        return FilterCatalog.FilterCatalog(params)

    def get_nih_filter(self) -> Any:
        """Get NIH filter."""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        return FilterCatalog.FilterCatalog(params)

    def get_brenk_filter(self) -> Any:
        """Get BRENK filter."""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        return FilterCatalog.FilterCatalog(params)

    def lipinski_heavy_atom_count(self, mol: Any) -> int:
        """Get number of heavy atoms."""
        return Lipinski.HeavyAtomCount(mol)

    def lipinski_num_rotatable_bonds(self, mol: Any) -> int:
        """Get number of rotatable bonds."""
        return Lipinski.NumRotatableBonds(mol)

    def lipinski_num_h_donors(self, mol: Any) -> int:
        """Get number of H donors."""
        return Lipinski.NumHDonors(mol)

    def lipinski_num_h_acceptors(self, mol: Any) -> int:
        """Get number of H donors."""
        return Lipinski.NumHAcceptors(mol)

    def rdmolops_get_sssr(self, mol: Any) -> int:
        """Get number of SSSR rings."""
        return Chem.rdmolops.GetSSSR(mol)

    def crippen_mol_log_p(self, mol: Any) -> float:
        """Get molar LogP."""
        return Crippen.MolLogP(mol)

    def crippen_mol_mr(self, mol: Any) -> float:
        """Get molar refractivity."""
        return Crippen.MolMR(mol)

    def renumber_atoms(self, mol: Any, full_tuple_order: List[int]) -> Any:
        """Renumber atoms."""
        return Chem.RenumberAtoms(mol, full_tuple_order)

    def replace_core(
        self,
        mol: Any,
        core: Any,
        label_by_index: bool = False,
        replace_dummies: bool = True,
        require_dummy_match: bool = False,
    ) -> Any:
        """Replace core."""
        return Chem.ReplaceCore(
            mol,
            core,
            labelByIndex=label_by_index,
            replaceDummies=replace_dummies,
            requireDummyMatch=require_dummy_match,
        )

    def get_editable_mol(self, mol: Any) -> Any:
        """Get editable molecule."""
        return Chem.RWMol(mol)

    def get_noneditable_mol(self, mol: Any) -> Any:
        """Get non-editable molecule."""
        return mol.GetMol()

    def combine_mols(self, mol1: Any, mol2: Any) -> Any:
        """Combine mols.
        
        NOTE: mol1 may need to be editable. Not certain.
        """
        return Chem.CombineMols(mol1, mol2)

    def get_atoms(self, mol: Any) -> List[Any]:
        """Get atoms from molecule."""
        return mol.GetAtoms()

    def get_atomic_num(self, atom: Any) -> int:
        """Get atomic number."""
        return atom.GetAtomicNum()

    def get_bonds(self, atom: Any) -> List[Any]:
        """Get bonds for atom."""
        return atom.GetBonds()

    def get_bond_type_as_double(self, bond: Any) -> float:
        """Get bond type as double value."""
        return bond.GetBondTypeAsDouble()

    def get_num_atoms(self, mol: Any) -> int:
        """Get number of atoms in molecule."""
        return mol.GetNumAtoms()

    def set_atom_map_num(self, atom: Any, num: int) -> None:
        """Set atom map number."""
        atom.SetAtomMapNum(num)

    def get_isotope(self, atom: Any) -> int:
        """Get isotope."""
        return atom.GetIsotope()

    def set_isotope(self, atom: Any, isotope: int) -> None:
        """Set isotope."""
        atom.SetIsotope(isotope)

    def get_atom_with_idx(self, mol: Any, idx: int) -> Any:
        """Get atom with index."""
        return mol.GetAtomWithIdx(idx)

    def get_idx(self, atom: Any) -> int:
        """Get index."""
        return atom.GetIdx()

    def get_neighbors(self, atom: Any) -> List[Any]:
        """Get neighbors."""
        return atom.GetNeighbors()

    def get_substruct_matches(
        self, mol: Any, query: Any, uniquify: bool = True, max_matches: int = 1000
    ) -> List[List[int]]:
        """Get substructure matches."""
        return mol.GetSubstructMatches(query, uniquify=uniquify, maxMatches=max_matches)

    def get_bond_between_atoms(self, mol: Any, atom1_idx: int, atom2_idx: int) -> Any:
        """Get bond between atoms."""
        return mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)

    def get_bond_type(self, bond: Any) -> Any:
        """Get bond type."""
        return bond.GetBondType()
    
    def is_in_ring(self, atom: Any) -> bool:
        """Check if atom is in ring."""
        return atom.IsInRing()

    def mol_from_pdb_file(self, pdb_file: str, sanitize: bool = True, remove_hs: bool = False) -> Any:
        """Create molecule from PDB file."""
        return Chem.MolFromPDBFile(pdb_file, sanitize=sanitize, removeHs=remove_hs)

    def is_aromatic(self, bond: Any) -> bool:
        """Check if a bond is aromatic."""
        return bond.GetIsAromatic()

    def brics_decompose(
        self,
        mol: Any,
        return_mols: bool = True,
        min_fragment_size: int = 1,
        keep_non_leaf_nodes: bool = False,
    ) -> List[Any]:
        """Decompose molecule using BRICS."""
        return list(
            BRICSDecompose(
                mol,
                returnMols=return_mols,
                minFragmentSize=min_fragment_size,
                keepNonLeafNodes=keep_non_leaf_nodes,
            )
        )

    def fragment_on_bonds(
        self, mol: Any, bond_indices: List[int], add_dummies: bool = True
    ) -> Any:
        """Fragment molecule on specified bonds."""
        return Chem.FragmentOnBonds(mol, bond_indices, addDummies=add_dummies)

    def mols_to_grid_image(
        self,
        mols: List[Any],
        mols_per_row: int,
        sub_img_size: Tuple[int, int],
        highlight_atom_lists: Optional[List[List[int]]] = None,
    ) -> Any:
        """Create a grid image of molecules."""
        return Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            highlightAtomLists=highlight_atom_lists,
        )

    def compute_2d_coords(self, mol: Any) -> int:
        """Compute 2D coordinates for a molecule."""
        return AllChem.Compute2DCoords(mol)

    def get_morgan_fingerprint_as_bit_vect(
        self, mol: Any, radius: int, n_bits: int
    ) -> Any:
        """Generate Morgan fingerprint as a bit vector."""
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

    def mols_from_sdf_file(
        self, sdf_file: str, sanitize: bool = False, remove_hs: bool = True
    ) -> List[Any]:
        """Load molecules from an SDF file."""
        # Note: returns an iterator in RDKit, but we'll return a list for simplicity
        supplier = Chem.SDMolSupplier(sdf_file, sanitize=sanitize, removeHs=remove_hs)
        return [mol for mol in supplier if mol is not None]

    def run_reactants(
        self, reaction: Any, reactants: Tuple[Any, ...]
    ) -> List[List[Any]]:
        """Run a reaction with the given reactants."""
        return reaction.RunReactants(reactants)

    def initialize_reaction(self, reaction: Any):
        """Initialize a reaction object."""
        reaction.Initialize()

    def get_prop(self, mol: Any, prop_name: str) -> Any:
        """Get a property from a molecule."""
        return mol.GetProp(prop_name)

    def mol_to_pdb_block(self, mol: Any) -> str:
        """Convert molecule to a PDB block string."""
        return Chem.MolToPDBBlock(mol)

    def get_rdk_fingerprint(self, mol: Any, max_path: int, fp_size: int) -> Any:
        """Generate RDKit fingerprint."""
        return rdmolops.RDKFingerprint(mol, maxPath=max_path, fpSize=fp_size)

    def bit_vect_to_bit_string(self, bit_vect: Any) -> str:
        """Convert a bit vector to a bit string."""
        return bit_vect.ToBitString()

    def mol_to_sdf_file(self, mol: Any, file_path: str) -> None:
        """Write a molecule to an SDF file."""
        writer = Chem.SDWriter(file_path)
        writer.write(mol)
        writer.close()

    def create_empty_editable_mol(self) -> Any:
        """Create an empty editable molecule."""
        return Chem.EditableMol(Chem.Mol())

    def copy_atom(self, atom: Any) -> Any:
        """Create a copy of an atom."""
        return Chem.Atom(atom)

    def add_atom_to_mol(self, editable_mol: Any, atom: Any) -> int:
        """Add an atom to an editable molecule."""
        return editable_mol.AddAtom(atom)

    def add_bond_to_mol(self, editable_mol: Any, begin_atom_idx: int, end_atom_idx: int, bond_type: Any) -> None:
        """Add a bond to an editable molecule."""
        editable_mol.AddBond(begin_atom_idx, end_atom_idx, bond_type)

    def create_conformer(self, num_atoms: int) -> Any:
        """Create a conformer."""
        return Chem.Conformer(num_atoms)

    def get_conformer(self, mol: Any, conf_id: int = -1) -> Any:
        """Get a conformer from a molecule."""
        return mol.GetConformer(conf_id)

    def get_atom_position(self, conformer: Any, atom_idx: int) -> Any:
        """Get the position of an atom in a conformer."""
        return conformer.GetAtomPosition(atom_idx)

    def set_atom_position_in_conformer(self, conformer: Any, atom_idx: int, pos: Any) -> None:
        """Set the position of an atom in a conformer."""
        conformer.SetAtomPosition(atom_idx, pos)

    def add_conformer_to_mol(self, mol: Any, conformer: Any) -> int:
        """Add a conformer to a molecule."""
        return mol.AddConformer(conformer)

    def get_num_conformers(self, mol: Any) -> int:
        """Get the number of conformers for a molecule."""
        return mol.GetNumConformers()

    def get_atom_degree(self, atom: Any) -> int:
        """Get the degree of an atom."""
        return atom.GetDegree()

    def mol_to_smarts(self, mol: Any) -> str:
        """Convert a molecule to a SMARTS string."""
        return Chem.MolToSmarts(mol)

    def get_atom_property(self, atom: Any, prop: str) -> Any:
        """Get a property from an atom."""
        return atom.GetProp(prop)

    def set_atom_property(self, atom: Any, prop: str, value: Any) -> None:
        """Set a property on an atom."""
        atom.SetProp(prop, value)

    def get_substruct_match(self, mol: Any, query: Any) -> Optional[Tuple[int, ...]]:
        """Get a single substructure match."""
        return mol.GetSubstructMatch(query)

    def remove_atom_from_editable_mol(self, editable_mol: Any, atom_idx: int) -> None:
        """Remove an atom from an editable molecule."""
        editable_mol.RemoveAtom(atom_idx)

    def get_single_bond_type(self) -> Any:
        """Get the single bond type object."""
        return Chem.BondType.SINGLE
    def create_atom(self, atomic_num: int) -> Any:
        """Create an atom."""
        return Chem.Atom(atomic_num)

    def filter_has_match(self, filter_obj: Any, mol: Any) -> bool:
        """Check if a molecule has a match in a filter."""
        return filter_obj.HasMatch(mol)

    def mol_from_sdf_file(self, sdf_file: str) -> Optional[Any]:
        """Create a single molecule from an SDF file."""
        supplier = Chem.SDMolSupplier(sdf_file, sanitize=False, removeHs=True)
        if supplier:
            return next(iter(supplier), None)
        return None

# TO CONSIDER
# HasSubstructMatch
# Could search for "mol."
# RunReactants
# .AddBond
# SetAtomMapNum
# GetIsotope
# GetNeighbors
# GetAtomWithIdx
# SetIsotope
# GetSubstructMatches
