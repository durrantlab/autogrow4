"""RDKit implementation of chemistry toolkit plugin."""
from typing import Any, Dict, List, Tuple
from autogrow.config.argparser import ArgumentVars
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

rdkit.RDLogger.DisableLog("rdApp.*")


class RDKitToolkit(ChemToolkitBase):
    """RDKit implementation of chemistry toolkit."""

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to the RDkit toolkit.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing:
                - The name of the argument group ("Chemistry Toolkit")
                - A list with one ArgumentVars object defining the argument to
                  enable the rdkit chemistry toolkit
        """
        return (
            "Chemistry Toolkit",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Fake conversion to 3D SDF. Converts only to 2D for speed. This is for testing purposes only. Use together with FakeDocking.",
                )
            ],
        )

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
        return Chem.Mol(mol)

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
