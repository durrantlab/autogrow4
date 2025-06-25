"""Plugin implementation of chemistry toolkit interface."""
from typing import Any, Dict, List, Optional, Tuple, Type, cast
from abc import abstractmethod
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase


class ChemToolkitBase(PluginBase):
    """Abstract base class for chemistry toolkit plugins."""

    def validate(self, params: dict):
        """Validate plugin parameters."""
        pass

    def run(self, **kwargs) -> Any:
        """Execute the chemistry toolkit method specified in kwargs."""
        method = kwargs.pop("method")
        return getattr(self, method)(**kwargs)

    @abstractmethod
    def mol_from_smiles(self, smiles: str, sanitize: bool = True) -> Any:
        """Create molecule from SMILES."""
        pass

    @abstractmethod
    def mol_to_smiles(
        self, mol: Any, canonical: bool = True, isomeric_smiles: bool = True
    ) -> str:
        """Convert molecule to SMILES."""
        pass

    @abstractmethod
    def mol_from_smarts(self, smarts: str) -> Any:
        """Create molecule from SMARTS."""
        pass

    @abstractmethod
    def reaction_from_smarts(self, smarts: str) -> Any:
        """Create reaction from SMARTS."""
        pass

    @abstractmethod
    def sanitize_mol(self, mol: Any, catch_errors: bool = False) -> Tuple[Any, Any]:
        """Sanitize a molecule.
        
        Args:
            mol (Any): Molecule to sanitize.
            catch_errors (bool): Whether to catch errors.

        Returns:
            Tuple[Any, Any]: The sanitized molecule and the response. The
                response should be an object with a property "name" that maps to
                a string, e.g., "SANITIZE_NONE". See rdkit documentation.
        """
        # NOTE: Second response must be object with property "name" that maps to
        # a string, e.g., "SANITIZE_NONE". See rdkit documentation.
        pass

    @abstractmethod
    def add_hs(self, mol: Any) -> Any:
        """Add hydrogens to molecule."""
        pass

    @abstractmethod
    def remove_hs(self, mol: Any, sanitize: bool = True) -> Any:
        """Remove hydrogens from molecule."""
        pass

    @abstractmethod
    def find_mcs(
        self,
        mols: List[Any],
        match_valences: bool = False,
        ring_matches_ring_only: bool = False,
        complete_rings_only: bool = False,
        timeout: int = 3600,
    ) -> Any:
        """Find maximum common substructure."""
        pass

    @abstractmethod
    def get_mol_frags(
        self, mol: Any, as_mols: bool = False, sanitize_frags: bool = True
    ) -> List[Any]:
        """Get molecular fragments."""
        pass

    @abstractmethod
    def get_morgan_fingerprint(
        self, mol: Any, radius: int, use_features: bool = False
    ) -> Any:
        """Generate Morgan fingerprints."""
        pass

    @abstractmethod
    def dice_similarity(self, fp1: Any, fp2: Any) -> float:
        """Calculate fingerprint similarity."""
        pass

    @abstractmethod
    def get_atoms(self, mol: Any) -> List[Any]:
        """Get atoms from molecule."""
        pass

    @abstractmethod
    def get_atomic_num(self, atom: Any) -> int:
        """Get atomic number."""
        pass

    @abstractmethod
    def get_formal_charge(self, atom: Any) -> int:
        """Get formal charge."""
        pass

    @abstractmethod
    def set_formal_charge(self, atom: Any, charge: int) -> None:
        """Set formal charge."""
        pass

    @abstractmethod
    def get_bonds(self, atom: Any) -> List[Any]:
        """Get bonds for atom."""
        pass

    @abstractmethod
    def get_bond_type_as_double(self, bond: Any) -> float:
        """Get bond type as double value."""
        pass

    @abstractmethod
    def remove_atoms(self, mol: Any, atoms_to_remove: List[int]) -> Any:
        """Remove atoms from molecule."""
        pass

    @abstractmethod
    def descriptors_exact_mol_wt(self, mol: Any) -> float:
        """Get exact molecular weight."""
        pass

    @abstractmethod
    def molsurf_tpsa(self, mol: Any) -> float:
        """Get topological polar surface area."""
        pass

    @abstractmethod
    def get_pains_a_filter(self) -> Any:
        """Get PAINS-A filter."""
        pass

    @abstractmethod
    def get_pains_b_filter(self) -> Any:
        """Get PAINS-B filter."""
        pass

    @abstractmethod
    def get_pains_c_filter(self) -> Any:
        """Get PAINS-C filter."""
        pass

    @abstractmethod
    def get_pains_filter(self) -> Any:
        """Get PAINS filter."""
        pass

    @abstractmethod
    def get_nih_filter(self) -> Any:
        """Get NIH filter."""
        pass

    @abstractmethod
    def get_brenk_filter(self) -> Any:
        """Get Brenk filter."""
        pass

    @abstractmethod
    def lipinski_num_rotatable_bonds(self, mol: Any) -> int:
        """Get number of rotatable bonds."""
        pass

    @abstractmethod
    def lipinski_heavy_atom_count(self, mol: Any) -> int:
        """Get number of heavy atoms."""
        pass
    @abstractmethod
    def lipinski_num_h_donors(self, mol: Any) -> int:
        """Get number of H donors."""
        pass

    @abstractmethod
    def lipinski_num_h_acceptors(self, mol: Any) -> int:
        """Get number of H acceptors."""
        pass

    @abstractmethod
    def rdmolops_get_sssr(self, mol: Any) -> int:
        """Get number of SSSR rings."""
        pass

    @abstractmethod
    def crippen_mol_log_p(self, mol: Any) -> float:
        """Get Crippen LogP."""
        pass

    @abstractmethod
    def crippen_mol_mr(self, mol: Any) -> float:
        """Get Crippen MR."""
        pass

    @abstractmethod
    def renumber_atoms(self, mol: Any, full_tuple_order: List[int]) -> Any:
        """Renumber atoms."""
        pass

    @abstractmethod
    def replace_core(
        self,
        mol: Any,
        core: Any,
        label_by_index: bool = False,
        replace_dummies: bool = True,
        require_dummy_match: bool = False,
    ) -> Any:
        """Replace core."""
        pass

    @abstractmethod
    def get_editable_mol(self, mol: Any) -> Any:
        """Get editable molecule."""
        pass

    @abstractmethod
    def combine_mols(self, mol1: Any, mol2: Any) -> Any:
        """Combine mols.
        
        NOTE: mol1 may need to be editable. Not certain.
        """
        pass

    @abstractmethod
    def get_noneditable_mol(self, mol: Any) -> Any:
        """Get non-editable molecule."""
        pass

    @abstractmethod
    def get_num_atoms(self, mol: Any) -> int:
        """Get number of atoms in molecule."""
        pass

    @abstractmethod
    def set_atom_map_num(self, atom: Any, num: int) -> None:
        """Set atom map number."""
        pass

    @abstractmethod
    def get_isotope(self, atom: Any) -> int:
        """Get isotope."""
        pass

    @abstractmethod
    def set_isotope(self, atom: Any, isotope: int) -> None:
        """Set isotope."""
        pass

    @abstractmethod
    def get_atom_with_idx(self, mol: Any, idx: int) -> Any:
        """Get atom with index."""
        pass

    @abstractmethod
    def get_neighbors(self, atom: Any) -> List[Any]:
        """Get neighbors."""
        pass

    @abstractmethod
    def get_idx(self, atom: Any) -> int:
        """Get index."""
        pass

    @abstractmethod
    def get_substruct_matches(
        self, mol: Any, query: Any, uniquify: bool = True, max_matches: int = 1000
    ) -> List[List[int]]:
        """Get substructure matches."""
        pass

    @abstractmethod
    def get_bond_between_atoms(self, mol: Any, atom1_idx: int, atom2_idx: int) -> Any:
        """Get bond between atoms."""
        pass

    @abstractmethod
    def get_bond_type(self, bond: Any) -> Any:
        """Get bond type."""
        pass

    @abstractmethod
    def is_in_ring(self, atom: Any) -> bool:
        """Check if atom is in ring."""
        pass


class ChemToolkitPluginManager(PluginManagerBase):
    """Plugin manager for chemistry toolkits."""

    def __init__(self, plugin_base_class: Type[PluginBase]):
        """
        Initialize the plugin manager.
        
        Args:
            plugin_base_class (Type[PluginBase]): The base class for plugins
                managed by this manager.
        """
        super().__init__(plugin_base_class)
        self._toolkit = None

    @property
    def toolkit(self):
        """Lazy load the toolkit."""
        if self._toolkit is None:
            toolkits = self.get_selected_plugins_from_params()
            if not toolkits:
                raise Exception("Must specify a chemistry toolkit!")
            if len(toolkits) > 1:
                raise Exception("Can only use one chemistry toolkit at a time!")
            self._toolkit = cast(ChemToolkitBase, self.plugins[toolkits[0]])
        return self._toolkit

    def on_plugin_manager_setup_done(self):
        """
        Perform any initialization tasks for the plugin manager.

        This method is called once during initialization of the plugin manager.
        Children can overwrite it.
        """
        # Just access the property to trigger loading
        _ = self.toolkit

    def execute(self, **kwargs) -> Any:
        """Execute selected chemistry toolkit method."""
        assert (
            False
        ), "For this plugin, please just access the toolkit directly (self.toolkit)."
