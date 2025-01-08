"""Module for Maximum Common Substructure (MCS) based crossover in AutoGrow."""

from typing import Any, List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.crossover import CrossoverBase
from autogrow.types import Compound
import autogrow.utils.mol_object_handling as MOH
from autogrow.utils.logging import log_debug
import autogrow.plugins.crossover.merge_functions.alignment_and_breaks as AnB
import autogrow.plugins.crossover.merge_functions.dict_and_r_groups as DnR
import autogrow.plugins.crossover.merge_functions.merge_w_core as MWC
from autogrow.plugins.plugin_manager_instances import plugin_managers


class MergeMCS(CrossoverBase):
    """
    Implements Maximum Common Substructure (MCS) based crossover for molecules.

    This class creates new molecules by merging two existing molecules based on
    their maximum common substructure. It handles molecule sanitization,
    protonation, MCS finding, and R-group selection and merging.
    """

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the plugin name
                and a list of ArgumentVars objects defining the required
                command-line arguments.
        """
        # TODO: These parameter names are not descriptive.

        return (
            "Maximum Common Substructure Crossover",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Run the Maximum Common Substructure Crossover plugin. Creates new molecules by merging two existing molecules.",
                ),
                ArgumentVars(
                    name="max_time_mcs_prescreen",
                    type=int,
                    default=1,
                    help="amount time the pre-screen MCS times out. Time out doesnt prevent \
                    mcs matching just takes what it has up to that point",
                ),
                ArgumentVars(
                    name="max_time_mcs_thorough",
                    type=int,
                    default=1,
                    help="amount time the thorough MCS times out. Time out doesnt prevent \
                    mcs matching just takes what it has up to that point",
                ),
                ArgumentVars(
                    name="min_atom_match_mcs",
                    type=int,
                    default=4,
                    help="Determines the minimum number of atoms in common for a substructurematch. \
                    The higher the more restrictive, but the more likely for two ligands not to match",
                ),
                ArgumentVars(
                    name="protanate_step",
                    action="store_true",
                    default=False,
                    help="Indicates if Smilesmerge uses protanated mols (if true) or deprot \
                    (if False) SmilesMerge is 10x faster when deprotanated",
                ),
            ],
        )

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass

    def run_crossover(
        self, predock_cmpd1: Compound, predock_cmpd2: Compound
    ) -> Optional[str]:
        """
        Run the main script for SmileMerge.

        Args:
            predock_cmpd1 (Compound): Compound of the first
                ligand.
            predock_cmpd2 (Compound):Compound of the second
                ligand.

        Returns:
            Optional[str]: SMILES string for the child ligand derived from mol1
                and mol2. Returns None if it failed at any point.

        Notes:
            Example inputs:
            lig_string_1 = "[N-]=[N+]=NCC(O)COc1cccc2ccccc12"
            lig_string_2 = "C#CCOc1ccc2ccccc2c1CO"
        """
        # lig_string_1 = "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"
        # lig_string_2 = "C# CCOc1ccc2ccccc2c1CO"
        # lig_string_1 = "C1 = CC = CC = C1"
        # Sanitize
        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        mol1 = chemtoolkit.mol_from_smiles(predock_cmpd1.smiles, sanitize=False)
        mol2 = chemtoolkit.mol_from_smiles(predock_cmpd2.smiles, sanitize=False)

        # Sanitize, deprotanate, and reprotanate both molecules
        mol1 = MOH.check_sanitization(mol1)
        mol2 = MOH.check_sanitization(mol2)
        if mol1 is None or mol2 is None:
            return None  # NOTE: JDD: Was False

        protanate_step = self.params["protanate_step"]
        mol1 = MOH.handleHs(mol1, protanate_step)
        mol2 = MOH.handleHs(mol2, protanate_step)

        # check that handleHs() worked for both molecules if fail move on to next
        # pair of molecules
        if mol1 is None or mol2 is None:
            return None

        # make a list of the two rdkit.Chem.rdchem.Mol objects
        mols = [mol1, mol2]

        # Use the below mcs_H function for Most Common Substructure searching.
        # This will prevent broken rings.
        mcs_results = chemtoolkit.find_mcs(
            mols,
            match_valences=False,
            ring_matches_ring_only=True,
            complete_rings_only=False,
            timeout=self.params["max_time_mcs_thorough"],
        )

        if mcs_results.canceled is True:
            return None

        # confirm that this meets the minimum number of matching atoms
        if mcs_results.numAtoms < self.params["min_atom_match_mcs"]:
            return None

        ### Convert mcs_res from into usable and referable forms
        mcs_mol = chemtoolkit.mol_from_smarts(mcs_results.smartsString)

        # handle_mcs_align_labeling_and_cyclicbreaks
        mol1, mol2, mcs_mol = AnB.handle_mcs_align_labeling_and_cyclicbreaks(
            mol1, mol2, mcs_mol
        )

        # confirm that this meets the minimum number of matching atoms
        if mol1 is None or mol2 is None or mcs_mol is None:
            return None
        # This is a very big step. This is where all of the atoms are fragmented
        # to determine what are the R-groups options. R-groups are consolidated
        # into B-groups for each molecule individually, the 2 dictionaries of
        # B-groups are consolidated into a master B-group which is fed into the
        # Mapping class which will randomly select a set of non-clashing B-groups,
        # such that no 2 chosen B's will have connections to the same anchors
        # (anchors are atoms in the common core). It returns a list of
        # smiless for the chosen R-groups. It returns None if a step failed.

        rs_chosen_smiles = DnR.handle_dicts_and_select_b_groups(mol1, mol2, mcs_mol)
        if rs_chosen_smiles is None:
            return None

        # MERGE ALL CHOSEN R-GROUPS WITH THE CORE
        ligand_new_mol = MWC.merge_smiles_with_core(rs_chosen_smiles, mcs_mol)
        if ligand_new_mol is None:
            return None

        # Log the merge
        ligand_new_mol_copy = chemtoolkit.get_editable_mol(ligand_new_mol)
        for atom in chemtoolkit.get_atoms(ligand_new_mol_copy):
            chemtoolkit.set_atom_map_num(atom, 0)
        ligand_new_mol_copy = chemtoolkit.get_noneditable_mol(ligand_new_mol_copy)
        clean_smiles = chemtoolkit.mol_to_smiles(
            ligand_new_mol_copy, canonical=True, isomeric_smiles=False
        )
        log_debug(
            f"Merge by MCS: {predock_cmpd1.smiles} + {predock_cmpd2.smiles} => {clean_smiles}"
        )

        # ligand_new_smiles is either a SMILES string if processing works
        # or None if processing fails
        return self._process_ligand_new_mol(ligand_new_mol)

    def _process_ligand_new_mol(self, ligand_new_mol: Any) -> Optional[str]:
        """
        Process the new ligand molecule.

        This function sanitizes the new molecule, removes isotope labels and
        fragments, checks for unassigned atoms, and converts the molecule to a
        SMILES string.

        Args:
            ligand_new_mol (rdkit.Chem.rdchem.Mol): The new ligand molecule to
                process.

        Returns:
            Optional[str]: The SMILES string of the processed ligand molecule,
                or None if processing failed at any step.
        """
        ligand_new_mol = MOH.check_sanitization(ligand_new_mol)
        if ligand_new_mol is None:
            return None

        chemtoolkit = plugin_managers.ChemToolkit.toolkit

        # REMOVE ALL THE ISOTOPES IN THE NEW MOLECULE
        ligand_new_mol_final = MWC.remove_all_isolabels(ligand_new_mol)

        # Remove any fragments incase 1 made it through
        ligand_new_mol_final = MOH.handle_frag_check(ligand_new_mol_final)
        if ligand_new_mol_final is None:
            return None

        # Make sure there are no unassigned atoms which made it through. These are
        # very unlikely but possible
        ligand_new_mol_final = MOH.check_for_unassigned_atom(ligand_new_mol_final)
        if ligand_new_mol_final is None:
            return None

        ligand_new_mol = MOH.check_sanitization(ligand_new_mol_final)
        if ligand_new_mol is None:
            return None

        return chemtoolkit.mol_to_smiles(ligand_new_mol, isomeric_smiles=True)
