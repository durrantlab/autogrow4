from typing import List, Optional, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.crossover import CrossoverBase
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH
from autogrow.utils.logging import log_debug
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
from rdkit.Chem import rdFMCS  # type: ignore
import autogrow.plugins.crossover.merge_functions.alignment_and_breaks as AnB
import autogrow.plugins.crossover.merge_functions.dict_and_r_groups as DnR
import autogrow.plugins.crossover.merge_functions.merge_w_core as MWC

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class MergeMCS(CrossoverBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""

        # TODO: These parameter names are not descriptive.

        return (
            "Maximum Common Substructure Crossover",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Run the Fragment Addition Mutation plugin. Creates new molecules by adding fragments to existing molecules, per user-specified reaction libraries.",
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

    def run_crossover(self, lig_string_1: str, lig_string_2: str) -> Optional[str]:
        """
        This runs the main script for SmileMerge.

        Inputs:
        :param dict params: User variables which will govern how the programs runs

        :param str lig_string_1: smile string for lig 1
        :param str lig_string_2: smile string for lig 2. example: lig_string_1 =
            "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"; example: lig_string_2 = "C#
            CCOc1ccc2ccccc2c1CO"

        Returns:
        :returns: str ligand_new_smiles: smile string for the child ligand derived
            from mol1 and mol2. Returns None if it failed at any point.
        """

        # lig_string_1 = "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"
        # lig_string_2 = "C# CCOc1ccc2ccccc2c1CO"
        # lig_string_1 = "C1 = CC = CC = C1"
        # Sanitize
        mol1 = Chem.MolFromSmiles(lig_string_1, sanitize=False)
        mol2 = Chem.MolFromSmiles(lig_string_2, sanitize=False)

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
        mcs_results = rdFMCS.FindMCS(
            mols,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            timeout=self.params["max_time_mcs_thorough"],
        )

        if mcs_results.canceled is True:
            return None

        # confirm that this meets the minimum number of matching atoms
        if mcs_results.numAtoms < self.params["min_atom_match_mcs"]:
            return None

        ### Convert mcs_res from into usable and referable forms
        mcs_mol = Chem.MolFromSmarts(mcs_results.smartsString)

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
        ligand_new_mol_copy = Chem.RWMol(ligand_new_mol)
        for atom in ligand_new_mol_copy.GetAtoms():
            atom.SetAtomMapNum(0)
        ligand_new_mol_copy = Chem.Mol(ligand_new_mol_copy)
        clean_smiles = AllChem.MolToSmiles(ligand_new_mol_copy, canonical=True, isomericSmiles=False)
        log_debug(f"Merge by MCS: {lig_string_1} + {lig_string_2} => {clean_smiles}")

        # ligand_new_smiles is either a SMILES string if processing works
        # or None if processing fails
        return self._process_ligand_new_mol(ligand_new_mol)

    def _process_ligand_new_mol(
        self, ligand_new_mol: rdkit.Chem.rdchem.Mol
    ) -> Optional[str]:
        """
        This function processes the ligand_new_mol.
        It either returns the SMILES string of ligand_new_mol (ligand_new_smiles)
        or None if it failed at any step.

        Inputs:
        :param str lig_string_1: smile string for lig 1

        Returns:
        :returns: str ligand_new_smiles: either returns the SMILES
            string of ligand_new_mol or None if it failed at any step.
        """

        ligand_new_mol = MOH.check_sanitization(ligand_new_mol)
        if ligand_new_mol is None:
            return None

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

        return Chem.MolToSmiles(ligand_new_mol, isomericSmiles=True)
