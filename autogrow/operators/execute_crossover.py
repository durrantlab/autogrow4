"""
Performs crossover operations on molecular structures using SMILES.

This module implements a crossover algorithm for molecular structures, using the
Most Common Substructure (MCS) method to identify suitable pairs of molecules
for crossover. It includes functions for selecting molecules, finding MCS, and
performing the actual crossover operation. TODO: Shouldn't be MCS specific here.
"""
import __future__

import random
import copy
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, cast

from autogrow.operators.mutant_crossover_parent import (
    CommonParallelResponse,
    CompoundGenerator,
)
from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.crossover import CrossoverPluginManager
from autogrow.types import Compound
from autogrow.utils.logging import LogLevel, log_debug, log_warning
import rdkit  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import rdFMCS  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


import autogrow.utils.mol_object_handling as MOH

# TODO: Lots of this code is specific to the MCS crossover. But that code should
# all be in the appropriate plugin, not here in this generic function.


def _test_for_mcs(
    params: Dict[str, Any], mol_1: rdkit.Chem.rdchem.Mol, mol_2: rdkit.Chem.rdchem.Mol
) -> Optional[rdkit.Chem.rdFMCS.MCSResult]:
    """
    Finds the Most Common Substructure (MCS) between two molecules.

    Args:
        params (Dict[str, Any]): User parameters governing the MCS search.
        mol_1 (rdkit.Chem.rdchem.Mol): The first RDKit molecule.
        mol_2 (rdkit.Chem.rdchem.Mol): The second RDKit molecule.

    Returns:
        Optional[rdkit.Chem.rdFMCS.MCSResult]: MCS result object if found,
        None otherwise.

    Note:
        - Recommended to use with molecules that have H's removed.
        - Implicit H's are recognized as part of MCS.
    """
    mols = [mol_1, mol_2]
    time_timeout = params["max_time_mcs_prescreen"]
    min_number_atoms_matched = params["min_atom_match_mcs"]

    try:
        result = rdFMCS.FindMCS(
            mols,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            timeout=time_timeout,
        )
    except Exception:
        return None

    # could be used for a theoretical timeout prefilter was to be implement
    # (but this isn't implemented) (ie. if it times out the prefilter dont use
    # in thorough MCS ligmerge) canceled: if True, the MCS calculation did not
    # finish

    # filter by UserDefined minimum number of atoms found. The higher the
    # number the more similar 2 ligands are but the more restrictive for
    # finding mergable ligands number of atoms in common found
    if result.numAtoms < min_number_atoms_matched:
        return None

    return None if result.canceled else result


def _find_sufficiently_similar_cmpd(
    params: Dict[str, Any], predock_cmpds: List[Compound], query_predock_cmpd: Compound,
) -> Optional[Compound]:
    """
    Selects a random molecule with satisfactory MCS to the given ligand.

    Args:
        params (Dict[str, Any]): User parameters governing the selection.
        predock_cmpds (List[Compound]): List of ligands to choose from.
        query_predock_cmpd (Compound): The reference (query) ligand.

    Returns:
        Optional[Compound]: A suitable second ligand if found,
        None otherwise.
    """
    count = 0
    shuffled_num_list = list(range(len(predock_cmpds) - 1))
    random.shuffle(shuffled_num_list)

    # Convert lig1 into an RDkit mol
    lig1_string = query_predock_cmpd.smiles
    lig1_mol = _convert_mol_from_smiles(lig1_string)

    while count < len(predock_cmpds) - 1:
        rand_num = shuffled_num_list[count]
        mol2_pair = predock_cmpds[rand_num]

        if mol2_pair.smiles == lig1_string:
            count += 1
            continue

        # Convert lig1 into an RDkit mol
        lig_2_string = mol2_pair.smiles
        lig2_mol = _convert_mol_from_smiles(lig_2_string)

        if lig2_mol is None:
            count += 1
            continue

        # it converts and it is not Ligand1. now lets test for a common
        # substructure
        if _test_for_mcs(params, lig1_mol, lig2_mol) is None:
            count += 1
            continue

        # We found a good pair of Ligands
        return mol2_pair

    return None


def _convert_mol_from_smiles(smiles: str) -> Union[rdkit.Chem.rdchem.Mol, bool, None]:
    """
    Converts a SMILES string to an RDKit molecule object.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        Union[rdkit.Chem.rdchem.Mol, bool, None]: RDKit molecule object if 
        conversion is successful, False otherwise.

    Note:
        The function also sanitizes and deprotonates the molecule.
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    except Exception:
        return None

    mol = MOH.check_sanitization(mol)
    if mol is None:
        return None

    mol = MOH.try_deprotanation(mol)
    return False if mol is None else mol


class CrossoverGenerator(CompoundGenerator):
    """Handles crossover-specific compound generation."""

    def prepare_params(self) -> Dict[str, Any]:
        return {k: v for k, v in self.params.items() if k != "parallelizer"}

    def prepare_job_inputs(
        self, compounds: List[Compound], num_to_process: int
    ) -> List[Tuple]:
        # Create a working copy for popping compounds
        working_compounds = compounds.copy()
        # Keep the original list intact for pairing
        available_compounds = compounds.copy()

        pairs = [
            working_compounds.pop() for _ in range(num_to_process) if working_compounds
        ]

        # Use the full available_compounds list for each pair
        return [
            (self.operation_params, compound, available_compounds) for compound in pairs
        ]

    def get_parallel_function(self) -> Callable:
        return _do_crossovers_smiles_merge

    def make_compound_id(self, result: CommonParallelResponse) -> str:
        parent1_id = result.parent_cmpds[0].id.split(")")[-1]
        parent2_id = result.parent_cmpds[1].id.split(")")[-1]
        random_id_num = random.randint(100, 1000000)
        return f"({parent1_id}+{parent2_id})Gen_{self.generation_num}_Cross_{random_id_num}"

    def get_operation_name(self) -> str:
        return "crossover"

    def get_operation_desc(self, result: CommonParallelResponse) -> str:
        """Get a description of the operation."""
        return f"{result.parent_cmpds[0].smiles} + {result.parent_cmpds[1].smiles} => {result.child_smiles} ({self.get_operation_name()})"

    def get_formatted_respose(self, results: Tuple) -> CommonParallelResponse:
        """Get a formatted response for the operation."""
        # (
        #     '[N-]=[N+]=NCCOc1c(C#CC(=O)[O-])ccc2ccccc12',
        #     Compound(smiles='[N-]=[N+]=NCCOc1cccc2ccccc12', id='naphthalene_112', additional_info='', docking_score=None, diversity_score=None, mol=None, fp=None, sdf_path=None, history=[]),
        #     Compound(smiles='O=C([O-])C#Cc1ccc2ccccc2c1', id='naphthalene_85', additional_info='', docking_score=None, diversity_score=None, mol=None, fp=None, sdf_path=None, history=[]))

        return CommonParallelResponse(
            child_smiles=results[0], parent_cmpds=[results[1], results[2]]
        )


# def _find_similar_cmpd(
#     params: Dict[str, Any],
#     ligands_list: List[Compound],
#     lig1_smile_pair: Compound,
# ) -> Optional[Compound]:
#     """
#     Finds a molecule with sufficient shared structure to the given ligand.

#     Args:
#         params (Dict[str, Any]): User parameters governing the process.
#         ligands_list (List[Compound]): List of ligands to choose from.
#         ligand1_pair (Compound): The reference ligand.

#     Returns:
#         Optional[Compound]: A suitable second ligand if found, None otherwise.
#     """
#     ligand_1_string = lig1_smile_pair.smiles

#     # check if ligand_1 can be converted to an rdkit mol
#     lig1 = _convert_mol_from_smiles(ligand_1_string)
#     if lig1 is False:
#         # Ligand1_string failed to be converted to rdkit mol format
#         return None

#     # GET TWO UNIQUE LIGANDS TO WITH A SHARED SUBSTRUCTURE
#     return _find_sufficiently_similar_cmpd(params, ligands_list, lig1_smile_pair)


def _do_crossovers_smiles_merge(
    params: Dict[str, Any],
    lig1_predock_cmpd: Compound,
    all_predock_cmpds: List[Compound],
) -> Optional[Tuple[str, Compound, Compound]]:
    """
    Performs a crossover operation between two ligands.

    Args:
        params (Dict[str, Any]): User parameters governing the process.
        lig1_predock_cmpd (Compound): Information for the first ligand.
        all_predock_cmpds (List[Compound]): List of all seed ligands.

    Returns:
        Optional[Tuple[str, Compound, Compound]]: 
        Tuple containing new ligand SMILES and parent ligand information,
        or None if crossover fails.

    Note:
        Attempts crossover up to 3 times before giving up.
    """
    # Run the run_smiles_merge_prescreen of the ligand. This gets a new a lig2
    # which passed the prescreen.
    # check if ligand_1 can be converted to an rdkit mol
    lig1_rdkit_mol = _convert_mol_from_smiles(lig1_predock_cmpd.smiles)
    if lig1_rdkit_mol is False:
        # Ligand1_string failed to be converted to rdkit mol format
        lig2_predock_cmpd = None

    # GET A UNIQUE LIGANDS TO WITH A SHARED SUBSTRUCTURE
    lig2_predock_cmpd = _find_sufficiently_similar_cmpd(
        params, all_predock_cmpds, lig1_predock_cmpd
    )

    if lig2_predock_cmpd is None:
        return None

    crossover_manager = cast(CrossoverPluginManager, plugin_managers.Crossover)

    counter = 0
    while counter < 3:
        # run SmilesMerge
        ligand_new_smiles = crossover_manager.run(
            predock_cmpd1=lig1_predock_cmpd, predock_cmpd2=lig2_predock_cmpd
        )

        if ligand_new_smiles is None:
            counter += 1
        else:
            # TODO: Filter accepts a list of PreDockedCompounds. So we need to
            # convert smiles string to that just for the purpose of filtering.
            # This is because crossover doesn't return a PreDockCompound object
            # yet. Need to refactor so that happens. As is, smiles getting
            # converted to PreDockCompound twice.
            tmp_predock_cmpd = Compound(smiles=ligand_new_smiles, id="tmp")

            # Filter Here
            pass_or_not = (
                len(plugin_managers.SmilesFilter.run(predock_cmpds=[tmp_predock_cmpd]))
                > 0
            )

            if not pass_or_not:
                counter += 1
            else:
                return (ligand_new_smiles, lig1_predock_cmpd, lig2_predock_cmpd)
    return None
