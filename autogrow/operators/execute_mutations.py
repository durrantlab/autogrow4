"""
Handles mutation operations for molecular structures.

This module contains functions for creating molecular mutations using various
reaction libraries. It includes functionality for parallel processing of
mutations and handling of mutation results.
"""
import __future__

import random
import copy
from typing import Any, Callable, Dict, List, Optional, Tuple, Union


from autogrow.operators.mutant_crossover_parent import (
    CommonParallelResponse,
    CompoundGenerator,
)
from autogrow.plugins.registry_base import plugin_managers
from autogrow.plugins.mutation import MutationBase, MutationPluginManager
from autogrow.types import Compound
from autogrow.utils.logging import LogLevel, log_debug, log_warning


class MutationGenerator(CompoundGenerator):
    """Handles mutation-specific compound generation."""

    def prepare_params(self) -> Dict[str, Any]:
        """
        Prepare operation-specific parameters.
        
        Returns:
            Dict[str, Any]: Dictionary containing the mutation plugin manager.
        """
        mutation_plugin_manager = plugin_managers.Mutation
        mutation_plugin_manager.setup_plugins()
        return {"plugin_manager": mutation_plugin_manager}

    def prepare_job_inputs(
        self, compounds: List[Compound], num_to_process: int
    ) -> List[Tuple]:
        """
        Prepare job inputs for parallel processing.

        Args:
            compounds (List[Compound]): List of compounds to mutate.
            num_to_process (int): Number of compounds to mutate.

        Returns:
            List[Tuple]: List of tuples containing the compound to mutate and the
                mutation plugin manager.
        """
        return [
            (compounds[i % len(compounds)], self.operation_params["plugin_manager"])
            for i in range(num_to_process)
        ]

    def get_parallel_function(self) -> Callable:
        """
        Get the function to run in parallel.
        
        Returns:
            Callable: Function to run in parallel.
        """
        return _run_mutation_for_multithread

    def make_compound_id(self, result: CommonParallelResponse) -> str:
        """
        Generate a unique compound ID.
        
        Args:
            result (CommonParallelResponse): The result of the mutation operation.
            
        Returns:
            str: Unique compound ID.
        """
        # _, reaction_id_number, zinc_id_comp_mol, _, parent_lig_id = result
        # ('[N-]=[N+]=Nc1ccc(N=[N+]=N)c2c(O)cccc12', '22', None, [Compound(smiles='N=[N+]=Nc1ccc(Cl)c2cccc(O)c12', id='naphthalene_43', additional_info='', docking_score=None, diversity_score=None, mol=None, fp=None, sdf_path=None, history=[])], 'naphthalene_43')

        parent_lig_id = result.parent_lig_id.split(")")[-1]
        random_id_num = random.randint(100, 1000000)

        if result.comp_mol_id is None:
            return f"({parent_lig_id})Gen_{self.generation_num}_Mutant_{result.reaction_id}_{random_id_num}"
        return f"({parent_lig_id}+{result.comp_mol_id})Gen_{self.generation_num}_Mutant_{result.reaction_id}_{random_id_num}"

    def get_operation_name(self) -> str:
        """
        Get the name of the operation.
        
        Returns:
            str: The name of the operation.
        """
        return "mutation"

    def get_operation_desc(self, result: CommonParallelResponse) -> str:
        """Get a description of the operation."""
        return f"{result.parent_cmpds[0].smiles} => {result.child_smiles}"

    def get_formatted_respose(self, results: Tuple) -> CommonParallelResponse:
        """Get a formatted response for the operation."""
        return CommonParallelResponse(
            child_smiles=results[0],
            reaction_id=results[1],
            comp_mol_id=results[2],
            parent_cmpds=results[3],
            parent_lig_id=results[4],
        )


def _run_mutation_for_multithread(
    cmpd: Compound, mutation_obj: MutationPluginManager
) -> Optional[Tuple[str, int, Union[str, None]]]:
    """
    Perform a single mutation operation on a Compound.

    This function is designed to be used in a multithreaded context, allowing
    for parallel processing of mutations.

    Args:
        smile (Compound): Compound of the molecule to mutate.
        mutation_obj (MutationPluginManager): Mutation object to perform the mutation.

    Returns:
        Optional[Tuple[str, int, Union[str, None]]]: Tuple containing the mutated SMILES
        string, reaction ID, and complementary molecule ID (if any), or None if
        the mutation fails.

    Note:
        This function is a wrapper around the mutation object's run method,
        making it suitable for use in multiprocessing contexts.
    """
    resp = mutation_obj.run(cmpd=cmpd)
    return None if resp is None else tuple(list(resp) + [cmpd.id])
