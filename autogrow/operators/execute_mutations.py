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


from autogrow.operators.mutant_crossover_parent import CompoundGenerator
from autogrow.plugins.plugin_managers import plugin_managers
from autogrow.plugins.mutation import MutationBase, MutationPluginManager
from autogrow.types import PreDockedCompound
from autogrow.utils.logging import LogLevel, log_debug, log_warning


class MutationGenerator(CompoundGenerator):
    """Handles mutation-specific compound generation."""

    def prepare_params(self) -> Dict[str, Any]:
        mutation_plugin_manager = plugin_managers.Mutation
        mutation_plugin_manager.setup_plugins()
        return {"plugin_manager": mutation_plugin_manager}

    def prepare_job_inputs(
        self, compounds: List[PreDockedCompound], num_to_process: int
    ) -> List[Tuple]:
        return [
            (compounds[i % len(compounds)], self.operation_params["plugin_manager"])
            for i in range(num_to_process)
        ]

    def get_parallel_function(self) -> Callable:
        return _run_mutation_for_multithread

    def make_compound_id(self, result: Tuple) -> str:
        _, reaction_id_number, zinc_id_comp_mol, parent_lig_id = result
        parent_lig_id = parent_lig_id.split(")")[-1]
        random_id_num = random.randint(100, 1000000)

        if zinc_id_comp_mol is None:
            return f"({parent_lig_id})Gen_{self.generation_num}_Mutant_{reaction_id_number}_{random_id_num}"
        return f"({parent_lig_id}+{zinc_id_comp_mol})Gen_{self.generation_num}_Mutant_{reaction_id_number}_{random_id_num}"

    def get_operation_name(self) -> str:
        return "mutation"


def _run_mutation_for_multithread(
    predock_cmpd: PreDockedCompound, mutation_obj: MutationPluginManager
) -> Optional[Tuple[str, int, Union[str, None]]]:
    """
    Performs a single mutation operation on a PreDockedCompound.

    This function is designed to be used in a multithreaded context, allowing
    for parallel processing of mutations.

    Args:
        smile (PreDockedCompound): PreDockedCompound of the molecule to mutate.
        mutation_obj (MutationPluginManager): Mutation object to perform the mutation.

    Returns:
        Optional[Tuple[str, int, Union[str, None]]]: Tuple containing the mutated SMILES
        string, reaction ID, and complementary molecule ID (if any), or None if
        the mutation fails.

    Note:
        This function is a wrapper around the mutation object's run method,
        making it suitable for use in multiprocessing contexts.
    """
    resp = mutation_obj.run(predock_cmpd=predock_cmpd)
    return None if resp is None else tuple(list(resp) + [predock_cmpd.name])
