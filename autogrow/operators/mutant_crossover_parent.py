from abc import ABC, abstractmethod
import copy
from dataclasses import dataclass
import random
from typing import Callable, Dict, List, Any, Tuple, Set, Optional

from autogrow.types import Compound
from autogrow.utils.logging import LogLevel, log_debug, log_warning


@dataclass
class CommonParallelResponse:
    """Common response for parallel processing operations."""

    child_smiles: str
    parent_cmpds: List[Compound]
    reaction_id: Optional[str] = None
    parent_lig_id: Optional[str] = None
    comp_mol_id: Optional[str] = None


class CompoundGenerator(ABC):
    """Abstract base class for compound generation operations."""

    def __init__(
        self,
        params: Dict[str, Any],
        generation_num: int,
        procs_per_node: int,
        num_compounds: int,
        predock_cmpds: List[Compound],
    ):
        self.params = params
        self.generation_num = generation_num
        self.procs_per_node = procs_per_node
        self.num_compounds = num_compounds
        self.cmpds = predock_cmpds
        self.operation_params = self.prepare_params()

    @abstractmethod
    def prepare_params(self) -> Dict[str, Any]:
        """Prepare operation-specific parameters."""
        pass

    @abstractmethod
    def prepare_job_inputs(
        self, compounds: List[Compound], num_to_process: int
    ) -> List[Tuple]:
        """Prepare job inputs for parallel processing."""
        pass

    @abstractmethod
    def get_parallel_function(self) -> Callable:
        """Get the function to run in parallel."""
        pass

    @abstractmethod
    def make_compound_id(self, result: CommonParallelResponse) -> str:
        """Generate a unique compound ID."""
        pass

    @abstractmethod
    def get_operation_name(self) -> str:
        """Get the name of the operation for logging."""
        pass

    @abstractmethod
    def get_operation_desc(self, result: CommonParallelResponse) -> str:
        """Get a description of the operation."""
        pass

    @abstractmethod
    def get_formatted_respose(self, results: Tuple) -> CommonParallelResponse:
        """Get a formatted response for the operation."""
        pass

    def generate(self) -> List[Compound]:
        """
        Generate new compounds using the specified operation.

        Returns:
            List[Compound]: List of new compounds.

        Note:
            Uses multiprocessing to generate compounds efficiently.
            Attempts to create unique compounds, avoiding duplicates.
        """
        new_cmpds: List[Compound] = []
        log_debug(
            f"Creating new compounds from selected compounds via {self.get_operation_name()}"
        )

        with LogLevel():
            cmpds_queue = copy.deepcopy(self.cmpds)

            smiles_already_generated = set()
            ids_already_generated = set()
            attempts_to_fill_queue = 0

            while len(new_cmpds) < self.num_compounds and attempts_to_fill_queue < 10:
                attempts_to_fill_queue += 1
                random.shuffle(cmpds_queue)

                # Prepare job inputs
                buffer_num = attempts_to_fill_queue * 2
                num_to_process = self.num_compounds - len(new_cmpds) + buffer_num
                num_to_process = max(num_to_process, self.procs_per_node)

                job_input_list = self.prepare_job_inputs(cmpds_queue, num_to_process)

                # Run parallel operation
                results = self.params["parallelizer"].run(
                    tuple(job_input_list), self.get_parallel_function()
                )

                for idx, res in enumerate(results):
                    if res is None:
                        continue

                    result = self.get_formatted_respose(res)

                    if result.child_smiles in smiles_already_generated:
                        continue

                    # Generate unique ID
                    new_lig_id = ""
                    while new_lig_id in ids_already_generated or not new_lig_id:
                        new_lig_id = self.make_compound_id(result)

                    ids_already_generated.add(new_lig_id)
                    smiles_already_generated.add(result.child_smiles)

                    # Create and store new compound
                    desc = self.get_operation_desc(result)

                    # Update history. If there is only one parent, just append.
                    # If there are multiple parents, append lists containing
                    # each lineage.
                    if sum(len(p._history) for p in result.parent_cmpds) == 0:
                        # Parents have no history. So start of a new history.
                        updated_history = []
                    elif len(result.parent_cmpds) == 1:
                        # Likely mutation. Only one parent, so just extend history.
                        updated_history = result.parent_cmpds[0]._history[:]
                    else:
                        # Likely crossover. Multiple parents, so extend each history.
                        updated_history = [
                            p._history[:] for p in result.parent_cmpds
                        ]

                    ligand_info = Compound(
                        smiles=result.child_smiles,
                        id=new_lig_id,
                        _history=updated_history,
                    )
                    ligand_info.add_history(self.get_operation_name().upper(), desc)
                    new_cmpds.append(ligand_info)

                # Replenish queue if needed
                if not cmpds_queue:
                    cmpds_queue = copy.deepcopy(self.cmpds)
                    random.shuffle(cmpds_queue)

            if len(new_cmpds) < self.num_compounds:
                log_warning(
                    f"Only able to create {len(new_cmpds)} of {self.num_compounds} "
                    f"requested {self.get_operation_name()}s."
                )

            return new_cmpds
