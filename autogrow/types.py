from dataclasses import dataclass
from typing import Any, List, Optional


@dataclass
class CompoundInfo:
    smiles: str
    name: str
    short_id: Optional[str] = None
    docking_score: Optional[float] = None
    diversity_score: Optional[float] = None
    fitness_score: Optional[float] = None

    # (ie. ligand name, SMILES string, docking score, diversity score...)
    # ["CCCC"  "zinc123"   1    -0.1  -0.1]

    def copy(self) -> "CompoundInfo":
        return CompoundInfo(
            smiles=self.smiles,
            name=self.name,
            short_id=self.short_id,
            docking_score=self.docking_score,
            diversity_score=self.diversity_score,
            fitness_score=self.fitness_score,
        )

    def to_list(self) -> List[str]:
        resp = [self.smiles, self.name]
        if self.docking_score is not None:
            resp.append(str(self.docking_score))
        if self.diversity_score is not None:
            resp.append(str(self.diversity_score))
        if self.fitness_score is not None:
            resp.append(str(self.fitness_score))
        return resp

    @property
    def score(self) -> float:
        if self.docking_score:
            return self.docking_score
        if self.diversity_score:
            return self.diversity_score
        if self.fitness_score:
            return self.fitness_score
        raise ValueError("No score available")

    def score_by_index_lookup(self, index: int) -> float:
        # For compatibility with old code, where this was a list, not a
        # dataclass.

        if index == -2 and self.docking_score is not None:
            return self.docking_score
        if index == -1 and self.diversity_score is not None:
            return self.diversity_score
        # No fitness score?
        raise ValueError("No score available")
