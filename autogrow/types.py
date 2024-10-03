from dataclasses import dataclass
from typing import Any, List, Optional
from rdkit import Chem  # type: ignore
from enum import Enum


class ScoreType(Enum):
    DOCKING = -2
    DIVERSITY = -1


@dataclass
class PostDockedCompoundInfo:  # Get new name when you figure out what context this is used in
    smiles: str
    id: str  # Like naphthalene_22
    short_id: str  # Like naphthalene_22
    additional_info: str  # Like naphthalene_22__1
    score: float  # Like -8.439
    diversity_score: Optional[float] = None
    mol: Optional[Chem.Mol] = None
    fp: Optional[Any] = None

    # fitness_score: float  # Like -8.439
    # diversity_score: Optional[float] = None
    # lig_efficieny: Optional[float] = None

    def to_list(self) -> List[str]:
        return [
            self.smiles,
            self.id,
            self.short_id,
            self.additional_info,
            str(self.score),
            str(self.diversity_score) if self.diversity_score is not None else "",
            # str(self.mol) if self.mol is not None else "",
            # str(self.fp) if self.fp is not None else "",
        ]


@dataclass
class PreDockedCompoundInfo:
    smiles: str
    name: str
    previous_docking_score: Optional[float] = None
    previous_diversity_score: Optional[float] = None

    # (ie. ligand name, SMILES string, docking score, diversity score...)
    # ["CCCC"  "zinc123"   1    -0.1  -0.1]

    def copy(self) -> "PreDockedCompoundInfo":
        return PreDockedCompoundInfo(
            smiles=self.smiles,
            name=self.name,
            previous_docking_score=self.previous_docking_score,
            previous_diversity_score=self.previous_diversity_score,
        )

    def to_list(self) -> List[str]:
        resp = [self.smiles, self.name]
        if self.previous_docking_score is not None:
            resp.append(str(self.previous_docking_score))
        if self.previous_diversity_score is not None:
            resp.append(str(self.previous_diversity_score))
        return resp

    def get_previous_score(self, score_type: ScoreType) -> float:
        if score_type == ScoreType.DOCKING:
            # NOTE: Used to be associated with index -2
            if self.previous_docking_score is not None:
                return self.previous_docking_score
            raise ValueError("No docking score available")
        if score_type == ScoreType.DIVERSITY:
            # NOTE: Used to be associated with index -1
            if self.previous_diversity_score is not None:
                return self.previous_diversity_score
            raise ValueError("No diversity score available")
        raise ValueError("Invalid score type")

    # @property
    # def score(self) -> float:
    #     if self.docking_score:
    #         return self.docking_score
    #     if self.diversity_score:
    #         return self.diversity_score
    #     raise ValueError("No score available")

    # def score_by_index_lookup(self, index: int) -> float:
    #     # For compatibility with old code, where this was a list, not a
    #     # dataclass.

    #     if index == -2 and self.previous_docking_score is not None:
    #         return self.previous_docking_score
    #     # if index == -1 and self.diversity_score is not None:
    #     #     return self.diversity_score
    #     # No fitness score?
    #     raise ValueError("No score available")
