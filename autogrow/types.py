from dataclasses import dataclass
from typing import Any, List, Optional
from rdkit import Chem  # type: ignore
from enum import Enum


class ScoreType(Enum):
    DOCKING = -2
    DIVERSITY = -1


@dataclass
class PostDockedCompound:  # Get new name when you figure out what context this is used in
    smiles: str
    id: str  # Like naphthalene_22
    short_id: str  # Like naphthalene_22
    additional_info: str  # Like naphthalene_22__1
    docking_score: float  # Like -8.439
    diversity_score: Optional[float] = None
    mol: Optional[Chem.Mol] = None
    fp: Optional[Any] = None
    docked_sdf_path: Optional[str] = None

    # fitness_score: float  # Like -8.439
    # diversity_score: Optional[float] = None
    # lig_efficieny: Optional[float] = None

    def to_list(self) -> List[str]:
        return [
            self.smiles,
            self.id,
            self.short_id,
            self.additional_info,
            str(self.docking_score),
            str(self.diversity_score) if self.diversity_score is not None else "",
            self.docked_sdf_path if self.docked_sdf_path is not None else "",
            # str(self.mol) if self.mol is not None else "",
            # str(self.fp) if self.fp is not None else "",
        ]
    
    @staticmethod
    def from_list(lst: List[str]) -> "PostDockedCompound":
        return PostDockedCompound(
            smiles=lst[0],
            id=lst[1],
            short_id=lst[2] if len(lst) > 2 else "",
            additional_info=lst[3] if len(lst) > 3 else "",
            docking_score=float(lst[4]) if len(lst) > 4 else 0.0,
            diversity_score=float(lst[5]) if len(lst) > 5 else None,
            docked_sdf_path=lst[6] if len(lst) > 6 else None,
        )


@dataclass
class PreDockedCompound:
    smiles: str
    name: str
    previous_docking_score: Optional[float] = None
    previous_diversity_score: Optional[float] = None
    sdf_3d_path: Optional[str] = None

    # (ie. ligand name, SMILES string, docking score, diversity score...)
    # ["CCCC"  "zinc123"   1    -0.1  -0.1]

    def copy(self) -> "PreDockedCompound":
        return PreDockedCompound(
            smiles=self.smiles,
            name=self.name,
            previous_docking_score=self.previous_docking_score,
            previous_diversity_score=self.previous_diversity_score,
            sdf_3d_path=self.sdf_3d_path,
        )

    def to_list(self) -> List[str]:
        resp = [self.smiles, self.name]
        if self.previous_docking_score is not None:
            resp.append(str(self.previous_docking_score))
        if self.previous_diversity_score is not None:
            resp.append(str(self.previous_diversity_score))
        if self.sdf_3d_path is not None:
            resp.append(self.sdf_3d_path)
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

    def to_post_docked_compound(self, score: float, docked_sdf: str) -> PostDockedCompound:
        return PostDockedCompound(
            smiles=self.smiles,
            id=self.name,  # Like naphthalene_22
            short_id=self.name,  # Like naphthalene_22
            additional_info="",  # Like naphthalene_22__1
            docking_score=score,  # Like -8.439
            docked_sdf_path=docked_sdf,
        )
