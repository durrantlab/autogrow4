"""Data structures for managing compounds throughout the AutoGrow docking
pipeline.

Defines classes representing compounds before and after docking, with utilities
for converting between different representations and handling scoring
information.
"""

from dataclasses import dataclass, field
from typing import Any, List, Optional, Type
from rdkit import Chem  # type: ignore
from enum import Enum


class ScoreType(Enum):
    """Enumeration of supported scoring types for compound evaluation.

    Attributes:
        DOCKING: Represents docking scores (typically negative values where
            lower is better).
        DIVERSITY: Represents diversity scores used to maintain population
            variety.
    """

    DOCKING = -2
    DIVERSITY = -1


@dataclass
class Compound:  # Get new id when you figure out what context this is used in
    """Represents a compound that has completed docking and scoring.

    A data container for compounds after docking, storing SMILES, identifiers,
    scores, and paths to related files. Provides methods for serialization
    to/from lists.

    Attributes:
        smiles (str): SMILES representation of the compound.
        id (str): Full identifier (e.g., 'naphthalene_22').
        additional_info (str): Extra information (e.g., 'naphthalene_22__1').
        docking_score (float): Docking score from simulation.
        diversity_score (Optional[float]): Score representing structural
            uniqueness.
        mol (Optional[Chem.Mol]): RDKit molecule object.
        fp (Optional[Any]): Molecular fingerprint.
        sdf_path (Optional[str]): Path to docked structure SDF file.
    """

    smiles: str
    id: str  # Like naphthalene_22
    additional_info: str = ""  # Like naphthalene_22__1
    docking_score: Optional[float]  = None # Like -8.439
    diversity_score: Optional[float] = None
    mol: Optional[Chem.Mol] = None
    fp: Optional[Any] = None
    sdf_path: Optional[str] = None
    history: List[str] = field(default_factory=list)

    # fitness_score: float  # Like -8.439
    # diversity_score: Optional[float] = None
    # lig_efficieny: Optional[float] = None

    @property
    def tsv_line(self) -> str:
        return f"{self.smiles}\t{self.id}\t{self.docking_score}\t{self.diversity_score}\t{self.sdf_path}\t{self.additional_info}\n"
    
    @staticmethod
    def from_tsv_line(tsv_line: str) -> "Compound":
        prts = tsv_line.strip().split()
        cmpd = Compound(
            smiles = prts[0],
            id = prts[1]
        )

        if len(prts) > 2:
            cmpd.docking_score = float(prts[2])
        if len(prts) > 3:
            cmpd.diversity_score = float(prts[3])
        if len(prts) > 4:
            cmpd.sdf_path = prts[4]
        if len(prts) > 5:
            cmpd.additional_info = prts[5]
        return cmpd

    def get_score_by_type(self, score_type: ScoreType) -> float:
        """Retrieves a previous score of the specified type.

        Args:
            score_type (ScoreType): Type of score to retrieve (DOCKING or
                DIVERSITY).

        Returns:
            float: The requested score value.

        Raises:
            ValueError: If the requested score is not available or score type is
                invalid.
        """
        if score_type == ScoreType.DOCKING:
            # NOTE: Used to be associated with index -2
            if self.docking_score is not None:
                return self.docking_score
            raise ValueError("No docking score available")
        if score_type == ScoreType.DIVERSITY:
            # NOTE: Used to be associated with index -1
            if self.diversity_score is not None:
                return self.diversity_score
            raise ValueError("No diversity score available")
        raise ValueError("Invalid score type")
