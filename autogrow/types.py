"""Data structures for managing compounds throughout the AutoGrow docking
pipeline.

Defines classes representing compounds before and after docking, with utilities
for converting between different representations and handling scoring
information.
"""

from dataclasses import dataclass
from typing import Any, List, Optional
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
class PostDockedCompound:  # Get new name when you figure out what context this is used in
    """Represents a compound that has completed docking and scoring.

    A data container for compounds after docking, storing SMILES, identifiers,
    scores, and paths to related files. Provides methods for serialization
    to/from lists.

    Attributes:
        smiles (str): SMILES representation of the compound.
        id (str): Full identifier (e.g., 'naphthalene_22').
        short_id (str): Shortened identifier (e.g., 'naphthalene_22').
        additional_info (str): Extra information (e.g., 'naphthalene_22__1').
        docking_score (float): Docking score from simulation.
        diversity_score (Optional[float]): Score representing structural
            uniqueness.
        mol (Optional[Chem.Mol]): RDKit molecule object.
        fp (Optional[Any]): Molecular fingerprint.
        docked_sdf_path (Optional[str]): Path to docked structure SDF file.
    """

    smiles: str
    id: str  # Like naphthalene_22
    short_id: str  # Like naphthalene_22
    additional_info: str  # Like naphthalene_22__1
    docking_score: float  # Like -8.439
    diversity_score: Optional[float] = None
    mol: Optional[Chem.Mol] = None
    fp: Optional[Any] = None
    sdf_path: Optional[str] = None

    # fitness_score: float  # Like -8.439
    # diversity_score: Optional[float] = None
    # lig_efficieny: Optional[float] = None

    def to_list(self) -> List[str]:
        """Converts compound data to a list representation for serialization.

        Returns:
            List[str]: Compound data as strings in a fixed order: SMILES, ID, 
                short_ID, additional info, docking score, diversity score, and
                SDF path.
        """
        return [
            self.smiles,
            self.id,
            self.short_id,
            self.additional_info,
            str(self.docking_score),
            str(self.diversity_score) if self.diversity_score is not None else "",
            self.sdf_path if self.sdf_path is not None else "",
            # str(self.mol) if self.mol is not None else "",
            # str(self.fp) if self.fp is not None else "",
        ]

    @staticmethod
    def from_list(lst: List[str]) -> "PostDockedCompound":
        """Creates a PostDockedCompound instance from a list of strings.

        Args:
            lst (List[str]): List containing compound data in the order: SMILES,
                ID, short_ID, additional info, docking score, diversity score,
                and SDF path.

        Returns:
            PostDockedCompound: New instance populated with the provided data.
        """
        return PostDockedCompound(
            smiles=lst[0],
            id=lst[1],
            short_id=lst[2] if len(lst) > 2 else "",
            additional_info=lst[3] if len(lst) > 3 else "",
            docking_score=float(lst[4]) if len(lst) > 4 else 0.0,
            diversity_score=float(lst[5]) if len(lst) > 5 else None,
            sdf_path=lst[6] if len(lst) > 6 else None,
        )


@dataclass
class PreDockedCompound:
    """Represents a compound before docking, with optional previous scoring data.

    A data container for compounds prior to docking, storing SMILES, name, and
    optional previous scores and file paths.

    Attributes:
        smiles (str): SMILES representation of the compound.
        name (str): Compound identifier.
        previous_docking_score (Optional[float]): Previous docking score if
            available.
        previous_diversity_score (Optional[float]): Previous diversity score if
            available.
        sdf_3d_path (Optional[str]): Path to 3D structure SDF file.
    """

    smiles: str
    name: str
    docking_score: Optional[float] = None
    diversity_score: Optional[float] = None
    sdf_path: Optional[str] = None

    # (ie. ligand name, SMILES string, docking score, diversity score...)
    # ["CCCC"  "zinc123"   1    -0.1  -0.1]

    def copy(self) -> "PreDockedCompound":
        """Creates a deep copy of the compound.

        Returns:
            PreDockedCompound: New instance with copied attributes.
        """
        return PreDockedCompound(
            smiles=self.smiles,
            name=self.name,
            docking_score=self.docking_score,
            diversity_score=self.diversity_score,
            sdf_path=self.sdf_path,
        )

    def to_list(self) -> List[str]:
        """Converts compound data to a list representation for serialization.

        Returns:
            List[str]: Available compound data as strings, including SMILES,
                name, and any available scores and paths.
        """
        resp = [self.smiles, self.name]
        if self.docking_score is not None:
            resp.append(str(self.docking_score))
        if self.diversity_score is not None:
            resp.append(str(self.diversity_score))
        if self.sdf_path is not None:
            resp.append(self.sdf_path)
        return resp

    def get_previous_score(self, score_type: ScoreType) -> float:
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

    def to_post_docked_compound(
        self, score: float, docked_sdf: str
    ) -> PostDockedCompound:
        """Converts to a PostDockedCompound with new docking results.

        Args:
            score (float): New docking score from simulation.
            docked_sdf (str): Path to new docked structure SDF file.

        Returns:
            PostDockedCompound: New instance with docking results and preserved
                compound identity.
        """
        return PostDockedCompound(
            smiles=self.smiles,
            id=self.name,  # Like naphthalene_22
            short_id=self.name,  # Like naphthalene_22
            additional_info="",  # Like naphthalene_22__1
            docking_score=score,  # Like -8.439
            sdf_path=docked_sdf,
        )
