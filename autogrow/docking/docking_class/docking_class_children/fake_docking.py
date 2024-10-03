"""
The child classes from ParentExample
"""
import __future__

import os
import sys
import glob
from typing import Any, Dict, List, Optional, Tuple, Union
import random
import autogrow.docking.delete_failed_mol as Delete
from autogrow.docking.docking_class.docking_class_children.vina_docking import (
    VinaDocking,
)
from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter
import autogrow.docking.ranking.ranking_mol as Ranking
from autogrow.docking.docking_class.parent_dock_class import ParentDocking
import autogrow.docking.scoring.execute_scoring_mol as Scoring
from autogrow.types import PreDockedCompoundInfo, PostDockedCompoundInfo

# Normally plugins should inherit from the ParentDocking class, but for this
# one, we will inherit from the VinaDocking class to keep it as close to that
# one as possible.


class FakeDocking(VinaDocking):
    """
    RUN FAKE DOCKING

    Inputs:
    :param class ParentDocking: Parent docking class to inherit from
    """

    def dock_ligand(self, lig_pdbqt_filename):
        """
        Generate a random docking score for the ligand and create a mock output file.

        Inputs:
        :param str lig_pdbqt_filename: the ligand pdbqt filename
        """
        # Generate a random docking score (lower is better, typically between -12 and 0)
        random_score = round(random.uniform(-12, 0), 2)

        # Create a mock Vina output file
        output_filename = f"{lig_pdbqt_filename}.vina"
        with open(output_filename, "w") as f:
            f.write(f"REMARK VINA RESULT:    {random_score}    0.000    0.000\n")
            f.write("MODEL 1\n")
            f.write("ENDMDL\n")

        print(f"\tGenerated random docking score for: {lig_pdbqt_filename}")
        print(f"\tRandom docking score: {random_score}")
