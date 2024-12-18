"""
This script contains the class VINA.

This is used to score Vina type docking such as QuickVina2 and Vina
"""
import __future__

import glob
import os
from typing import Dict, List, Optional, Union

from autogrow.docking.scoring.scoring_classes.parent_scoring_class import ParentScoring
from autogrow.types import Compound, Compound


class VINA(ParentScoring):
    """
    Score a given ligand for its binding affinity using VINA, QuickVina02, etc.

    Inputs: :param class ParentFilter: a parent class to initialize off of.
    """

    def __init__(
        self,
        params: Optional[Dict[str, Union[str, int, float, bool]]] = None,
        smiles_dict: Optional[Dict[str, Compound]] = None,
        test_boot: bool = True,
    ) -> None:
        """
        Initialize the class with the given parameters.

        Inputs:
        :param dict params: Dictionary of User variables
        :param dict smiles_dict: a dict of ligand info of SMILES, IDS, and
            short ID
        :param bool test_boot: used to initialize class without objects for
            testing purpose
        """
        if not test_boot:
            self.params = params

            self.smiles_dict = smiles_dict

    #######################
    # Executed by the Execute_Scoring.py script
    #######################
    def find_files_to_score(self, file_path: str) -> List[str]:
        """
        Find all files of the appropriate file format within the dir.
        
        For this class its .pdbqt.vina files.

        ALL SCORING FUNCTIONS MUST HAVE THIS FUNCTION.

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list list_of_files: list of all files to be scored within
            the dir
        """
        self.file_path = file_path
        return glob.glob(f"{file_path}*.pdbqt.vina")

    def run_rescoring(self, vina_output_file: str) -> str:
        """
        Ignore this function.
        
        It is not applicable but is kept because the other rescoring functions
        require this function.

        Inputs:
        :param str vina_output_file: Path to a vina output file to be rescored

        Returns:
        :returns: str "Not Applicable": Because this doesn't need to rescore
            the docking results
        """
        return "Not Applicable"

    def run_scoring(self, file_path: str) -> Optional[Compound]:
        """
        Get all relevant scoring info and return as a list.

        This is required for all Scoring Functions. Additional manipulations
        may go here but there are none for this script..

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list list_of_lig_data: information about the scored ligand.
            Score is last index (ie. [lig_id_shortname, any_details,
            fitness_score_to_use] )
        """
        return self.get_score_from_a_file(file_path)

    def get_score_from_a_file(self, file_path: str) -> Optional[Compound]:
        """
        Make a list of a ligands information including its docking score.

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list lig_info: a list containing all info from
            self.smiles_dict for a given ligand and the ligands id_name and
            the docking score from the best pose.
        """
        # grab the index of the ligand for the score
        basefile = os.path.basename(file_path)
        basefile_strip = basefile.replace(".pdbqt.vina", "")
        basefile_split = basefile.split("__")
        ligand_short_name = basefile_split[0]

        affinity = None

        with open(file_path, "r") as f:
            for line in f:
                if "REMARK VINA" in line:
                    line_stripped = line.replace("REMARK VINA RESULT:", "").replace(
                        "\n", ""
                    )
                    line_split = line_stripped.split()

                    if affinity is None or affinity > float(line_split[0]):
                        affinity = float(line_split[0])

        if affinity is None:
            # This file lacks a pose to use
            return None

        # Obtain additional file

        lig_info_vals = [ligand_short_name, basefile_strip, affinity]

        return self.merge_smile_info_w_affinity_info(lig_info_vals)

    def merge_smile_info_w_affinity_info(self, lig_info: List) -> Optional[Compound]:
        """
        Get info from self.smiles_dict get and merge that with affinity info.

        This will also replace the original SMILES string with that of the
        SMILES string in the PDB which conservers stereoChem

        Inputs:
        :param list lig_info: list containing [ligand_name, affinity]

        Returns:
        :returns: list ligand_full_info: a list containing all info from
            self.smiles_dict for a given ligand and the ligands id_name and
            the docking score from the best pose. returns None if
            ligand_short_name isn't in the self.smiles_dict which should never
            happen
        """
        ligand_short_name, basefile_strip, affinity = lig_info

        # Get SMILES String of PDB
        pdb_path = self.file_path + basefile_strip + ".pdb"
        if not os.path.exists(pdb_path):
            return None

        new_smiles = None
        with open(pdb_path, "r") as f:
            for line in f:
                if "REMARK Final SMILES string: " in line:
                    new_smiles = line.replace(
                        "REMARK Final SMILES string: ", ""
                    ).replace("\n", "")
                    break

        if new_smiles is None:
            # If the REMARK SECTION IS NOT THERE raise except. Avoid this
            # if possible as rdkit can missinterpret bonds because pdbs
            # dont specify bond types
            raise Exception(f"Could not get SMILES string from PDB file: {pdb_path}")

        assert self.smiles_dict is not None, "smiles_dict is None"
        if ligand_short_name in self.smiles_dict:
            return Compound(
                smiles=new_smiles,
                id=ligand_short_name,
                additional_info=lig_info[0],
                docking_score=lig_info[1],
                # diversity_score=lig_info[2],
            )

            # ligand_full_info = [
            #     new_smiles,
            #     self.smiles_dict[ligand_short_name].name,
            #     *lig_info,
            # ]

            # Convert all to strings
            # ligand_full_info = [str(x) for x in ligand_full_info]

            # return ligand_full_info
        # Return None because lig name isn't in dictionary.
        # This is precautionary to prevent key errors later.
        # This should not occur
        return None
