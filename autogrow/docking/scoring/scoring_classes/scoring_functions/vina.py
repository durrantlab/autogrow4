"""
This script contains the class VINA.
This is used to score Vina type docking such as QuickVina2 and Vina
"""
import __future__

import glob
import os
from typing import Dict, List, Optional, Union

from autogrow.docking.scoring.scoring_classes.parent_scoring_class import ParentScoring


class VINA(ParentScoring):
    """
    This will Score a given ligand for its binding affinity based on VINA or
    QuickVina02 type docking.

    Inputs:
    :param class ParentFilter: a parent class to initialize off of.
    """

    def __init__(
        self,
        params: Optional[Dict[str, Union[str, int, float, bool]]] = None,
        smiles_dict: Optional[Dict[str, List[str]]] = None,
        test_boot: bool = True,
    ) -> None:
        """
        This will take params and a list of smiles.

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
        Find all files of the appropriate file format within the dir. For this
        class its .pdbqt.vina files.

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
        This is not applicable but is kept because the other rescoring
        functions require this function.

        Inputs:
        :param str vina_output_file: Path to a vina output file to be rescored

        Returns:
        :returns: str "Not Applicable": Because this doesn't need to rescore
            the docking results
        """

        return "Not Applicable"

    def run_scoring(self, file_path: str) -> Optional[List[str]]:
        """
        Get all relevant scoring info and return as a list

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

    def get_score_from_a_file(self, file_path: str) -> Optional[List[str]]:
        """
        Make a list of a ligands information including its docking score.

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list lig_info: a list containing all info from
            self.smiles_dict for a given ligand and the ligands short_id_name and
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

        lig_info = [ligand_short_name, basefile_strip, affinity]

        lig_info = self.merge_smile_info_w_affinity_info(lig_info)
        if lig_info is None:
            return None

        lig_info = [str(x) for x in lig_info]

        return lig_info

    def merge_smile_info_w_affinity_info(
        self, lig_info: List[str]
    ) -> Optional[List[str]]:
        """
        From the info in self.smiles_dict get that info and merge that with
        the affinity info

        This will also replace the original SMILES string with that of the
        SMILES string in the PDB which conservers stereoChem

        Inputs:
        :param list lig_info: list containing [ligand_short_name, affinity]

        Returns:
        :returns: list ligand_full_info: a list containing all info from
            self.smiles_dict for a given ligand and the ligands short_id_name and
            the docking score from the best pose. returns None if
            ligand_short_name isn't in the self.smiles_dict which should never
            happen
        """

        ligand_short_name = lig_info[0]

        # Get SMILES String of PDB
        pdb_path = self.file_path + lig_info[1] + ".pdb"
        if not os.path.exists(pdb_path):
            return None

        new_smiles_string = None
        with open(pdb_path, "r") as f:
            for line in f:
                if "REMARK Final SMILES string: " in line:
                    new_smiles_string = line.replace(
                        "REMARK Final SMILES string: ", ""
                    ).replace("\n", "")
                    break
        if new_smiles_string is None:
            # If the REMARK SECTION IS NOT THERE raise except. Avoid this
            # if possible as rdkit can missinterpret bonds because pdbs
            # dont specify bond types
            raise Exception(f"Could not get SMILES string from PDB file: {pdb_path}")

        assert self.smiles_dict is not None, "smiles_dict is None"
        if ligand_short_name in list(self.smiles_dict.keys()):
            ligand_full_info = [
                new_smiles_string,
                self.smiles_dict[ligand_short_name][1],
                *lig_info,
            ]
            # Convert all to strings
            ligand_full_info = [str(x) for x in ligand_full_info]

            return ligand_full_info
        # Return None because lig name isn't in dictionary.
        # This is precautionary to prevent key errors later.
        # This should not occur
        return None
