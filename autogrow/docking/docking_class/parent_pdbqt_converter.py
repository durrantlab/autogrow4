"""
This script holds the parent class for file conversion for docking.
This is used as the basis for all file conversion classes.
"""
import __future__
from abc import ABC, abstractmethod
import os
from typing import Any, Dict, Optional, Tuple


class ParentPDBQTConverter(ABC):
    """
    Docking

    Inputs:
    :param class object: a class to initialize on
    """

    def __init__(
        self, params: Optional[Dict[str, Any]] = None, test_boot: bool = True,
    ):
        """
        Require to initialize any pdbqt conversion class.

        Inputs:
        :param dict params: Dictionary of User variables
        :param bool test_boot: used to initialize class without objects for
            testing purpose
        """

        pass

    def get_name(self) -> str:
        """
        Returns the current class name.

        Returns:
        :returns: str self.__class__.__name__: the current class name.
        """

        return self.__class__.__name__

    @abstractmethod
    def convert_receptor_pdb_files_to_pdbqt(
        self,
        receptor_file: str,
        run_path: str,
        aux_script_path: Optional[str],
        number_of_processors: int,
    ) -> None:
        """
        Make sure a PDB file is properly formatted for conversion to pdbqt

        Inputs:
        :param str receptor_file:  the file path of the receptor
        :param str run_path: file path of the script or executable used for
            conversion.
        :param str aux_script_path: the receptor4.py file path from mgl tools.
        :param int number_of_processors: number of processors to multithread
        """

        # raise NotImplementedError(
        #     "convert_receptor_pdb_files_to_pdbqt() not implemented"
        # )

        pass

    @abstractmethod
    def convert_ligand_pdb_file_to_pdbqt(
        self, pdb_file: str
    ) -> Tuple[bool, Optional[str]]:
        """
        Convert the ligands of a given directory from pdb to pdbqt format

        Inputs:
        :param str pdb_file: the file name, a string.
        Returns:
        :returns: bool bool: True if it worked; False if its the gypsum param
            file or if it failed to make PDBQT
        :returns: str smile_name: name of the SMILES string from a pdb file
            None if its the param file
        """

        # raise NotImplementedError("rank_and_save_output_smi() not implemented")
        pass

    #######################################
    # Handle Failed PDBS                  #
    #######################################
    def get_smile_name_from_pdb(self, pdb_file: str) -> str:
        """
        This will return the unique identifier name for the compound

        Inputs:
        :param str pdb_file: pdb file path
        Returns:
        :returns: str line_stripped: the name of the SMILES string
                                with the new lines and COMPND removed
        """
        line_stripped = "unknown"
        if os.path.exists(pdb_file):
            with open(pdb_file, "r") as f:
                for line in f:
                    if "COMPND" in line:
                        # Need to remove whitespaces on both ends
                        line_stripped = line.replace("COMPND", "").strip()
                        # Need to remove whitespaces on both ends
                        line_stripped = line_stripped.replace("\n", "").strip()

            # line_stripped is now the name of the smile for this compound

        return line_stripped
