"""
The child classes from ParentExample
"""

import __future__

import contextlib
import os
import subprocess
import datetime
from typing import Any, Dict, Optional, Union

import rdkit.Chem as Chem  # type: ignore

import autogrow.docking.delete_failed_mol as Delete
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH

from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter


class ObabelConversion(ParentPDBQTConverter):
    """
    This is a class to convert ligands from PDB to PDBQT format using
    commandline obabel

    Openbabel citations:
        - N M O'Boyle, M Banck, C A James, C Morley, T Vandermeersch, and G R
          Hutchison. "Open Babel: An open chemical toolbox." J. Cheminf. (2011),
          3, 33. DOI:10.1186/1758-2946-3-33
        - The Open Babel Package, version 2.3.1 http://openbabel.org (accessed
          Oct 2011)

    Inputs:
    :param class ParentPDBQTConverter: Parent PDBQTConverter class to inherit
      from
    """

    def __init__(
        self, params: Optional[Dict[str, Any]] = None, test_boot: bool = True,
    ) -> None:
        """
        get the specifications for Vina from vars load them into the self
        variables we will need and convert the receptor to the proper file
        format (ie pdb-> pdbqt)

        Inputs:
        :param dict params: Dictionary of User variables
        :param bool test_boot: used to initialize class without objects for
            testing purpose
        """

        if not test_boot:

            assert params is not None, "params must be passed to ObabelConversion"

            self.params = params
            self.debug_mode = params["debug_mode"]

            # VINA SPECIFIC VARS
            receptor_file = params["filename_of_receptor"]
            obabel_path = self.params["obabel_path"]
            number_of_processors = params["number_of_processors"]
            # docking_executable = params["docking_executable"]

            ###########################

            # convert Receptor from PDB to PDBQT
            self.convert_receptor_pdb_files_to_pdbqt(
                receptor_file, obabel_path, None, number_of_processors
            )

            self.receptor_pdbqt_file = f"{receptor_file}qt"

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
        :param str run_path: file path of the obabel_path executable
        :param int number_of_processors: number of processors to multithread
        """

        count = 0
        while not os.path.exists(f"{receptor_file}qt"):

            count = count + 1
            if count > 10000:
                print(
                    "ERROR: I've tried 10,000 times to convert the file \""
                    + receptor_file
                    + '" to the PDBQT format. Aborting program...'
                )
                raise ValueError(
                    "ERROR: I've tried 10,000 times to convert the file \""
                    + receptor_file
                    + '" to the PDBQT format. Aborting program...'
                )

            # make sure the receptors have been converted to PDBQT. If not, do
            # the conversion.
            receptors = [receptor_file]
            need_to_covert_receptor_to_pdbqt = [
                filename
                for filename in receptors
                if not os.path.exists(f"{filename}qt")
            ]
            # This should only be 1 receptor. If Autogrow is expanded to
            # handle multiple receptor, one will need Multiprocess this. Fix
            # this to 1 processor, so no overwriting issues, but could be
            # expanded if we ever wanted to do multiple receptors

            # create a file to run the pdbqt
            for i in need_to_covert_receptor_to_pdbqt:

                self.prepare_receptor_multiprocessing(run_path, i)

    def prepare_receptor_multiprocessing(self, obabel_path, mol_filename):
        """
        This prepares the receptor for multiprocessing.

        Inputs:
        :param str obabel_path: file path of the obabel_path executable
        :param str mol_filename:  the file path of the receptor
        """

        # This converts PDB to PDBQT with rigid converter. The -xrp option
        # preserves the atom index and prevents root and BRANCHES in the
        # output file. Unfortunately these may not be the best converted pdbqt
        # files. We recommend using MGLTools for converting the receptor file.
        # Developer note should replace this when better option becomes
        # available
        print("Converting receptor PDB file to PDBQT using obabel")
        command = f"{obabel_path} -ipdb {mol_filename} -opdbqt -xrp > {mol_filename}qt"

        try:
            os.system(command)
        except Exception as e:
            raise Exception("Could not convert receptor with obabel") from e
        if os.path.exists(f"{mol_filename}qt") is False:
            raise Exception("Could not convert receptor with obabel")

        # Run Clean up on pdbqt receptor
        # obabel leaves a ROOT, ENDROOT, TORSDOF comments in file which
        # need to be removed prior to docking
        # We will comment them out with REMARK
        # We will also add in a header
        printout = "REMARK Receptor file prepared using obabel on: "
        printout = printout + str(datetime.datetime.now()) + "\n"
        printout += f"REMARK Filename is: {mol_filename}qt\n"
        printout += f"REMARK Prepared by running: {command}\n"
        with open(f"{mol_filename}qt", "r") as fil_path:
            for line in fil_path:
                if "ROOT" in line:
                    line = f"REMARK {line}"
                if "TORSDOF" in line:
                    line = f"REMARK {line}"
                printout = printout + line
        with open(f"{mol_filename}qt", "w") as fil_path:
            fil_path.write(printout)

    ###################################################
    # Convert the Ligand from PDB to PDBQT DockingModel
    ###################################################
    def convert_ligand_pdb_file_to_pdbqt(self, pdb_file):
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

        smile_name = self.get_smile_name_from_pdb(pdb_file)

        # gypsum makes 1 files labeled params which is not a valid pdb, but is
        # actually a log Do not convert the params files
        if "params" in pdb_file:
            return False, None

        obabel_path = self.params["obabel_path"]
        # if the file already has been converted to a .pbdqt skip this file
        # take all other pdb file names
        if not os.path.exists(f"{pdb_file}qt"):

            # make sure its in proper format
            self.convert_pdb_to_pdbqt_acceptable_format(pdb_file)
            self.prepare_ligand_processing(obabel_path, pdb_file)
            if not os.path.exists(f"{pdb_file}qt"):
                # FILE FAILED TO CONVERT TO PDBQT DELETE PDB AND RETURN FALSE
                if self.debug_mode is False:
                    print(
                        f"PDBQT not generated: Deleting {os.path.basename(pdb_file)}..."
                    )

                    # REMOVED FOR LIGANDS WHICH FAILED TO CONVERT TO PDBQT
                    Delete.delete_all_associated_files(pdb_file)
                    return False, smile_name
                # In debug mode but pdbqt file does not exist
                print(f"PDBQT not generated: {os.path.basename(pdb_file)}...")
                return False, smile_name

        return True, smile_name

    # Convert Ligand from PDB to PDBQT conversion
    def prepare_ligand_processing(self, obabel_path, mol_filename):
        """
        This function will convert a single ligand from PDB to PDBQT using
        obabel. It has 10seconds to successfully convert this. It will try to
        convert the ligand up to 3 times If it fails to do so 3 times, whether
        because it timed out or because obabel failed or because of an obabel
        Glitch, it will stop and the ligand won't be docked.

        It will print the ligand if it fails 3 times. It will also fail if the
        molecule is unable to be imported into rdkit and sanitized. This is
        because obabel is sensitive to issues like atoms replaced with *,
        formatting errors, and improper valences. Because obabel will crash
        with these issues the RDKit check is especially useful to prevent hard
        crashes.

        Inputs:
        :param str obabel_path: file path of the obabel_path executable
        :param str mol_filename:  the file path of the ligand
        """

        params = self.params
        timeout_option = params["timeout_vs_gtimeout"]

        # Check that the PDB is a valid PDB file in RDKIT
        try:
            mol = Chem.MolFromPDBFile(mol_filename, sanitize=False, removeHs=False)
            if mol is not None:
                mol = MOH.check_sanitization(mol)
        except Exception:
            mol = None

        temp_file = f"{mol_filename}_temp"
        if mol is not None:
            count = 0
            # timeout or gtimeout

            command = f"{timeout_option} 10 {obabel_path} -ipdb {mol_filename} -opdbqt > {mol_filename}qt"

            while not os.path.exists(f"{mol_filename}qt"):
                if count < 3:
                    # We will try up to 3 times
                    try:
                        subprocess.check_output(f"{command} 2> {temp_file}", shell=True)
                    except Exception:
                        with contextlib.suppress(Exception):
                            os.system(f"{command} 2> {temp_file}")
                        if not os.path.exists(f"{mol_filename}qt"):
                            printout = (
                                f"Failed to convert {count} times: {mol_filename}"
                            )
                            print(printout)

                    count = count + 1

                else:
                    printout = f"COMPLETELY FAILED TO CONVERT: {mol_filename}"
                    print(printout)
                    break

        if os.path.exists(temp_file) and self.debug_mode is False:
            command = f"rm {temp_file}"
            try:
                os.system(command)
            except Exception:
                print(f"Check permissions. Could not delete {temp_file}")

    # Convert PDB to acceptable PDBQT file format before converting
    def convert_pdb_to_pdbqt_acceptable_format(self, filename):
        """
        Make sure a PDB file is properly formatted for conversion to pdbqt

        Inputs:
        :param str filename: the file path of the pdb file to be converted
        """
        # read in the file
        output_lines = []
        with open(filename, "r") as f:
            for line in f:
                line = line.replace("\n", "")
                if line[:5] == "ATOM " or line[:7] == "HETATM ":
                    # fix things like atom names with two letters
                    first = line[:11]
                    middle = (
                        line[11:17].upper().strip()
                    )  # Need to remove whitespaces on both ends
                    last = line[17:]

                    middle_firstpart = ""
                    middle_lastpart = middle

                    for _ in range(len(middle_lastpart)):
                        if middle_lastpart[:1].isupper() is not True:
                            break  # you reached the first number

                        middle_firstpart = middle_firstpart + middle_lastpart[:1]
                        middle_lastpart = middle_lastpart[1:]
                    # now if there are more than two letters in
                    # middle_firstpart, keep just two
                    if len(middle_firstpart) > 2:
                        middle_lastpart = middle_firstpart[2:] + middle_lastpart
                        middle_firstpart = middle_firstpart[:2]

                    if middle_firstpart not in ["BR", "ZN", "FE", "MN", "CL", "MG"]:
                        # so just keep the first letter for the element part
                        # of the atom name
                        middle_lastpart = middle_firstpart[1:] + middle_lastpart
                        middle_firstpart = middle_firstpart[:1]

                    middle = middle_firstpart.rjust(3) + middle_lastpart.ljust(3)

                    line = first + middle + last

                    # make sure all parts of the molecule belong to the same
                    # chain and resid
                    line = f"{line[:17]}LIG X 999{line[26:]}"

                output_lines.append(line)
        with open(filename, "w") as f:

            for line in output_lines:
                f.write(line + "\n")
