"""
Run file type conversions from .smi to .sdf to .pdb
"""
import __future__

import glob
import sys
import os
from os.path import basename

import rdkit
import rdkit.Chem as Chem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

# Some pathing to allow for importing Gypsum Functions
current_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
gypsum_dir = str(current_dir) + os.sep + "convert_files" + os.sep + "gypsum_dl"
gypsum_gypsum_dir = (
    str(current_dir)
    + os.sep
    + "convert_files"
    + os.sep
    + "gypsum_dl"
    + os.sep
    + "gypsum_dl"
)
sys.path.extend([gypsum_dir, current_dir, gypsum_gypsum_dir])

import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH



class StdoutRedirection:
    """Standard output redirection context manager

    Class taken from
    https://stackoverflow.com/questions/4110891/how-to-redirect-the-output-of-print-to-a-txt-file
    """

    def __init__(self, path):
        """
        Inputs:
        :param str path: the path
        """
        self._path = path

    def __enter__(self):
        """
        This will return the class self object
        and flush the print statements.

        Returns:
        :returns: self self: class self object
        """
        sys.stdout.flush()
        sys.stdout = open(self._path, mode="w")
        sys.stdout.flush()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit the output piping

        Inputs:
        :param obj exc_type: exc_type
        :param obj exc_val: exc_val
        :param obj exc_tb: exc_tb
        """
        sys.stdout.flush()
        sys.stdout.close()
        sys.stdout = sys.__stdout__


def convert_to_3d(vars, smi_file, smile_file_directory):
    """
    This function converts SMILES from 1D to 3D using gypsum Gypsum converts
    SMILES in an .smi file to 3D .sdf files Then rdkit converts the sdfs to
    PDB files.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param str smi_file: the file name of the .smi file
    :param srt smile_file_directory: the directory path which contains the
        .smi file
    """

    print("CONVERTING SMILES TO SDF")
    # convert smiles in an .SMI file to sdfs using gypsum
    gypsum_output_folder_path = convert_smi_to_sdfs_with_gypsum(
        vars, smi_file, smile_file_directory
    )
    print("CONVERTING SMILES TO SDF COMPLETED")

    print("CONVERTING SDF TO PDB")
    # convert sdf files to PDBs using rdkit
    convert_sdf_to_pdbs(vars, smile_file_directory, gypsum_output_folder_path)
    print("CONVERTING SDF TO PDB COMPLETED")


def convert_smi_to_sdfs_with_gypsum(vars, gen_smiles_file, smile_file_directory):
    """
    Convert a file of SMILES to a set of 3d .sdf files using Gypsum. This does
    so by making a set of .json files for running Gypsum for every ligand in
    the .smi file.

    This will print out the list of ligands which failed to convert to 3D.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param str gen_smiles_file: the file name of the .smi file to be converted
        to 3D sdf's
    :param srt smile_file_directory: the directory path which contains the
        .smi file

    Returns:
    :returns: str gypsum_output_folder_path: a path to the folder with all of
        the 3D sdf's created by gypsum.
    """

    max_variance = vars["max_variants_per_compound"]
    number_of_processors = vars["number_of_processors"]
    max_variants_per_compound = vars["max_variants_per_compound"]
    gypsum_thoroughness = vars["gypsum_thoroughness"]
    min_ph = vars["min_ph"]
    max_ph = vars["max_ph"]
    pka_precision = vars["pka_precision"]

    # Make a new folder to put gypsum .smi's and json. Name folder
    # gypsum_submission_files.
    folder_path = "{}gypsum_submission_files{}".format(smile_file_directory, os.sep)
    if os.path.exists(folder_path) is False:
        os.makedirs(folder_path)

    # Make Output for Gypsum folder (where .sdf's go)
    gypsum_output_folder_path = "{}3D_SDFs{}".format(smile_file_directory, os.sep)
    if os.path.exists(gypsum_output_folder_path) is False:
        os.makedirs(gypsum_output_folder_path)

    # Make a folder to put the log files into within the 3D_SDFs folder
    gypsum_log_path = "{}log{}".format(gypsum_output_folder_path, os.sep)
    if os.path.exists(gypsum_log_path) is False:
        os.makedirs(gypsum_log_path)

    # Make All of the json files to submit to gypsum
    list_of_jsons = make_smi_and_gyspum_submit(
        gen_smiles_file,
        folder_path,
        gypsum_output_folder_path,
        max_variance,
        gypsum_thoroughness,
        min_ph,
        max_ph,
        pka_precision,
    )

    timeout_option = vars["timeout_vs_gtimeout"]
    gypsum_timeout_limit = vars["gypsum_timeout_limit"]
    if "python_path" not in vars.keys():
        python_path = "python"
    else:
        if (
                vars["python_path"] != "python"
                and os.path.exists(vars["python_path"]) is False
        ):
            printout = "python path provided was not found."
            printout = printout + "\n\t{}".format(vars["python_path"])
            printout = (
                printout
                + "\nPlease either leave blank or provide proper \
                path to python enviorment."
            )
            print(printout)
            raise Exception(printout)
        else:
            python_path = vars["python_path"]

    # create a the job_inputs to run gypsum in multithread
    job_input = tuple(
        [
            tuple(
                [
                    gypsum_log_path,
                    json_path,
                    timeout_option,
                    gypsum_timeout_limit,
                    python_path,
                ]
            )
            for json_path in list_of_jsons
        ]
    )

    if vars["parallelizer"].return_mode() == "mpi":
        failed_to_convert = vars["parallelizer"].run(
            job_input, run_gypsum_multiprocessing_mpi
        )
        sys.stdout.flush()
    else:
        failed_to_convert = vars["parallelizer"].run(
            job_input, run_gypsum_multiprocessing
        )

    lig_failed_to_convert = [x for x in failed_to_convert if x is not None]
    lig_failed_to_convert = list(set(lig_failed_to_convert))
    if len(lig_failed_to_convert) > 0:
        print("The Following ligands Failed to convert in Gypsum")
        print("Likely due to a Timeout")
        print(lig_failed_to_convert)
    sys.stdout.flush()
    return gypsum_output_folder_path


def make_smi_and_gyspum_submit(gen_smiles_file, folder_path,
                               gypsum_output_folder_path, max_variance,
                               gypsum_thoroughness, min_ph, max_ph,
                               pka_precision):
    """
    Make an individual .json file to submit to Gypsum for every ligand in the
    .smi file. It then executes gypsum for every .json file. This also makes a
    .smi file for each ligand which will be required for Gypsum.

    The .smi file for each ligand will be noted within the .json file to
    instruct Gypsum where to find the SMILE.

    Inputs:
    :param str gen_smiles_file: the file name of the .smi file to be converted
        to 3D sdf's
    :param srt folder_path: the directory path which will contain the inputs
        and outputs from Gypsum
    :param str gypsum_output_folder_path: a path to the folder with all of the
        3D sdf's created by gypsum.
    :param int max_variance: User variable for how many conformers per ligand
        should be made by Gypsum
    :param int gypsum_thoroughness: User variable for How widely Gypsum-DL
        will search for low-energy conformers. Larger values increase run times
        but can produce better results
    :param float min_ph: User variable for Minimum pH to consider by
        Dimorphite-DL
    :param float max_ph: User variable for Maximum pH to consider by
        Dimorphite-DL
    :param float pka_precision: User variable for Size of pH substructure
        ranges by Dimorphite-DL

    Returns:
    :returns: list list_of_jsons: a list of the paths to all .jsons to be
        submitted to Gypsum.
    """
    list_of_jsons = []

    with open(gen_smiles_file) as smiles_file:
        for line in smiles_file:
            line_full = line
            if line == "\n":
                continue
            line = line.replace("\n", "")
            line = line.replace("    ", "\t")
            parts = line.split("\t")  # split line into parts seperated by 4-spaces
            if len(parts) == 0 or len(parts) == 1:
                print(parts)
            smile = parts[0]
            # ligand_name example
            # (Gen_30_Cross_639427+Gen_31_Cross_717928)Gen_34_Cross_709666 But
            # bash doesn't like + or () for file names so we will abridge
            # lig_name_short name for above example becomes
            # Gen_34_Cross_709666 if ligand is from the source files we wont
            # split the name

            ligand_name = parts[1]
            if len(ligand_name.split(")")) == 2:
                lig_name_short = ligand_name.split(")")[1]
            elif len(ligand_name.split(")")) == 1:
                lig_name_short = ligand_name
            else:
                printout = "Ligand name failed to abridge. Smiles may be \
                            named in improper format please seperate with _ \
                            or camelcase. Our formatting is: \
                            (Gen_2_Cross_631+Gen_3_Cross_744)Gen_4_Cross_702 \
                            which reads as Gen_34_Cross_702 (aka ligand 702) \
                            was produced by crossover using ligands: \
                            Gen_2_Cross_631 and Gen_3_Cross_744. \
                            This will abridge to Gen_4_Cross_702 for saving \
                            files.\nThe failed ligand name was \
                            {}".format(ligand_name)

                print(printout)
                raise Exception(printout)

            smi_line = "{}\t{}".format(smile, lig_name_short)

            json_path = "{}{}.json".format(folder_path, lig_name_short)
            smi_path = "{}{}.smi".format(folder_path, lig_name_short)
            gypsum_out_sdf_path = "{}{}.sdf".format(
                gypsum_output_folder_path, lig_name_short
            )

            # make .smi file
            with open(smi_path, "w") as smi_file:
                smi_file.write(smi_line)

            # Make .json file

            json_params = """{{"source": "{}",
            "output_folder": "{}",
            "num_processors": 1,
            "use_durrant_lab_filters": true,
            "job_manager": "serial",
            "separate_output_files":true,
            "max_variants_per_compound": {},
            "thoroughness": {},
            "min_ph": {},
            "max_ph": {},
            "pka_precision": {}
            }}""".format(
                smi_path,
                gypsum_output_folder_path,
                max_variance,
                gypsum_thoroughness,
                min_ph,
                max_ph,
                pka_precision,
            )

            # make .json file
            with open(json_path, "w") as json_file:
                json_file.write(json_params)

            list_of_jsons.append(json_path)

    return list_of_jsons


def run_gypsum_multiprocessing_mpi(gypsum_log_path, json_path, timeout_option,
                                   gypsum_timeout_limit, python_path):
    """
    This converts the a single ligand from a SMILE to a 3D SDF using Gypsum.
    This is used within a multithread.

    This uses a Try statement and a bash script to be able to create a timeout
    if it requires more than 5 minutes to convert a single ligands

    Inputs:
    :param str gypsum_log_path: a path to the folder to place the log files
        produced when running gypsum.
    :param str json_path: the path to the json file which will be run to
        convert to 3D sdf for a single ligand
    :param str timeout_option: this is taken from vars["timeout_vs_gtimeout"].
        This tells the Bash system whether to use "timeout" or "gtimeout".
        gtimeout is used on most MacOS, while timeout is used on most Linux OS.
    :param int gypsum_timeout_limit: this is taken from
        vars["gypsum_timeout_limit"]. It determines the maximum amount of time to
        run Gypsum per ligand
    :param str python_path: Taken from vars["python_path"]. Not used here but
        for symetry with run_gypsum_multiprocessing

    Returns:
    :returns: str lig_id: the name of the ligand if it failed or None if it
        successfully converted to 3D sdf.
    """
    current_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    gypsum_dir = (
        str(current_dir) + os.sep + "convert_files" + os.sep + "gypsum_dl" + os.sep
    )
    gypsum_gypsum_dir = str(gypsum_dir) + os.sep + "gypsum_dl" + os.sep
    sys.path.extend([current_dir, gypsum_dir, gypsum_gypsum_dir])

    from func_timeout import func_timeout, FunctionTimedOut

    from autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Start import (
        prepare_molecules,
    )

    json_file_name = os.path.basename(json_path)
    lig_id = json_file_name.replace(".json", "")
    log_file = "{}{}_log.txt".format(gypsum_log_path, lig_id)

    # The following information in INPUTS will be overriden once in gypsum_dl.
    # The only information that is import is the json_path, which will be used
    # to define all the parameters.
    INPUTS = {
        "json": json_path,
        "source": None,
        "output_folder": None,
        "job_manager": "serial",
        "num_processors": 1,
        "max_variants_per_compound": None,
        "thoroughness": None,
        "separate_output_files": False,
        "add_pdb_output": False,
        "add_html_output": False,
        "min_ph": None,
        "max_ph": None,
        "pka_precision": None,
        "skip_optimize_geometry": False,
        "skip_alternate_ring_conformations": False,
        "skip_adding_hydrogen": False,
        "skip_making_tautomers": False,
        "skip_enumerate_chiral_mol": False,
        "skip_enumerate_double_bonds": False,
        "let_tautomers_change_chirality": False,
        "use_durrant_lab_filters": False,
        "2d_output_only": False,
        "cache_prerun": False,
        "test": False,
    }

    try:
        with StdoutRedirection(log_file):
            # prepare_molecules(INPUTS)
            func_timeout(gypsum_timeout_limit, prepare_molecules, args=(INPUTS,))

        sys.stdout.flush()
    except:
        # This Ligand Timed out
        return lig_id

    # Check if it worked if it failed return lig_id if it works return None
    did_gypsum_complete = check_gypsum_log_did_complete(log_file)
    if did_gypsum_complete is False:
        return lig_id
    elif did_gypsum_complete is None:
        return lig_id
    
    return None


def run_gypsum_multiprocessing(gypsum_log_path, json_path, timeout_option,
                               gypsum_timeout_limit, python_path):
    """
    This converts the a single ligand from a SMILE to a 3D SDF using Gypsum.
    This is used within a multithread.

    This uses a Try statement and a bash script to be able to create a timeout
    if it requires more than 5 minutes to convert a single ligands

    Inputs:
    :param str gypsum_log_path: a path to the folder to place the log files
        produced when running gypsum.
    :param str json_path: the path to the json file which will be run to
        convert to 3D sdf for a single ligand
    :param str timeout_option: this is taken from vars["timeout_vs_gtimeout"].
        This tells the Bash system whether to use "timeout" or "gtimeout".
        gtimeout is used on most MacOS, while timeout is used on most Linux OS.
    :param int gypsum_timeout_limit: this is taken from
        vars["gypsum_timeout_limit"]. It determines the maximum amount of time to
        run Gypsum per ligand
    :param str python_path: Taken from vars["python_path"]. It is the path to
        the python enviorment to use. If not provided it will defer to using just
        python

    Returns:
    :returns: str lig_id: the name of the ligand if it failed or None if it
        successfully converted to 3D sdf.
    """

    json_file_name = os.path.basename(json_path)
    lig_id = json_file_name.replace(".json", "")
    log_file = "{}{}_log.txt".format(gypsum_log_path, lig_id)

    current_dir_executables = os.path.abspath(os.path.dirname(__file__))
    run_single_gypsum_executable = "{}{}run_single_gypsum_{}.sh".format(
        current_dir_executables, os.sep, timeout_option
    )
    gypsum_dl_py_executable = "{}{}gypsum_dl{}run_gypsum_dl.py".format(
        current_dir_executables, os.sep, os.sep
    )

    # run as a .sh bash script to be able to use timeout feature and pipe log.
    # There is a timeout limit for gypsum defined by users defined within
    # run_single_gypsum_executable_timeout but we add an extra 30 second
    # timeout limit wrapped around for security since we are using the
    # os.system function anyway timeout or gtimeout.
    command = (
        timeout_option
        + " "
        + str(gypsum_timeout_limit + 30)
        + " bash {} {} {} {} {} >> {} 2>> {}".format(
            run_single_gypsum_executable,
            gypsum_dl_py_executable,
            json_path,
            gypsum_timeout_limit,
            python_path,
            log_file,
            log_file,
        )
    )
    try:
        os.system(command)
    except:
        return lig_id
    # Check if it worked if it failed return lig_id if it works return None
    did_gypsum_complete = check_gypsum_log_did_complete(log_file)
    if did_gypsum_complete is False:
        return lig_id
    elif did_gypsum_complete is None:
        return lig_id
    else:
        return None


def check_gypsum_log_did_complete(log_file_path):
    """
    This function checks a log_file_path to see if the last line reads
    "TIMEOUT". If it does then gypsum timed out before converting a .smi. If
    it timedout return False. If it completed return True

    Inputs:
    :param str log_file_path: the path to the log from converting a ligand
        with Gypsum.

    Returns:
    :returns: bool bol:  Returns True if the conversion worked (or if it was a
        blank .smi with no name or SMILES). Returns False if it timedout. Returns
        None if it failed to convert due to errors with the SMILES string, .json,
        or .smi.

        -examples would be
                1) if the SMILE was invalid chemically
                2) There was SMILE name but no SMILE string
                3) The SMILES wasn't a SMILE but rather a different data type
                    ie None, False, non-chemical string...
    """

    if os.path.exists(log_file_path) is False:
        return False

    sys.stdout.flush()
    with open(log_file_path) as log:
        data = log.readlines()
    if len(data) == 0:
        # For whatever reason it didn't write to the file.
        return None

    lastline = data[-1]
    lastline = lastline.replace("\n", "")
    lastline = lastline.replace("\t", "")
    lastline = lastline.replace(" ", "")

    if str(lastline) == "TIMEOUT":
        return False
    elif str(lastline) == "FAILED ERRORS WITH THE LIGAND":
        return None
    else:
        return True


def convert_sdf_to_pdbs(vars, gen_folder_path, sdfs_folder_path):
    """
    It will find any .sdf files within the folder_path and convert them to
    .pdb types using rdkit.Chem. It also makes a subfolder to store the pdb
    files if one doesn't already exist in the folder_path.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param str gen_folder_path: Path of the folder for the current generation
    :param str sdfs_folder_path: Path of the folder with all of the 3D .sdf
        files to convert
    """

    files = []

    if os.path.isdir(sdfs_folder_path):
        # so it's a directory, go through the directory and find all the sdf files
        if sdfs_folder_path[-1:] != os.sep:
            sdfs_folder_path = (
                sdfs_folder_path + os.sep
            )  # so add a / to the end of the directory

        files.extend(glob.glob(sdfs_folder_path + "*.sdf"))
        files.extend(glob.glob(sdfs_folder_path + "*.SDF"))
    files = list(set(files))
    if len(files) == 0:
        printout = "\nThere are no sdf's to convert to PDB's. There may be an issue with Gypsum.\n"
        print(printout)
        raise Exception(printout)

    # create a new subfolder if one doesn't already exist. folder will be with
    # the generation and will be titled PDBs pdb_subfolder_path will become
    # the the output folder
    pdb_subfolder_path = gen_folder_path + "PDBs" + os.sep
    if not os.path.isdir(pdb_subfolder_path):
        os.makedirs(pdb_subfolder_path)

    job_inputs = []
    for file_path in files:
        if "params" in file_path:
            continue
        job_inputs.append(tuple([pdb_subfolder_path, file_path]))
    job_inputs = tuple(job_inputs)

    # Convert sdf files to pdbs in multithread
    vars["parallelizer"].run(job_inputs, convert_single_sdf_to_pdb)


def convert_single_sdf_to_pdb(pdb_subfolder_path, sdf_file_path):
    """
    This will convert a given .sdf into seperate .pdb files.

    Inputs:
    :param str pdb_subfolder_path: Path of the folder to place all created pdb
        files
    :param str sdf_file_path: Path of the sdf_file_path to convert to pdb
        files
    """

    if os.path.exists(sdf_file_path) is True:

        file_basename = basename(sdf_file_path)
        file_basename = file_basename.split("__input1")[0]

        file_output_name = "{}{}_".format(pdb_subfolder_path, file_basename)

        try:
            mols = Chem.SDMolSupplier(
                sdf_file_path, sanitize=False, removeHs=False, strictParsing=False
            )
        except:
            mols = None
        try:
            mols_noH = Chem.SDMolSupplier(
                sdf_file_path, sanitize=True, removeHs=True, strictParsing=False
            )
        except:
            mols_noH = [None for x in range(0, len(mols))]

        if mols is None:
            pass
        elif len(mols) == 0:
            pass
        else:
            # if len(mols)==0 gypsum output a blank file by accident
            # if mols is None rdkit couldn't import the sdf
            if len(mols) != 0:
                counter = 0
                for i in range(0, len(mols)):
                    mol = mols[i]
                    # Extra precaution to prevent None's within a set of good
                    # mols
                    if mol is None:
                        continue

                    mol = MOH.check_sanitization(mol)
                    # Filter out any which failed
                    if mol is None:
                        continue

                    # pdb_name indexed to 1
                    pdb_name = "{}_{}.pdb".format(file_output_name, counter + 1)
                    if mol is not None:  # For extra precaution...
                        Chem.MolToPDBFile(mol, pdb_name, flavor=32)
                        # Add header to PDB file with SMILES containing
                        # protanation and stereochem

                        no_H_Smiles = mols_noH[i]
                        if no_H_Smiles is None:
                            no_H_Smiles = Chem.MolToSmiles(mol)

                        if no_H_Smiles is None:
                            print("SMILES was None for: ", pdb_name)
                            printout = "REMARK Final SMILES string: {}\n".format("None")
                        elif type(no_H_Smiles) == str:
                            printout = "REMARK Final SMILES string: {}\n".format(
                                no_H_Smiles
                            )
                        elif type(no_H_Smiles) == type(Chem.MolFromSmiles("C")):
                            printout = "REMARK Final SMILES string: {}\n".format(
                                Chem.MolToSmiles(no_H_Smiles)
                            )

                        with open(pdb_name) as f:
                            printout = printout + f.read()
                        with open(pdb_name, "w") as f:
                            f.write(printout)
                        printout = ""

                    counter = counter + 1
            else:
                pass
