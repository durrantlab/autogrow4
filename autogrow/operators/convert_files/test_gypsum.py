import json
import os
import glob

test_imports()


def test_imports():
    """
    run imports
    """
    # DEMO:
    print("TRIAL")
    # Can't import anything from run_gypsum_dl because parser is open code
    # note in function from
    # autogrow.operators.convert_files.gypsum_dl.run_gypsum_dl import
    # print_gypsum_citation
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.dimorphite_dl.dimorphite_dl
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.ChemUtils
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolContainer
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MyMol
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Start
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Utils
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.ThreeD.Convert2DTo3D
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.ChemUtils
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.MakeTautomers
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Test.Tester
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.EnumerateDoubleBonds
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.ThreeD.GenerateAlternate3DNonaromaticRingConfs
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.EnumerateChiralMols
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.ThreeD.Convert2DTo3D
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.IO.SaveToSDF
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.AddHydrogens
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Start
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolContainer
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.IO.LoadFiles
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.PrepareSmiles
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.ThreeD.PrepareThreeD
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Utils
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.DeSaltOrigSmiles
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MyMol
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.IO.Web2DOutput
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.IO.ProcessOutput
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.ThreeD.Minimize3D
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.IO.SaveToPDB
    import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Steps.SMILES.DurrantLabFilter

    print("SUCESS")


def run_gypsum(json_path, gypsum_log_path):
    """
    run gypsum

    Inputs:
    :param str json_path: the path to the json file which will be run to
        convert to 3D sdf for a single ligand
    :param str gypsum_log_path: a path to the folder to place the log files
        produced when running gypsum.

    Returns:
    :returns: str lig_id: the name of the ligand if it failed or None if it
        successfully converted to 3D sdf.
    """

    json_file_name = os.path.basename(json_path)
    lig_id = json_file_name.replace(".json", "")
    log_file = "{}{}_log.txt".format(gypsum_log_path, lig_id)

    current_dir_executables = os.path.abspath(os.path.dirname(__file__))
    run_single_gypsum_executable = "{}{}run_single_gypsum.sh".format(
        current_dir_executables, os.sep
    )
    gypsum_py_executable = "{}{}gypsum{}run_gypsum.py".format(
        current_dir_executables, os.sep, os.sep
    )

    # run as a .sh bash script to be able to use timeout feature and pipe log
    command = "bash {} {} {} > {}".format(
        run_single_gypsum_executable, gypsum_py_executable, json_path, log_file
    )
    try:
        os.system(command)
    except:
        return lig_id
