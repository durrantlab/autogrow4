import json
import os
import glob

test_imports()

def test_imports():
    """
    run imports
    """
    #DEMO:
    print("TRIAL")
    # Can't import anything from run_gypsum_dl because parser is open code note in function
    # from autogrow.Operators.ConvertFiles.gypsum_dl.run_gypsum_dl import print_gypsum_citation
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.dimorphite_dl.dimorphite_dl
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.ChemUtils
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.MolContainer
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.MolObjectHandling
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.MyMol
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Parallelizer
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Start
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Utils
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.ThreeD.Convert2DTo3D
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.ChemUtils
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.MolObjectHandling
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Parallelizer
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.MakeTautomers
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Test.Tester
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.EnumerateDoubleBonds
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.ThreeD.GenerateAlternate3DNonaromaticRingConfs
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.EnumerateChiralMols
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.ThreeD.Convert2DTo3D
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.IO.SaveToSDF
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.AddHydrogens
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Start
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.MolContainer
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.IO.LoadFiles
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.PrepareSmiles
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.ThreeD.PrepareThreeD
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Utils
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.DeSaltOrigSmiles
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.MyMol
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.IO.Web2DOutput
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.IO.ProcessOutput
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.ThreeD.Minimize3D
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.IO.SaveToPDB
    import autogrow.Operators.ConvertFiles.gypsum_dl.gypsum_dl.Steps.SMILES.DurrantLabFilter

    print("SUCESS")




def run_gypsum(json_path,gypsum_log_path):

    """
    run gypsum

    Inputs:
    :param str json_path: the path to the json file which will be run to convert to 3D sdf for a single ligand
    :param str gypsum_log_path: a path to the folder to place the log files produced when running gypsum.
    
    Returns:
    :returns: str lig_id: the name of the ligand if it failed or None if it successfully converted to 3D sdf.
    """ 
    json_file_name = os.path.basename(json_path)
    lig_id = json_file_name.replace(".json","")
    log_file = "{}{}_log.txt".format(gypsum_log_path,lig_id)

    current_dir_executables = os.path.abspath(os.path.dirname(__file__))
    run_single_gypsum_executable  = "{}{}run_single_gypsum.sh".format(current_dir_executables, os.sep)
    gypsum_py_executable  = "{}{}gypsum{}run_gypsum.py".format(current_dir_executables, os.sep, os.sep)

    # run as a .sh bash script to be able to use timeout feature and pipe log
    command = "bash {} {} {} > {}".format(run_single_gypsum_executable, gypsum_py_executable, json_path, log_file)
    try:
        os.system(command)
    except:
        return lig_id




if __name__ == "__main__":
    folder = "/home/jacob/Desktop/Outputfolder/Run_10/generation_0/gypsum_submission_files/"
    gypsum_output_folder_path = "/home/jacob/Desktop/Outputfolder/Run_9/generation_0/test_gypsum/"
    gypsum_log_path = "{}log{}".format(gypsum_output_folder_path, os.sep)

    for file_path in glob.glob(folder+"*.json"):
        print("")
        print("############")
        print(file_path)
        run_gypsum(file_path,gypsum_log_path)
        print("############")
        print("")