"""
This script will convert a docked .pdbqt.vina file into seperate .pdb file. 

This is done by splitting up a single .pdbqt.vina into seperate .pdbqt files for each docked pose.
Then it removes a column of the .pdbqt and saves as a .pdb file. 

If variable --top_num_of_poses is not set it will convert all poses.
    If --top_num_of_poses == 1 it will only convert the top docked pose to .pdb
    If --top_num_of_poses == 2 it will only convert the top 2 docked poses to .pdb
    If --top_num_of_poses == 10 but there only 8 poses it will convert the 8 poses and stop
    

If --docking_threshold is not set it will convert all poses to .pdb;
    If --docking_threshold == -10.0 it will only convert poses with docking scores less than 
        or equal to -10.0 (Remember docking scores are better when more negative)
    
--docking_threshold and --top_num_of_poses work as AND type operators. 
    Example:
        --docking_threshold == -11.4 and --top_num_of_poses==5
        It will take the top 5 poses as long as they also have docking scores <=-11.4
        

"""

import os
import sys

def convert_pdbqt_to_pdb(pdbqt_file,pdb_file):

    os.system("cut -c 1-60,70-79 {} > {}".format(pdbqt_file,pdb_file))



#######
if __name__ == "__main__":
    try:
        pdbqt_file=sys.argv[1]
    except:
        print("Need a path to a pdbqt to be converted to pdb")
    if os.path.exists(pdbqt_file)==False:
        printout= "pdbqt file does not exists: {}".format(pdbqt_file)
        print(printout)
        raise Exception(printout)

    if ".pdbqt" not in pdbqt_file and ".PDBQT" not in pdbqt_file:
        printout= "pdbqt file must be a .pdbqt file.: {}".format(pdbqt_file)
        print(printout)
        raise Exception(printout)
    if pdbqt_file[-2].lower()!="q" or pdbqt_file[-1].lower()!="t":
        printout= "pdbqt file must be a .pdbqt file.: {}".format(pdbqt_file)
        print(printout)
        raise Exception(printout)

    try:
        pdb_file=sys.argv[2]
        print("converting {} to .pdb".format(pdbqt_file))
        print("\tThe new file will be {}".format(pdb_file))
    except:
        print("PDB outputfile will be the same as pdbqt file minus qt")
        pdb_file = pdbqt_file.replace(".pdbqt",".pdb")
        pdb_file = pdbqt_file.replace(".PDBQT",".PDB")
        print("\t pdbqt_file = {} \t pdb_file = {}".format(pdbqt_file,pdb_file))
    

    convert_pdbqt_to_pdb(pdbqt_file,pdb_file)
    print("finished")