"""
This script is used to decompress or recompress AutoGrow data.

If you use the reduce_files_sizes option AutoGrow will convert concatinate and compress
all files in the PDBs directory of each generation. This is useful when doing larger runs as
data transfer is faster and data storage is reduced when files are merged and compressed.
    -The concatination script that is run in AutoGrow 4 can be found at:
            autogrow4/autogrow/docking/concatinate_files.py
This script will either:
    1) Return the files back to their original uncompressed and deconcatinated formating 
                or 
    2) Concatinate and then compress the files into a single file.


The formating of the concatination is:
    "\n##############################File_name: {}\n".format(os.path.basename(file_name_1))
    ... Content of the 1st file...    
    "\n##############################$$END_FILE$$ {}".format(os.path.basename(file_name_1))
    "\n##############################File_name: {}\n".format(os.path.basename(file_name_2))
    ... Content of the 2nd file...    
    "\n##############################$$END_FILE$$ {}".format(os.path.basename(file_name_2))

This concatinated file is tar.gz compressed.
"""

import __future__

import support_scripts.Multiprocess as mp


import glob
import os
import sys
import copy
import gzip
import shutil

def compress_file(file_name):
    """
    Compress the concatinated file
    """
    with open(file_name, "r") as f:
        printout = f.read()
    printout = printout.encode('utf-8')
    with gzip.open(file_name + ".gz", 'wb') as f:
        f.write(printout)

#######
def decompress_file(decompressed_file):
    """
    Decompress a file
    """
    out_file = decompressed_file.replace('.gz',"")
    with gzip.open(decompressed_file, 'rb') as f_comp:
        with open(out_file, 'wb') as f_decomp:
            shutil.copyfileobj(f_comp, f_decomp)
    return out_file

#######
def seperate_files(compressed_file, outfolder):
    """
    This will seperate a compressed and concatinated file into seperated decompressed files.

    """
    directory = os.path.abspath(compressed_file.split(os.path.basename(compressed_file))[0]) + os.sep

    compressed_file = os.path.abspath(compressed_file)
    if compressed_file.split(".")[-1] == "gz":
        decompressed_file = decompress_file(compressed_file)
    else:
        decompressed_file= compressed_file
    if os.path.exists(decompressed_file) == False:
        raise Exception("Failed to decompress the file")

    printout = ""
    list_of_new_files = []
    out_file = None    
    with open(decompressed_file, "r") as f:
        for line in f.readlines():
            if "$$END_FILE$$" in line:
                if out_file != None and os.path.exists(out_file)==False:
                    with open(out_file, "w") as f:
                        f.write(printout + "\n")
                out_file = None
                printout = ""
                continue
            elif "File_name:" in line:

                printout = ""

                # Split the line up and grab the relative file path
                # convert to absolute path
                out_file = outfolder + os.sep + line.split("##############################File_name: ")[1].replace("\n","")
                out_file = os.path.abspath(out_file)
                list_of_new_files.append(out_file)
                continue
            else:
                printout = printout + line
                continue
    

    all_are_made = True
    for f in list_of_new_files:
        if os.path.exists(f)==False:
            print("file failed to decompress: {}".format(f))
            all_are_made = False
    if all_are_made == True:
        torun = "rm {}".format(decompressed_file)
        os.system(torun)
#######
def get_file_info(file_name):
    file_name_insert = "\n##############################File_name: {}\n".format(os.path.basename(file_name))
    file_termination_insert = "\n##############################$$END_FILE$$ {}".format(os.path.basename(file_name))
    concat = file_name_insert + open(file_name).read() + file_termination_insert
    return concat
#######
def del_files(file_or_folder):
    """
    This deletes all temporary files.

    Inputs:
    :param str file_or_folder: the file or folder to delete
    
    """
    if os.path.exists(file_or_folder) == True:
        if os.path.isdir(file_or_folder) == True:
            try:
                shutil.rmtree(file_or_folder)
            except:
                pass
        else:
            try:
                os.remove(command)
            except:
                pass

        # If it failed to delete try via bash command
        if os.path.exists(file_or_folder) == True:
            command = "rm -rf {}".format(file_or_folder)
            try:
                os.system(command)
            except:
                pass
    else:
        pass
    if os.path.exists(file_or_folder):
        print("couldn't delete file: {}".format(file_name))
#######
def run_concatination(directory):
    concat_file = directory + os.sep + "compresed_PDBS.txt"
    print("Start Concatination: To seperate files use the file_concatination_and_compression.py in the Utility script folder.")
    file_list = glob.glob(directory + os.sep + "*")
    file_list = [os.path.abspath(x) for x in file_list]

    with open(concat_file, "a+") as f:
        for file_name in file_list:
            f.write(get_file_info(file_name))

    job_list = [(file_path,) for file_path in file_list]
    print("\tFinish Concatination")
    print("\tRemoving files that were concatinated")
    mp.multi_threading(job_list, -1, del_files)
    print("\tCompressing file")
    compress_file(concat_file)
    if os.path.exists(concat_file  + ".gz"):
        del_files(concat_file)
    print("Finished Compression")
#######
if __name__ == "__main__":
    try:
        todo = sys.argv[1] 
        directory = sys.argv[2]
    except:
        print("Use this script to concatinate the files in a folder into a single outputfolder\
        or deconcatinate files into their original format. This is useful for transfering files\
        more efficiently. example use:\
                1) python concatinate_file.py concat $PATH/Dir/ \
                2) python concatinate_file.py concat $PATH/Dir/")

    directory = os.path.abspath(directory)
    lig_list = glob.glob(directory +os.sep+ "*")
    lig_list = [os.path.abspath(x) for x in lig_list]

    print("BEFORE")
    print(os.path.getsize(directory))


    concat_file = directory + os.sep + "PDB_concated.txt"


    if sys.argv[1] == "concat":
        run_concatination(directory)
        print("FINISH CONCAT")
        print("After concate")
        print(os.path.getsize(directory))
    elif sys.argv[1] == "deconcat":
        compressed_file = sys.argv[2]
        if os.path.exists(compressed_file)==False:
            raise Exception("File to Decompress doesn't exist")
        try:
            outfolder = sys.argv[3]
        except:
            outfolder = os.path.abspath(compressed_file.split(os.path.basename(compressed_file))[0]) + os.sep
        if os.path.exists(outfolder)==False:
            os.mkdir(outfolder)

        directory = os.path.abspath(compressed_file.split(os.path.basename(compressed_file))[0]) + os.sep

        print("BEFORE")
        print(os.path.getsize(directory))

        seperate_files(compressed_file,outfolder)
        print("After deconcate")
        print(os.path.getsize(directory))

        del_files(compressed_file)
        print("After deconcate")
        print(os.path.getsize(directory))
