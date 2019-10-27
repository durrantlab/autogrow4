import glob
import sys
import os

import support_scripts.Multiprocess as mp

def gzip(input_tar_file):

    os.system("gzip {}".format(input_tar_file))
    print("FINISHED gzip {}".format(input_tar_file))


def tar(input_folder, output_file):
    
    os.system("tar -cf {} {}".format(output_file, input_folder))
    print("FINISHED TAR {}".format(input_folder))

def run_comp(outfolder, subfolders_to_individually_comp, num_processors):
    
    job_list = []
    output_files = []
    for x in subfolders_to_individually_comp:
        
        output_file = outfolder + x.split(os.sep)[-2] + ".tar"
        output_files.append(output_file)
        temp = tuple([x, output_file])
        job_list.append(temp)

    print(job_list)

    print("\nRunning tar on: ")
    print("\t Infolders: ", subfolders_to_individually_comp)
    print("\t output_file: ", output_files)
    print("")


    mp.MultiThreading(job_list, num_processors, tar)

    failed_tar = []
    passed_tar = []
    print(type(output_files))
    for x in output_files:
        print(x)
        if os.path.exists(x)==False:
            failed_tar.append(x)
        else:
            temp = tuple([x])
            passed_tar.append(temp)
        
    if len(failed_tar)>0:
        print("the following folders failed to tar:", failed_tar)

    print("\nRunning gzip on: ")
    print(passed_tar)
    print("")
    mp.MultiThreading(passed_tar, num_processors, gzip)

    

def copy_file(infile, outfile):
    os.system("cp {} {}".format(infile, outfile))

    print("FINISHED copying {}".format(infile))


def run_stuff():
    
    infolder =sys.argv[1]
    outfolder = sys.argv[2]

    try:
        skip = sys.argv[3]
    except:
        skip = None

    if os.path.exists(outfolder) == False:
        os.mkdir(outfolder)
    if os.path.exists(infolder) == False:
        raise Exception("infolder given does not exist: "+infolder)

    subfolders_to_individually_comp = glob.glob(infolder+"/*/")

    # COPY OVER ANYTHING IN THE TOP DIR BUT NOT IN A FOLDER
    all_infolder = glob.glob(infolder+"/*")
    non_folders = [x for x in all_infolder if x not in subfolders_to_individually_comp]
    jobs = []
    for i in non_folders:
        outfile = outfolder + os.path.basename(i)
        temp = tuple([i, outfile])
        jobs.append(temp)

    print("Copying files in top directory not in folder to final folder")
    print(jobs)
    mp.MultiThreading(jobs, len(jobs), copy_file)
    print("FINISHED CP")


    if skip != None:
        print("HELLLLLLO")
        raise Exception("FINISH")

    num_processors = len(subfolders_to_individually_comp)


    run_comp(outfolder, subfolders_to_individually_comp, num_processors)
    
    
if __name__ == "__main__":
    run_stuff()