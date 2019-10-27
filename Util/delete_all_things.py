import glob
import sys
import os
# import raw_input
import support_scripts.Multiprocess as mp

def check_response():
    prompt = "ARE YOU SURE YOU WANT TO DELETE EVERYTHING?"
    ask = str(input(prompt+' (y/n): ')).lower().strip()
    if ask[0] == 'y':
        return True
    if ask[0] == 'n':
        raise Exception("TIME TO TERMINATE")
    else:
        return None


def check_layer(layer_top):
    stuff = glob.glob(layer_top +"/*")

    dirs = [x for x in stuff if os.path.isdir(x)==True]
    files = [x for x in stuff if os.path.isfile(x)==True]
    if len(dirs) == 0:
        return True, dirs, files
    else:
        return False,dirs, files

def del_file(file_name):

    todo = "timeout 1 rm -rf {}".format(file_name)
    try:os.system(todo)
    except: pass
    if os.path.exists(file_name)==True:
        todo = "timeout 1 echo '' > {}".format(file_name)
        try:os.system(todo)
        except: pass
        todo = "timeout 1 rm -rf {}".format(file_name)
        try: os.system(todo)
        except: pass


    
def del_folder(folder):
    if os.path.exists(folder)==True:
        os.system("rm -rf {}".format(folder))

    

def run_find(jobs):
    dirs = []
    files = []
    
    bottom_local = []
    dir_list_of_list = []
    files_list_of_list = []
    results  = mp.MultiThreading(jobs, -1, check_layer)
    for x in results:
        bottom_local.append(x[0])
        dir_list_of_list.extend(x[1])
        files_list_of_list.extend(x[2])
    
    bottom_local = list(set(bottom_local))
    
    for i in dir_list_of_list:
        dirs.append(i)
    for i in files_list_of_list:
        files.append(i)

    for bottom in bottom_local:
        if bottom == False:
            return False, dirs, files


    return True, dirs, files

if __name__ == "__main__":
    top_dir = str(sys.argv[1])
    try:response = str(sys.argv[2])
    except:response =None
    
    while response == None:
        response = check_response()

    list_files = []
    list_folders = [[top_dir]]
    layers_to_check = [top_dir]
    bottom =False
    counter =0 
    while bottom==False and counter <100:
        bottom = True
        new_layers_to_check = []
        print("")
        print("Counter: ", counter)
        print("layers_to_check: ", layers_to_check)
        print("")

        jobs = [tuple([x]) for x in layers_to_check]

        bottom_local, dirs, files = run_find(jobs)

        if bottom_local == False:
            bottom = False

            new_layers_to_check.extend(dirs)
        
        list_files.extend(files)
        layers_to_check = new_layers_to_check
        if layers_to_check != []:
            list_folders.append(layers_to_check)
        new_layers_to_check = []
        counter = counter+1
        if bottom == True:
            print("BOTTOM")
            break

    print("")
    print("")
    print("")
    print("deleting Files: ", len(list_files))
    list_files.sort()
    print("1st file in list =: ",list_files[0])

    jobs = [tuple([x]) for x in list_files]
    mp.MultiThreading(jobs, -1, del_file)

    print("")
    print("")
    print("")
    print("deleting Folders: ")


    for i in range(len(list_folders)-1,-1,-1):

        layer = list_folders[i]
        
        jobs = [tuple([x,]) for x in layer]
        mp.MultiThreading(jobs, -1, del_folder)


    print("\nFINSIHED!!!")
