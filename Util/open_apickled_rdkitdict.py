import glob
import os
import sys
import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import support_scripts.Multiprocess as mp
import pickle
import copy



def get_mols_dict_from_pickle(file_path):
    print("")
    print(file_path)
    with open(file_path, 'rb') as handle:
        mols_dict = pickle.load(handle)

    return mols_dict

def write_pickle_to_file(file_path, obj):
    with open(file_path, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
        

def write_final_outputs(out_file, pickle_file):
    mols_dict = get_mols_dict_from_pickle(pickle_file)
    
    printout = ""
    with open(out_file, "w") as f:
        f.write(printout)

    with open(out_file, "a") as f:
        for zinc_id in mols_dict:
            printout = ""
            temp = [mols_dict[zinc_id][0], mols_dict[zinc_id][0]]
            printout = printout + "\t".join(temp) + "\n"
            f.write(printout)
            printout = ""
    mols_dict = None


if __name__ == "__main__":
    file_name = sys.argv[1]
    if file_name != "pass":
        if os.path.isdir(file_name)==False:
            print(file_name)
            print("$$$$$$$$$$")
            mols_dict = get_mols_dict_from_pickle(file_name)
            keys_list = list(mols_dict.keys())
            print(mols_dict)
            mols_dict = None

            print(len(keys_list))
            mols_dict = None
        elif os.path.isdir(file_name)==True:
            folder = file_name
            dict_count={}
            for x in glob.glob(folder+"*"):
                print(x)
                mols_dict = get_mols_dict_from_pickle(x)
                # print(mols_dict)
                keys_list = list(mols_dict.keys())
                # print(mols_dict)
                mols_dict = None

                dict_count[x] = len(keys_list)
                    
            print("###########")
            print(dict_count)

            print("###########")


            less_than_100 = []
            more_than_50k = []
            for x in list(dict_count.keys()):
                if dict_count[x] < 100:

                    less_than_100.append(os.path.basename(x))
                elif dict_count[x] > 50000:
                    more_than_50k.append(os.path.basename(x))
            print("###########")
            print("less_than_100")
            print(less_than_100)
            print("###########")
            print("more_than_50k")
            print(more_than_50k)


        else:
            print("")
            print(file_name)
            print(os.path.isdir(file_name))
            print("")
            raise Exception("File/Folder Doesn't Exist")
    # keys_list = list(mols_dict.keys())
    else:
        # print(len(keys_list))
        # /home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/splitup/
        
        # 
        folder = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/splitup/"
        folder = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/splitup_filtered/"
        folder = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/filtered_props/"
        folder = "/home/jspiegel/DataB/jspiegel/FILTER_FOR_AUTO/Large_filter_one/"
        dict_count={}
        for x in glob.glob(folder+"*"):
            print(x)
            mols_dict = get_mols_dict_from_pickle(x)
            # print(mols_dict)
            keys_list = list(mols_dict.keys())
            # print(mols_dict)
            mols_dict = None

            dict_count[x] = len(keys_list)
            
        print("###########")
        print(dict_count)

        print("###########")


        less_than_1000 = []
        for x in list(dict_count.keys()):
            if dict_count[x] < 2000:

                less_than_1000.append(os.path.basename(x))

        print("less_than_1000")
        print(less_than_1000)