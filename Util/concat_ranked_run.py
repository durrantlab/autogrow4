
import glob
import sys
import os

infolder = "/home/jspiegel/projects/AUTOGROW_PARP/Run_1/"
all_folders_list = [f for f in sorted(os.listdir(infolder)) if os.path.isdir(infolder+f)]

folder_list = []
for folder in all_folders_list:
    if folder!="Data" and len(folder.split('_'))==2:
        folder_list.append(folder)

folder_list.sort(key=lambda x: int(x.split('_')[1]))

def concat_ranked_smis(infolder, folder_list):
    concat_file = "/home/jspiegel/projects/AUTOGROW_PARP/Run_1/concat_ranked.smi"
    data_all = []   
    data_dict = {}
    for gen_folder in folder_list:
        gen_folder_name = infolder + gen_folder + "/"
        ranked_file = glob.glob(gen_folder_name + "*_ranked.smi")

        for rank_file in ranked_file:
            # Check number of lines
            num_lines = sum(1 for line in open(rank_file,'r')) 
        
        # with open(rank_file, 'r') as rf:
            
        #     line = rf.readlines()

        #     data_all.append(line)
            
        #     parts = line.split('\t')      # split line into parts seperated by 4-spaces

        #     for part in parts:

        #         if part is "\n":
        #             continue
                

        #         print(part)
        with open(rank_file, 'r') as rf:
            
            for line in rf.readlines():
                line = line.replace("\n","")
                parts = line.split('\t')      # split line into parts seperated by 4-spaces
                
                name = parts[1]
                if name in data_dict.keys():
                    previous_entry = data_dict[name]
                    previous_score = previous_entry[-2]
                    current_score = parts[-2]
                    if float(current_score) < float(previous_score):
                        data_dict[name] = parts
                else:
                    data_dict[name] = parts


    list_of_lig = [data_dict[x] for x in data_dict.keys()]
    list_of_lig.sort(key = lambda x: float(x[-2]),reverse = False)



    with open(concat_file, "w") as CF:
        for line in list_of_lig:


            output_line = "\t".join(line) + "\n"
            CF.write(output_line)

if __name__ == "__main__":
    concat_ranked_smis(infolder, folder_list)


