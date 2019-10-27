import __future__

import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt

# SPECIFY WHICH RUN NUMBER
# infolder = "/home/jacob/Desktop/Outputfolder/Run_47/"
infolder = "/home/jspiegel/projects/AUTOGROW_PARP/Run_1/"
infolder = "/home/jacob/Desktop/PARP_AUTOGROW_DATA/Run_0/Run_0/"
infolder = sys.argv[1]

all_folders_list = [f for f in sorted(os.listdir(infolder)) if os.path.isdir(infolder+f)]

folder_list = []
for folder in all_folders_list:
    if folder!="Data" and len(folder.split('_'))==2:
        folder_list.append(folder)

folder_list.sort(key=lambda x: int(x.split('_')[1]))

top_score_per_gen_fifty = 50
top_score_per_gen_twenty = 20
top_score_per_gen_ten = 10
top_score_per_gen_one = 1


def get_average_score_per_gen(infolder, folder_list):
    average_score_list = []
    for gen_folder in folder_list:
        gen_folder_name = infolder + gen_folder + os.sep
        ranked_file = glob.glob(gen_folder_name + "*_ranked.smi")
        print("")
        print("ranked_files: ", ranked_file)
        print("")
        for rank_file in ranked_file:
            # write as a tab delineated .smi file
            with open(rank_file, 'r') as f:
                gen_affinity_sum = float(0.0)
                num_lines_counter = float(0.0)
                for line in f:
                    line = line.replace("\n","")
                    parts = line.split('\t')      # split line into parts seperated by 4-spaces

                    choice_list = []
                    for i in range(0,len(parts)):
                        choice_list.append(parts[i])

                    
                    gen_affinity_sum = gen_affinity_sum + float(choice_list[-2])
                    num_lines_counter = num_lines_counter + float(1.0)


            gen_affinity_average = gen_affinity_sum/num_lines_counter

            average_score_list.append(gen_affinity_average)

    average_affinity_dict = {}

    for i in range(0,len(average_score_list)):
        gen_name = "generation_{}".format(i)
        average_affinity_dict[gen_name] = average_score_list[i]

    print_gens(average_affinity_dict)
    return average_affinity_dict

def get_average_top_score_per_gen(infolder, folder_list, top_score_per_gen):
    average_score_list = []

    for gen_folder in folder_list:
        gen_folder_name = infolder + gen_folder + "/"
        ranked_file = glob.glob(gen_folder_name + "*_ranked.smi")

        for rank_file in ranked_file:
            # Check number of lines
            num_lines = sum(1 for line in open(rank_file,'r')) 
       
            if num_lines >= top_score_per_gen:
                # read as a tab delineated .smi file
                with open(rank_file, 'r') as f:
                    gen_affinity_sum = float(0.0)
            
                    for i,line in enumerate(f.readlines()):
                        if i >= top_score_per_gen:
                            break
                        line = line.replace("\n","")
                        parts = line.split('\t')      # split line into parts seperated by 4-spaces

                        choice_list = []
                        for i in range(0,len(parts)):
                            choice_list.append(parts[i])

                        gen_affinity_sum = gen_affinity_sum + float(choice_list[-2])
                        
                    gen_affinity_average = gen_affinity_sum/top_score_per_gen
                    average_score_list.append(gen_affinity_average)
                
            else:
                average_score_list.append("N/A")

    print("out of for loop")

    average_affinity_dict = {}
    for i in range(0,len(average_score_list)):
        gen_name = "generation_{}".format(i)
        average_affinity_dict[gen_name] = average_score_list[i]

    print_gens(average_affinity_dict)
    return average_affinity_dict

def print_gens(average_affinity_dict):
    print("generation_number              average affinity score")
    affinity_keys = list(average_affinity_dict.keys())
    affinity_keys.sort(key=lambda x: int(x.split('_')[1]))
    for gen in affinity_keys:
        print(gen,"                  ",average_affinity_dict[gen])



def make_graph(dictionary):

    list_generations = []
    list_of_gen_names = []
    list_of_scores = []
    print(dictionary)

    for key in dictionary.keys():
        print(key)
        list_of_gen_names.append(key)    

        score = dictionary[key]
        list_of_scores.append(score)

        gen = key.replace("generation_","")
        print(gen)
        gen = int(gen)
        list_generations.append(gen)
        list_of_gen_names.append(key)    

    enough=True
    for i in list_of_scores:
        if i is "N/A":
            enough=False
            break

    print("")
    print("")
    print("list_generations")
    print(list_generations)
    print("list_of_scores")
    print(list_of_scores)
    print("")
    print("PLOT STUFF:")

    plt.plot(list_generations, list_of_scores, line=2.0)
    plt.show()
    
if __name__ == "__main__":
    print("Overall Scoring Average for all Compounds")
    average_affinity_dict = get_average_score_per_gen(infolder, folder_list)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", top_score_per_gen_fifty)
    top_fifty_dict = get_average_top_score_per_gen(infolder, folder_list, top_score_per_gen_fifty)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", top_score_per_gen_twenty)
    top_twenty_dict = get_average_top_score_per_gen(infolder, folder_list, top_score_per_gen_twenty)
    print("") 
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", top_score_per_gen_ten)
    top_ten_dict = get_average_top_score_per_gen(infolder, folder_list, top_score_per_gen_ten)
    print("") 
    print("Best Score per generation")
    print("Number of top scoring compounds: ", 1)
    top_one_dict = get_average_top_score_per_gen(infolder, folder_list, 1)
    print("")
    print("")
    print("Graphing")
    make_graph(top_one_dict)
    print("")
