""" Plots AutoGrow Run"""
import __future__

import os
import glob
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt  # type: ignore


def get_usable_format(infile: str) -> List[List[str]]:
    """
    This code takes a string for an file which is formatted as an .smi file. It
    opens the file and reads in the components into a usable list.

    The .smi must follow the following format for each line:
        MANDATORY INFO
            part 1 is the SMILES string
            part 2 is the SMILES name/ID

        Optional info
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating
                - For default setting this is the Docking score
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)

    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """
    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_smiles = []

    if os.path.exists(infile) is False:
        print(f"\nFile of Source compounds does not exist: {infile}\n")
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts separated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts separated by 4-spaces

            choice_list = [parts[i] for i in range(len(parts))]
            usable_smiles.append(choice_list)

    return usable_smiles


def get_average_score_per_gen(
    infolder: str, folder_list: List[str]
) -> Dict[str, Union[float, str]]:
    """
    This script will get the average docking score from the ranked .smi file
    from each generation.

    Inputs:
    :param str infolder: the path of the folder which has all of the
        generation folders
    :param list folder_list: a list of generation folders for each generation
        within infolder

    Returns:
    :returns: list usable_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """
    average_affinity_dict = {}
    for gen_folder in folder_list:
        gen_folder_name = infolder + gen_folder + os.sep
        ranked_file = glob.glob(f"{gen_folder_name}*_ranked.smi")

        for rank_file in ranked_file:
            # write as a tab delineated .smi file
            with open(rank_file, "r") as f:
                gen_affinity_sum = 0.0
                num_lines_counter = 0.0
                for line in f:
                    line = line.replace("\n", "")
                    parts = line.split(
                        "\t"
                    )  # split line into parts separated by 4-spaces

                    choice_list = [parts[i] for i in range(len(parts))]
                    gen_affinity_sum = gen_affinity_sum + float(choice_list[-2])
                    num_lines_counter = num_lines_counter + 1.0

            gen_affinity_average = gen_affinity_sum / num_lines_counter

            gen_num = os.path.basename(rank_file).split("_")[1]
            gen_name = f"generation_{gen_num}"
            average_affinity_dict[gen_name] = gen_affinity_average

    print_gens(average_affinity_dict)
    return average_affinity_dict


def get_average_top_score_per_gen(
    infolder: str, folder_list: List[str], top_score_per_gen: int
) -> Dict[str, Union[float, str]]:
    """
    This script will get the average docking score of the top N number of
    ligands ranked .smi file from each generation.

    Inputs:
    :param str infolder: the path of the folder which has all of the
        generation folders
    :param list folder_list: a list of generation folders for each generation
        within infolder
    :param int top_score_per_gen: the number of ligands to determine the
        average score. ie) if top_score_per_gen=50 it will return the average of
        the top 50 scores.

    Returns:
    :returns: dict average_affinity_dict: dictionary of average affinity
        scores for top_score_per_gen number of ligands
    """
    average_affinity_dict = {}

    for gen_folder in folder_list:
        gen_folder_name = infolder + gen_folder + os.sep
        ranked_file = glob.glob(f"{gen_folder_name}*_ranked.smi")

        for rank_file in ranked_file:
            # Check number of lines
            num_lines = 0
            with open(rank_file, "r") as rf:
                for _ in rf:
                    num_lines = num_lines + 1

            if num_lines >= top_score_per_gen:
                # read as a tab delineated .smi file
                with open(rank_file, "r") as f:
                    gen_affinity_sum = 0.0

                    for i, line in enumerate(f.readlines()):
                        if i >= top_score_per_gen:
                            break
                        line = line.replace("\n", "")
                        parts = line.split(
                            "\t"
                        )  # split line into parts separated by 4-spaces

                        choice_list = [parts[j] for j in range(len(parts))]
                        gen_affinity_sum = gen_affinity_sum + float(choice_list[-2])

                    gen_affinity_average = gen_affinity_sum / top_score_per_gen

                    gen_num = os.path.basename(rank_file).split("_")[1]
                    gen_name = f"generation_{gen_num}"
                    average_affinity_dict[gen_name] = gen_affinity_average

            else:
                gen_num = os.path.basename(rank_file).split("_")[1]
                gen_name = f"generation_{gen_num}"
                average_affinity_dict[gen_name] = "N/A"

    print_gens(average_affinity_dict)
    return average_affinity_dict


def print_gens(average_affinity_dict: Dict[str, Union[float, str]]) -> None:
    """
    This prints out the average scores for each generation

    Inputs:
    :param dict average_affinity_dict: dictionary of average affinity scores
        for top_score_per_gen number of ligands
    """
    print("generation_number              average affinity score")
    affinity_keys = list(average_affinity_dict.keys())
    affinity_keys.sort(key=lambda x: int(x.split("_")[1]))
    for gen in affinity_keys:
        print(gen, "                  ", average_affinity_dict[gen])


def make_graph(
    dictionary: Dict[str, Union[float, str]]
) -> Tuple[Optional[List[int]], Optional[List[Union[float, str]]]]:
    """
    Because some generations may not have 50 ligands this basically checks to
    see if theres enough ligands and prepares lists to be plotted

    Inputs:
    :param dict dictionary: dictionary of average affinity scores for
        top_score_per_gen number of ligands

    Returns:
    :returns: list list_generations: list of ints for each generation to be
        plotted. if a generation lacks ligands to generate the average it will
        return "N/A"
    :returns: list list_of_scores: list of averages for each generation; if a
        generation lacks ligands to generate the average it will return "N/A"
    """
    list_generations = []
    list_of_gen_names = []
    list_of_scores = []

    for key, score in dictionary.items():
        list_of_gen_names.append(key)

        list_of_scores.append(score)

        gen = key.replace("generation_", "")

        gen = int(gen)
        list_generations.append(gen)
        list_of_gen_names.append(key)

    # Check that there are no N/A in list_of_scores
    # If any gen has an N/A we will not plot the entire line
    for i in list_of_scores:
        if i == "N/A":
            return None, None

    return list_generations, list_of_scores


def run_plotter(
    params: Dict[str, Any],
    dict_of_averages: Dict[str, Dict[str, Union[float, str]]],
    outfile: str,
) -> None:
    """
    This plots the averages into a matplotlib figure. It will require you to
    answer questions about titles and labels

    Inputs:
    :param dict params: dict of user variables which will govern how the
        programs runs
    :param dict dict_of_averages: a dictionary of dictionaries containing the
        average of each generation for the top 50,20, 10, and 1 ligand(s) and the
        overall average for each generation.
    :param str outfile: Path for the output file for the plot
    """
    average_affinity_dict = dict_of_averages["average_affinity_dict"]
    top_fifty_dict = dict_of_averages["top_fifty_dict"]
    top_twenty_dict = dict_of_averages["top_twenty_dict"]
    top_ten_dict = dict_of_averages["top_ten_dict"]
    top_one_dict = dict_of_averages["top_one_dict"]

    # print("Graphing Overall Average")
    list_generations_average, list_of_scores_average = make_graph(average_affinity_dict)
    print_fifty = all(top_fifty_dict[key] != "N/A" for key in top_fifty_dict.keys())
    list_generations_fifty = []
    list_of_scores_fifty = []
    if print_fifty:
        list_generations_fifty, list_of_scores_fifty = make_graph(top_fifty_dict)

    print_twenty = all(top_twenty_dict[key] != "N/A" for key in top_twenty_dict.keys())
    list_generations_twenty = []
    list_of_scores_twenty = []
    if print_twenty:
        list_generations_twenty, list_of_scores_twenty = make_graph(top_twenty_dict)

    print_ten = all(top_ten_dict[key] != "N/A" for key in top_ten_dict.keys())
    list_generations_ten, list_of_scores_ten = make_graph(top_ten_dict)
    # print("Graphing top_one_dict")
    list_generations_one, list_of_scores_one = make_graph(top_one_dict)
    # print("")

    ax = plt.subplot(111)

    ax.plot(
        list_generations_average, list_of_scores_average, color="b", label="Average"
    )
    if print_fifty:
        ax.plot(list_generations_fifty, list_of_scores_fifty, color="c", label="Top 50")

    if print_twenty:
        ax.plot(
            list_generations_twenty, list_of_scores_twenty, color="m", label="Top 20"
        )

    if print_ten:
        ax.plot(list_generations_ten, list_of_scores_ten, color="g", label="Top 10")

    ax.plot(list_generations_one, list_of_scores_one, color="r", label="Top 1")

    ax.set_ylim()

    receptor_name = os.path.basename(params["receptor_path"])
    scoring_type = params["scoring_choice"]
    docking_type = params["scoring_choice"]
    num_lig = (
        int(params["number_of_mutants"])
        + int(params["number_of_crossovers"])
        + int(params["number_elitism_advance_from_previous_gen"])
    )
    number_of_conf_per_lig = str(params["max_variants_per_compound"])

    # Get Customizations
    title_of_figure = f"{scoring_type} Scores for {receptor_name} using {docking_type}"
    plt.title(title_of_figure, fontweight="semibold")

    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.274), fontsize="small")
    number_of_lig_per_gen = str(num_lig)

    output = (
        f"{number_of_lig_per_gen} lig/gen"
        + "\n"
        + number_of_conf_per_lig
        + " variants/lig"
    )

    plt.text(
        5.4, -8.5, output, bbox=dict(facecolor="white", alpha=0.5), fontsize="small"
    )

    # legend1 = plt.legend([lines[i].get_label() for i in range(lines_leg)],
    #           loc='center left', bbox_to_anchor=(1, 0.274),fontsize='small')
    # legend2 = plt.legend([output],loc='center left',
    #           bbox_to_anchor=(1, 0.774),fontsize='small')
    # # help(plt.legend)
    # ax.add_artist(legend1)
    # ax.add_artist(legend2)

    ax.set_ylim()

    if "VINA" in str(scoring_type):
        y_label = "Docking Affinity (kcal/mol"
    else:
        y_label = "Fitness Score"

    plt.ylabel(y_label, fontweight="semibold")

    plt.xlabel("Generation Number", fontweight="semibold")

    try:
        if print_ten:
            plt.savefig(outfile, bbox_inches="tight", dpi=1000)
        else:
            # Remove the bbox_inches="tight" is necessary if
            # the plot is too small for
            plt.savefig(outfile, dpi=500)
    except Exception:
        printout = (
            "\nUNABLE TO CREATE PLOT: \n"
            + "Population size or number of generations was "
        )
        printout += "too small to effectively plot. \n"
        print(printout)


def print_data_table(
    infolder: str, folder_list: List[str]
) -> Dict[str, Dict[str, Union[float, str]]]:
    """
    This function takes a folder of an Autogrow Run and a list of all folders
    within the infolder, and finds the average of each generation, the average
    of the top 50,20, 10, and 1 ligand(s) in each generation.

    It prints the average docking score values in a table and returns that
    information as a dictionary of dictionaries.

    Inputs:
    :param str infolder: a string for the file path to a directory containing
        an Autogrow run. ie) "PATH/Run_0/"
    :param list folder_list: a list of every generation folders within the
        infolder

    Returns
    :returns: dict dict_of_averages: a dictionary of dictionaries containing
        the average of each generation for the top 50,20, 10, and 1 ligand(s) and
        the overall average for each generation.
    """
    print("Overall Scoring Average for all Compounds")
    average_affinity_dict = get_average_score_per_gen(infolder, folder_list)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", 50)
    top_fifty_dict = get_average_top_score_per_gen(infolder, folder_list, 50)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", 20)
    top_twenty_dict = get_average_top_score_per_gen(infolder, folder_list, 20)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", 10)
    top_ten_dict = get_average_top_score_per_gen(infolder, folder_list, 10)
    print("")
    print("Best Score per generation")
    print("Number of top scoring compounds: ", 1)
    top_one_dict = get_average_top_score_per_gen(infolder, folder_list, 1)
    print("")
    print("")
    return {
        "average_affinity_dict": average_affinity_dict,
        "top_fifty_dict": top_fifty_dict,
        "top_twenty_dict": top_twenty_dict,
        "top_ten_dict": top_ten_dict,
        "top_one_dict": top_one_dict,
    }


# Run Everything
def generate_figures(params: Dict[str, Any]) -> None:
    """
    This runs everything to make a line plot of the results of an Autogrow
    simulation.

    Inputs:
    :param dict params: dict of user variables which will govern how the
        programs runs
    """
    for _ in range(10):
        print("")
    infolder = params["output_directory"]

    outfile = f"{infolder}data_line_plot.png"

    all_folders_list = [
        f for f in sorted(os.listdir(infolder)) if os.path.isdir(infolder + f)
    ]
    folder_list = [
        folder
        for folder in all_folders_list
        if folder != "Data" and len(folder.split("_")) == 2
    ]
    folder_list.sort(key=lambda x: int(x.split("_")[1]))

    dict_of_averages = print_data_table(infolder, folder_list)
    run_plotter(params, dict_of_averages, outfile)
