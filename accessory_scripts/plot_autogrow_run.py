"""
Plots a line plot of the average score for each generation of AutoGrow run.

Example submit:
    python autogrow4/accessory_scripts/plot_autogrow_run.py\
        -i $PATH/Run_1/Run_0/ \
        --plot_reference_lines [['Olaparib Score',-12.8,'y'],\
            ['Niraparib',-10.7,'k'],['NAD/NADH',-10.3,'purple'],\
                ['ADP-ribose',-9.3,'maroon']]
"""
import __future__

import os
import glob
import json
import copy
import argparse
from typing import Any, Dict, List, Optional, Tuple, Union
import matplotlib  # type: ignore
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt  # type: ignore
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.Chem import Lipinski
from rdkit import DataStructs
import prolif
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE


def read_smi_file(source_file: str):
    rdkit_molecules = []
    with open(source_file) as smiles_file:
        for tsv_line in smiles_file:
            prts = tsv_line.replace("    ", "\t").strip().split("\t")
            rdkit_molecules.append(Chem.MolFromSmiles(prts[0]))

    if None in rdkit_molecules:
        raise Exception("An input molecule was not successfully read from " + source_file)

    return rdkit_molecules


def read_sdf_file(source_file: str):
    rdkit_molecules = []
    r = Chem.SDMolSupplier(source_file, sanitize=False)
    for compound in r:
        rdkit_molecules.append(compound)
    r.reset()
    return rdkit_molecules


def calc_diversity_scores(reference_comp, new_comps):
    """
    Calculate diversity scores between a reference compound and a list of new
    compounds.

    This function computes Morgan Fingerprints for each molecule and calculates
    pairwise Dice Similarity scores (1.0 means a perfect match, 0.0 means no
    match at all).

    Args:
        reference_comp: RDKit molecule.
        new_comps: list of RDKit molecules.

    Returns:
        List[int]: similarity values.

    Note:
        - Lower diversity scores indicate more diverse molecules.
        - Removes any None entries from the input list.
        - Uses Morgan Fingerprints with radius 10 and feature-based encoding.
    """
    similarity_values = []

    reference_comp_fp = GetMorganFingerprint(reference_comp, radius=10, useFeatures=True)
    for mol in new_comps:
        if mol is not None:
            mol_fp = GetMorganFingerprint(mol, radius=10, useFeatures=True)

            # if DiceSimilarity=1.0 it is a perfect match, the smaller the
            # number, the more diverse it is.
            diversity_score = DataStructs.DiceSimilarity(reference_comp_fp, mol_fp)
            similarity_values.append(diversity_score)

    return similarity_values


def calc_interaction_fp_per_generation(params: Dict[str, Any], infolder: str, interactions: List[str], analyze_gen_0: bool):
    """
    Calculate percent of interaction types per generation.

    Args:
    :param dict params: dict of user variables which will govern how the
        program runs.
    :param str infolder: Path to the folder containing all generation folders.
    :param str interactions: List of interaction types to calculate.
    :param bool analyze_gen_0: If True, the generation 0 related to the input compounds is analyzed.

    Returns:
        np.array: percent of interaction types per generation.
        list: list containing labels of the rows
        list: list containing labels of the columns
    """
    receptor = Chem.MolFromPDBFile(params["receptor_path"], sanitize=True)
    fingerprints = prolif.Fingerprint(interactions)

    all_interactions = []
    number_interactions_per_gen = {}
    num_compounds_per_generation = {}
    for gen_folder in os.listdir(infolder):
        if "generation" in gen_folder and os.path.isdir(os.path.join(infolder, gen_folder)):
            gen_folder_name = infolder + os.sep + gen_folder + os.sep
            rank_file = glob.glob(f"{gen_folder_name}*_ranked.smi")[0]
            gen_num = os.path.basename(rank_file).split("_")[1]
            num_compounds_per_generation[gen_num] = 0

            with open(rank_file, "r") as f:
                for line in f:
                    line = line.replace("\n", "")
                    parts = line.split(
                        "\t"
                    )  # split line into parts separated by 5-spaces

                    docked_comp_path = [parts[i] for i in range(len(parts))][4]
                    if docked_comp_path != "None":
                        num_compounds_per_generation[gen_num] = num_compounds_per_generation[gen_num] + 1
                        docked_comp = read_sdf_file(docked_comp_path)[0]
                        prot = prolif.Molecule.from_rdkit(receptor)
                        lig = prolif.Molecule.from_rdkit(docked_comp)
                        ifp = fingerprints.generate(lig, prot, metadata=True)
                        df = prolif.to_dataframe({0: ifp}, fingerprints.interactions)

                        for column_name in df.columns:
                            interaction = column_name[1].split(".")[0] + "_" + column_name[2]
                            if interaction not in number_interactions_per_gen:
                                number_interactions_per_gen[interaction] = {}
                                all_interactions.append(interaction)
                            if gen_num not in number_interactions_per_gen[interaction]:
                                number_interactions_per_gen[interaction][gen_num] = 0
                            number_interactions_per_gen[interaction][gen_num] = \
                                number_interactions_per_gen[interaction][gen_num] + 1

    all_generations = []
    result_matrix = np.zeros((len(num_compounds_per_generation), len(all_interactions)))
    for gen_num_idx in range(len(num_compounds_per_generation) - (1 if not analyze_gen_0 else 0)):
        gen_num = (gen_num_idx + 1) if not analyze_gen_0 else gen_num_idx
        all_generations.append(f"generation_{gen_num}")
        for interaction_num in range(len(all_interactions)):
            interaction = all_interactions[interaction_num]
            if str(gen_num) in number_interactions_per_gen[interaction]:
                percent_of_interactions = number_interactions_per_gen[interaction][str(gen_num)]
                percent_of_interactions = round(float(percent_of_interactions /
                                                      num_compounds_per_generation[str(gen_num)]), 2)
                result_matrix[gen_num_idx, interaction_num] = percent_of_interactions

    return result_matrix, all_generations, all_interactions


def generate_tSNE_scatterplot(infolder: str, params: Dict[str, Any], outfile: str):
    """
    Generate a t-SNE scatterplot for the compounds generated by Autogrow and the input compounds.

    This function computes Morgan Fingerprints for each molecule to be used as the features for the
    t-SNE analsyis.

    Args:
    :param str infolder: Path to the folder containing all generation folders.
    :param dict params: dict of user variables which will govern how the
        program runs.
    :param str outfile: Path for the output file for the plot
    """
    label_list = []
    result_matrix = []
    for gen_folder in os.listdir(infolder):
        if "generation" in gen_folder and os.path.isdir(os.path.join(infolder, gen_folder)):
            gen_folder_name = infolder + os.sep + gen_folder + os.sep
            rank_file = glob.glob(f"{gen_folder_name}*_ranked.smi")[0]
            gen_num = os.path.basename(rank_file).split("_")[1]

            new_comps = []
            with open(rank_file, "r") as f:
                for line in f:
                    line = line.replace("\n", "")
                    parts = line.split(
                        "\t"
                    )  # split line into parts separated by 4-spaces

                    choice_list = [parts[i] for i in range(len(parts))]
                    new_comps.append(Chem.MolFromSmiles(choice_list[0]))

                for mol in new_comps:
                    if mol is not None:
                        mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=10, nBits=2048)
                        result_matrix.append(list(mol_fp))
                        label_list.append("generation_" + str(gen_num))

    result_matrix = np.array(result_matrix)
    tsne_model = TSNE(n_components=2, random_state=0)
    tsne_results = tsne_model.fit_transform(result_matrix)
    tsne_results = np.vstack((tsne_results.T, np.array(label_list))).T

    tsne_df = pd.DataFrame(data=tsne_results, columns=("Dim_1", "Dim_2", "Generations"))
    tsne_df["Dim_1"] = pd.to_numeric(tsne_df["Dim_1"])
    tsne_df["Dim_2"] = pd.to_numeric(tsne_df["Dim_2"])
    grouped_tsne_df = tsne_df.groupby('Generations')

    x = grouped_tsne_df['Dim_1'].apply(lambda x: x.values)
    y = grouped_tsne_df['Dim_2'].apply(lambda y: y.values)

    plt.clf()
    fig, ax = plt.subplots()
    for gen_id in range(len(x)):
        gen_name = "generation_" + str(gen_id)
        x_gen = x[gen_id]
        y_gen = y[gen_id]
        ax.scatter(x_gen, y_gen, label=gen_name, c='black' if gen_id == 0 else None)

    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

    # Add titles and labels
    plt.xlabel("Dimension 1", fontweight="semibold")
    plt.ylabel("Dimension 2", fontweight="semibold")

    plt.savefig(outfile, bbox_inches="tight", format=params["outfile_format"], dpi=1000)


def calc_diversity_scores_per_generation(infolder: str):
    """
    Calculate diversity scores between the compounds corresponding to each generation.

    This function computes Morgan Fingerprints for each molecule and calculates
    pairwise Dice Similarity scores (1.0 means a perfect match, 0.0 means no
    match at all).

    Args:
    :param str infolder: Path to the folder containing all generation folders.

    Returns:
        Dict[str, List[float]]: Dictionary of similarity values, with generation ID as keys.

    Note:
        - Lower diversity scores indicate more diverse molecules.
        - Removes any None entries from the input list.
        - Uses Morgan Fingerprints with radius 10 and feature-based encoding.
    """
    average_dict = {}
    for gen_folder in os.listdir(infolder):
        if "generation" in gen_folder and os.path.isdir(os.path.join(infolder, gen_folder)):
            gen_folder_name = infolder + os.sep + gen_folder + os.sep
            rank_file = glob.glob(f"{gen_folder_name}*_ranked.smi")[0]
            gen_num = os.path.basename(rank_file).split("_")[1]

            new_comps = []
            with open(rank_file, "r") as f:
                for line in f:
                    line = line.replace("\n", "")
                    parts = line.split(
                        "\t"
                    )  # split line into parts separated by 4-spaces

                    choice_list = [parts[i] for i in range(len(parts))]
                    new_comps.append(Chem.MolFromSmiles(choice_list[0]))

            similarity_values = []
            for i in range(len(new_comps)):
                reference_comp = new_comps[i]
                if reference_comp is not None:
                    reference_comp_fp = GetMorganFingerprint(reference_comp, radius=10, useFeatures=True)
                    for j in range(i + 1, len(new_comps)):
                        mol = new_comps[j]
                        if mol is not None:
                            mol_fp = GetMorganFingerprint(mol, radius=10, useFeatures=True)

                            # if DiceSimilarity=1.0 it is a perfect match, the smaller the number,
                            # the more diverse it is.
                            diversity_score = DataStructs.DiceSimilarity(reference_comp_fp, mol_fp)
                            if diversity_score >= 0.9:
                                print("Compounds " + str(i + 1) + " and " + str(j + 1) + " of the generation "
                                      + str(gen_num) + " are similar in: " + str(diversity_score))
                            similarity_values.append(diversity_score)

            gen_name = f"generation_{gen_num}"
            average_dict[gen_name] = similarity_values

    return average_dict


def get_similarity_list_per_input_comp(infolder: str, source_file: str) -> Dict[str, List[float]]:
    """
    Calculate the similarity between the input compounds and the compounds generated by Autogrow.

    Args:
        infolder (str): Path to the folder containing all generation folders.
        source_file (str): Path to the input files containing the compounds.

    Returns:
        Dict[str, List[float]]: Dictionary of similarity values, with compound ID as keys.
    """
    similarity_dict = {}

    if ".smi" in source_file:
        input_rdkit_molecules = read_smi_file(source_file)
    elif ".sdf" in source_file:
        input_rdkit_molecules = read_sdf_file(source_file)

    ranked_file = infolder + os.sep + "summary_ranked.sdf"
    ranked_molecules = read_sdf_file(ranked_file)

    for i in range(len(input_rdkit_molecules)):
        reference_mol = input_rdkit_molecules[i]
        similarity_dict["compound_" + str(i + 1)] = calc_diversity_scores(reference_mol, ranked_molecules)

    return similarity_dict


def get_ave_similarity_per_generated_comp(infolder: str, source_file: str) -> Dict[str, float]:
    """
    Calculate the average similarity between the compounds generated by Autogrow and the input compounds.

    Args:
        infolder (str): Path to the folder containing all generation folders.
        source_file (str): Path to the input files containing the compounds.

    Returns:
        Dict[str, float]: Dictionary of average similarity values, with compound ID as keys.

    Notes:
        if a generated molecule is not successfully read, a similarity equal to 1 is assigned.
    """
    average_similarity_dict = {}

    if ".smi" in source_file:
        input_rdkit_molecules = read_smi_file(source_file)
    elif ".sdf" in source_file:
        input_rdkit_molecules = read_sdf_file(source_file)

    ranked_file = infolder + os.sep + "summary_ranked.sdf"
    ranked_molecules = read_sdf_file(ranked_file)

    for i in range(len(ranked_molecules)):
        reference_mol = ranked_molecules[i]
        if reference_mol is None:
            average_similarity_dict["compound_" + str(i + 1)] = -0.05
        else:
            similarity_values = calc_diversity_scores(reference_mol, input_rdkit_molecules)
            similarity = 0.0;
            for value in similarity_values:
                similarity = similarity + value
            average_similarity_dict["compound_" + str(i + 1)] = (similarity / len(input_rdkit_molecules))

    return average_similarity_dict


def get_efficiency_per_generated_comp(infolder: str) -> Dict[str, float]:
    """
    Calculate the ligand efficiency for the compounds generated by Autogrow.

    Args:
        infolder (str): Path to the folder containing all generation folders.

    Returns:
        Dict[str, float]: Dictionary of average similarity values, with compound ID as keys.

    Notes:
        if a generated molecule is not successfully read, an efficiency equal to 0 is assigned.
    """
    average_similarity_dict = {}

    ranked_file = infolder + os.sep + "summary_ranked.sdf"
    ranked_molecules = read_sdf_file(ranked_file)

    for i in range(len(ranked_molecules)):
        mol = ranked_molecules[i]
        if mol is None:
            average_similarity_dict["compound_" + str(i + 1)] = 0.0
        else:
            efficiency = float(mol.GetProp('Docking_Score')) / Lipinski.HeavyAtomCount(mol)
            average_similarity_dict["compound_" + str(i + 1)] = efficiency

    return average_similarity_dict


def get_score_list_per_gen(infolder: str, ligand_efficiency: bool) -> Dict[str, List[float]]:
    """
    Get the docking scores (or ligand efficiencies) for each generation.

    This function reads ranked .smi files from each generation folder and
    gets the docking scores (or ligand efficiencies) for each generation.

    Inputs:
    :param str infolder: Path to the folder containing all generation folders.
    :param bool ligand_efficiency: If True, the ligand efficiency will be calculated instead
        of the affinity score.

    Returns:
        Dict[str, List[float]]: Dictionary of lists of docking scores (or ligand efficiencies)
        for each generation, with generation names as keys.
    """
    average_dict = {}
    for gen_folder in os.listdir(infolder):
        if os.path.isdir(os.path.join(infolder, gen_folder)):
            gen_folder_name = infolder + os.sep + gen_folder + os.sep
            ranked_file = glob.glob(f"{gen_folder_name}*_ranked.smi")

            for rank_file in ranked_file:
                gen_list = []
                with open(rank_file, "r") as f:
                    for line in f:
                        line = line.replace("\n", "")
                        parts = line.split(
                            "\t"
                        )  # split line into parts separated by 4-spaces

                        choice_list = [parts[i] for i in range(len(parts))]
                        gen_list.append(float(choice_list[2]) if not ligand_efficiency else float(
                            choice_list[2]) / Lipinski.HeavyAtomCount(Chem.MolFromSmiles(
                                choice_list[0], sanitize=False)))

                gen_num = os.path.basename(rank_file).split("_")[1]
                gen_name = f"generation_{gen_num}"
                average_dict[gen_name] = gen_list

    return average_dict


def get_average_score_per_gen(infolder: str, ligand_efficiency: bool) -> Dict[str, List[float]]:
    """
    Calculate the average docking score (or ligand efficiency) for each generation.

    This function reads ranked .smi files from each generation folder and
    calculates the average docking score (or ligand efficiency) for each generation.

    Inputs:
    :param str infolder: Path to the folder containing all generation folders.
    :param bool ligand_efficiency: If True, the ligand efficiency will be calculated instead
        of the average affinity score.

    Returns:
        Dict[str, List[float]]: Dictionary of average docking scores (or ligand efficiencies)
        for each generation, with generation names as keys.
    """
    average_affinity_dict = {}
    for gen_folder in os.listdir(infolder):
        if os.path.isdir(os.path.join(infolder, gen_folder)):
            gen_folder_name = infolder + os.sep + gen_folder + os.sep
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
                        gen_affinity_sum = gen_affinity_sum + float(choice_list[2]) \
                            if not ligand_efficiency \
                            else gen_affinity_sum + (float(choice_list[2]) /
                                                     Lipinski.HeavyAtomCount(
                                                         Chem.MolFromSmiles(choice_list[0], sanitize=False)))
                        num_lines_counter = num_lines_counter + 1.0

                gen_affinity_average = gen_affinity_sum / num_lines_counter

                gen_num = os.path.basename(rank_file).split("_")[1]
                gen_name = f"generation_{gen_num}"
                average_affinity_dict[gen_name] = gen_affinity_average

    print_gens(average_affinity_dict)
    return average_affinity_dict


def get_average_top_score_per_gen(infolder: str, top_score_per_gen: int, ligand_efficiency: bool) \
        -> Dict[str, List[float]]:
    """
    This script will get the average docking score (or ligand efficiency) for the top N number of
    ligands ranked .smi file from each generation.

    Inputs:
    :param str infolder: the path of the folder which has all the
        generation folders
    :param int top_score_per_gen: the number of ligands to determine the
        average score. i.e.,) if top_score_per_gen=50 it will return the average of
        the top 50 scores.
    :param bool ligand_efficiency: if True, the ligand efficiency will be calculated instead
        of the average affinity score.

    Returns:
    :returns: dict average_affinity_dict: dictionary of average affinity
        scores for top_score_per_gen number of ligands
    """

    average_affinity_dict = {}

    for gen_folder in os.listdir(infolder):
        if os.path.isdir(os.path.join(infolder, gen_folder)):
            gen_folder_name = infolder + os.sep + gen_folder + os.sep
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
                            gen_affinity_sum = gen_affinity_sum + float(choice_list[2]) \
                                if not ligand_efficiency \
                                else gen_affinity_sum + (float(choice_list[2]) /
                                                         Lipinski.HeavyAtomCount(
                                                             Chem.MolFromSmiles(choice_list[0], sanitize=False)))

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
    :param dict average_affinity_dict: dictionary of average values
        for top_score_per_gen number of ligands
    """

    print("generation_number              average")
    affinity_keys = list(average_affinity_dict.keys())
    affinity_keys.sort(key=lambda x: int(x.split("_")[1]))
    for gen in affinity_keys:
        print(gen, "                  ", average_affinity_dict[gen])


def make_graph(dictionary: Dict[str, Union[float, str]], analyze_gen_0: bool) -> Tuple[Optional[List[int]], Optional[List[float]]]:
    """
    Because some generations may not have 50 ligands this basically checks to see if
    there are enough ligands and prepares lists to be plotted

    Inputs:
    :param dict dictionary: dictionary of average affinity scores for
        top_score_per_gen number of ligands
    :param bool analyze_gen_0: if True, the compounds of the generation 0, which are the input compounds, will be
    analyzed.

    Returns:
    :returns: list list_generations: list of ints for each generation to be plotted.
        if a generation lacks ligands to generate the average it will return "N/A"
    :returns: list list_of_scores: list of averages for each generation;
        if a generation lacks ligands to generate the average it will return "N/A"
    """
    list_generations = []
    list_of_gen_names = []
    list_of_scores = []

    for gen in range(len(dictionary) - (1 if not analyze_gen_0 else 0)):
        key = "generation_" + (str(gen + 1) if not analyze_gen_0 else str(gen))
        score = dictionary[key]

        # print(key)
        list_of_gen_names.append(key)
        list_of_scores.append(score)

        list_generations.append((gen + 1) if not analyze_gen_0 else gen)
        list_of_gen_names.append(key)

    for i in list_of_scores:
        if i == "N/A":
            return None, None

    return list_generations, list_of_scores


def run_score_plotter(
    params: Dict[str, Any],
    dict_of_averages: Dict[str, Dict[str, Union[float, str]]],
    outfile: str,
    ligand_efficiency: bool,
    analyze_gen_0: bool
) -> None:
    """
    This plots the averages into a matplotlib figure. It will require you to
    answer questions about titles and labels

    Inputs:
    :param dict params: dict of user variables which will govern how the
        program runs.
    :param dict dict_of_averages: A dictionary of dictionaries containing the
        average of each generation for the top 50,20, 10, and 1 ligand(s) and the
        overall average for each generation.
    :param str outfile: Path for the output file for the plot
    :param bool ligand_efficiency: If True, the ligand efficiency will be calculated instead
        of the average affinity score.
    :param bool analyze_gen_0: If True, the compounds of the generation 0, which are the input compounds, will be
        analyzed.
    """

    top_15 = dict_of_averages["top_15"]
    top_10 = dict_of_averages["top_10"]
    top_5 = dict_of_averages["top_5"]
    top_1 = dict_of_averages["top_1"]

    list_generations_15 = []
    list_of_scores_15 = []
    list_generations_10 = []
    list_of_scores_10 = []

    average_affinity_dict = dict_of_averages["average_dict"]
    list_generations_average, list_of_scores_average = make_graph(average_affinity_dict, analyze_gen_0)

    print_top_15 = all(top_15[key] != "N/A" for key in top_15.keys())
    if print_top_15:
        list_generations_15, list_of_scores_15 = make_graph(top_15, analyze_gen_0)

    print_top_10 = all(top_10[key] != "N/A" for key in top_10.keys())
    if print_top_10:
        list_generations_10, list_of_scores_10 = make_graph(top_10, analyze_gen_0)

    list_generations_ten, list_of_scores_ten = make_graph(top_5, analyze_gen_0)
    list_generations_one, list_of_scores_one = make_graph(top_1, analyze_gen_0)

    plt.clf()
    ax = plt.subplot(111)

    ax.plot(
        list_generations_average, list_of_scores_average, color="b", label="Average"
    )

    if print_top_15:
        ax.plot(list_generations_15, list_of_scores_15, color="c", label="Top 15")

    if print_top_10:
        ax.plot(list_generations_10, list_of_scores_10, color="m", label="Top 10")

    ax.plot(list_generations_ten, list_of_scores_ten, color="g", label="Top 5")
    ax.plot(list_generations_one, list_of_scores_one, color="r", label="Top 1")

    if params["plot_reference_lines"] is not None:
        for ref_info in params["plot_reference_lines"]:
            ax.axhline(
                y=ref_info[1], color=ref_info[2], linestyle=":", label=ref_info[0]
            )

    ax.set_ylim()

    # Get Customizations
    receptor_name = os.path.basename(params["receptor_path"])
    title_of_figure = f"Docking scores for {receptor_name}" if not ligand_efficiency else \
        f"Ligand efficiencies for {receptor_name}"
    plt.title(title_of_figure, fontweight="semibold")

    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.274), fontsize="small")
    ax.set_ylim()

    scoring_type = params["vina_like_executable"]
    if "vina" in str(scoring_type):
        y_label = "Docking Affinity (kcal/mol)" if not ligand_efficiency else "Ligand Efficiency"
    else:
        y_label = "Fitness Score" if not ligand_efficiency else "Ligand Efficiency"
    plt.ylabel(y_label, fontweight="semibold")

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel("Generation Number", fontweight="semibold")

    plt.savefig(outfile, bbox_inches="tight", format=params["outfile_format"], dpi=1000)


def run_boxplot(
    params: Dict[str, Any],
    dictionary_of_values: Dict[str, List[float]],
    outfile: str,
    key_start_with: str,
    x_label: str,
    y_label: str,
    title_of_figure: str = None,
    analyze_gen_0: bool = False,
) -> None:
    data = []
    yticklabels = []
    for i in range(len(dictionary_of_values) - (1 if not analyze_gen_0 else 0)):
        key = key_start_with + "_" + (str(i + 1) if not analyze_gen_0 else str(i))
        data.append(dictionary_of_values[key])
        yticklabels.append(key)

    plt.clf()
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)

    # Creating plot with customizations
    bp = ax.boxplot(data, patch_artist=True, notch=False, vert=0,
                    showmeans=True,
                    meanprops={"markerfacecolor": "black", "markeredgecolor": "black"})
    for median in bp['medians']:
        median.set_color('black')

    # Setting x-axis labels

    plt.xlabel(x_label, fontweight="semibold")

    plt.ylabel(y_label, fontweight="semibold")
    ax.set_yticklabels(yticklabels)

    # Adding title
    title_of_figure = title_of_figure if title_of_figure is not None else ""
    plt.title(title_of_figure, fontweight="semibold")

    # Show plot
    plt.savefig(outfile, bbox_inches="tight", format=params["outfile_format"], dpi=1000)


def run_plotter(
    params: Dict[str, Any],
    dictionary_of_values: Dict[str, Any],
    outfile: str,
    key_start_with: str,
    x_label: str,
    y_label: str,
    title_of_figure: str = None,
) -> None:

    x = []
    y = []
    for i in range(len(dictionary_of_values)):
        y.append(dictionary_of_values[key_start_with + "_" + str(i + 1)])
        x.append(i + 1)

    # Create the plot
    plt.clf()
    plt.plot(x, y, marker='o', linestyle='-', color='b')

    # Add titles and labels
    plt.title(title_of_figure if title_of_figure is not None else "", fontweight="semibold")
    plt.xlabel(x_label if x_label is not None else "", fontweight="semibold")
    plt.ylabel(y_label if y_label is not None else "", fontweight="semibold")

    # Add a legend
    plt.legend()

    # Show plot
    plt.savefig(outfile, bbox_inches="tight", format=params["outfile_format"], dpi=1000)


def run_heatmap(
        params: Dict[str, Any],
        outfile: str,
        result_matrix,
        x_label_list,
        y_label_list,
        title_of_figure
) -> None:
    plt.clf()
    fig, ax = plt.subplots()
    ax.imshow(result_matrix, cmap="Blues")

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(range(len(x_label_list)), labels=x_label_list, rotation=45, ha="right", rotation_mode="anchor",
                  size="x-small")
    ax.set_yticks(range(len(y_label_list)), labels=y_label_list, size="x-small")

    # Loop over data dimensions and create text annotations.
    for i in range(len(y_label_list)):
        for j in range(len(x_label_list)):
            value = result_matrix[i, j]
            if value == 0.0:
                text = "0"
            elif value == 1.0:
                text = "1"
            else:
                text = "." + str(result_matrix[i, j]).split(".")[1]

            ax.text(j, i, text, ha="center", va="center", color="black", size="xx-small")

    ax.set_title(title_of_figure, size="small")
    plt.savefig(outfile, bbox_inches="tight", format=params["outfile_format"], dpi=1000)


def print_data_table(infolder: str, ligand_efficiency: bool) -> Dict[str, Any]:
    """
    This function takes a folder of an Autogrow Run and a list of all folders
    within the infolder, and finds the average of each generation, the average
    of the top 50,20, 10, and 1 ligand(s) in each generation.

    It prints the average docking score values in a table and returns that
    information as a dictionary of dictionaries.

    Inputs:
    :param str infolder: a string for the file path to a directory containing
        an Autogrow run.
    :param bool ligand_efficiency: If True, the ligand efficiency will be calculated instead
        of the average affinity score.

    Returns
    :returns: dict dict_of_averages: a dictionary of dictionaries containing
        the average of each generation for the top 50,20, 10, and 1 ligand(s) and
        the overall average for each generation.
    """

    print("Overall " + ("Scoring" if not ligand_efficiency else "Ligand Efficiency") + " Average for all Compounds")
    average_affinity_dict = get_average_score_per_gen(infolder, ligand_efficiency)
    print("")
    print("Average for Top " + ("Scoring Compounds " if not ligand_efficiency else "Ligand Efficiencies"))
    print("Number of top compounds: ", 15)
    top_15 = get_average_top_score_per_gen(infolder, 15, ligand_efficiency)
    print("")
    print("Average for Top " + ("Scoring Compounds " if not ligand_efficiency else "Ligand Efficiencies"))
    print("Number of top compounds: ", 10)
    top_10 = get_average_top_score_per_gen(infolder, 10, ligand_efficiency)
    print("")
    print("Average for Top " + ("Scoring Compounds " if not ligand_efficiency else "Ligand Efficiencies"))
    print("Number of top compounds: ", 5)
    top_5 = get_average_top_score_per_gen(infolder, 5, ligand_efficiency)
    print("")
    print("Best " + ("Scoring" if not ligand_efficiency else "Ligand Efficiency") + " per generation")
    print("Number of top compounds: ", 1)
    top_1 = get_average_top_score_per_gen(infolder, 1, ligand_efficiency)
    print("")
    print("")
    return {
        "average_dict": average_affinity_dict,
        "top_15": top_15,
        "top_10": top_10,
        "top_5": top_5,
        "top_1": top_1,
    }


def generate_figures(params: Dict[str, Any], analyze_gen_0: bool) -> None:
    """
    This runs everything to make a line plot for the results of an Autogrow simulation.

    Inputs:
    :param dict params: dict of user variables which will govern how the
        program runs.
    :param bool analyze_gen_0: If True, the compounds of the generation 0, which are the input compounds, will be
        analyzed.
    """

    infolder = params["infolder"]
    outfile = params["outfile"]

    dict_of_averages = print_data_table(infolder, False)
    run_score_plotter(params, dict_of_averages,
                      outfile + os.sep + "plotter_by_generation_for_docking_scores." + params["outfile_format"],
                      ligand_efficiency=False,
                      analyze_gen_0=analyze_gen_0)

    dict_of_averages = print_data_table(infolder, True)
    run_score_plotter(params, dict_of_averages,
                      outfile + os.sep + "plotter_by_generation_for_ligand_efficiencies." + params["outfile_format"],
                      ligand_efficiency=True,
                      analyze_gen_0=analyze_gen_0)

    dict_of_score_lists = get_score_list_per_gen(infolder, False)
    scoring_type = params["vina_like_executable"]
    if "vina" in str(scoring_type):
        x_label = "Docking Affinity (kcal/mol)"
    else:
        x_label = "Fitness Score"
    receptor_name = os.path.basename(params["receptor_path"])
    run_boxplot(params, dict_of_score_lists,
                outfile + os.sep + "boxplot_by_generation_for_docking_scores." + params["outfile_format"],
                key_start_with="generation",
                x_label=x_label,
                y_label="Number of Generations",
                title_of_figure="Docking Scores for " + receptor_name,
                analyze_gen_0=analyze_gen_0)

    dict_of_score_lists = get_score_list_per_gen(infolder, True)
    x_label = "Ligand Efficiency"
    receptor_name = os.path.basename(params["receptor_path"])
    run_boxplot(params, dict_of_score_lists,
                outfile + os.sep + "boxplot_by_generation_for_ligand_efficiencies." + params["outfile_format"],
                key_start_with="generation",
                x_label=x_label,
                y_label="Number of Generations",
                title_of_figure="Ligand_efficiencies for " + receptor_name,
                analyze_gen_0=analyze_gen_0)

    source_file = str(params["source_compound_file"])
    dict_of_similarity_lists = get_similarity_list_per_input_comp(infolder, source_file)
    run_boxplot(params, dict_of_similarity_lists,
                outfile + os.sep + "boxplot_of_similarity_between_input_and_new_compounds." + params["outfile_format"],
                key_start_with="compound",
                x_label="Dice Similarity Values",
                y_label="Input Compounds")

    dict_of_score_lists = calc_diversity_scores_per_generation(infolder)
    run_boxplot(params, dict_of_score_lists,
                outfile + os.sep + "boxplot_by_generation_for_similarities." + params["outfile_format"],
                key_start_with="generation",
                x_label="Dice Similarity Values",
                y_label="Number of Generations",
                analyze_gen_0=analyze_gen_0)

    dict_of_similarity_lists = get_ave_similarity_per_generated_comp(infolder, source_file)
    run_plotter(params, dict_of_similarity_lists,
                outfile + os.sep + "plotter_of_ave_similarities_between_each_new_compound_and_input_compounds." + params["outfile_format"],
                key_start_with="compound",
                x_label="New compounds (ID) sorted from the highest to lowest affinity",
                y_label="Average similarity")

    dict_of_efficiency_lists = get_efficiency_per_generated_comp(infolder)
    run_plotter(params, dict_of_efficiency_lists,
                outfile + os.sep + "plotter_of_ligand_efficiency_for_every_new_compound." + params[
                    "outfile_format"],
                key_start_with="compound",
                x_label="New compounds (ID) sorted from the highest to lowest affinity",
                y_label="Ligand efficiency")

    generate_tSNE_scatterplot(infolder=infolder,
                              outfile=outfile + os.sep + "tsne_for_input_and_new_compounds." + params["outfile_format"],
                              params=params)

    interactions = ["Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "Anionic", "Cationic", "CationPi", "PiCation"]
    result_matrix, y_label_list, x_label_list = calc_interaction_fp_per_generation(params, infolder, interactions,
                                                                                   analyze_gen_0)
    run_heatmap(params,
                outfile + os.sep + "heatmap_of_all_interactions_by_generation." + params["outfile_format"],
                result_matrix,
                x_label_list,
                y_label_list,
                "Analysis of Interactions")

    for interaction in interactions:
        pos_for_interaction = [i for i in range(len(x_label_list)) if interaction in x_label_list[i]]
        if len(pos_for_interaction) > 0:
            result_matrix_c = result_matrix[:, pos_for_interaction]
            x_label_list_c = [x_label_list[i] for i in pos_for_interaction]
            run_heatmap(params,
                        outfile + os.sep + f"heatmap_of_{interaction}_interactions_by_generation." + params["outfile_format"],
                        result_matrix_c,
                        x_label_list_c,
                        y_label_list,
                        f"Analysis of {interaction} Interactions")


def retrieve_vars_dict(autogrow_vars_json: str) -> Dict[str, Any]:
    """
    This will retrieve a variable dictionary from a AutoGrow vars json file.

    Inputs:
    :param str autogrow_vars_json: path to AutoGrow json variable file
    Returns:
    :returns: dict params: a dictionary of variable to use
    """
    if os.path.exists(autogrow_vars_json) is False:
        raise Exception(
            "variable file could not be found. It should be the \
            vars.json file written by AutoGrow in the output folder of the run."
        )
    try:
        with open(autogrow_vars_json, "r") as f:
            params = json.load(f)
    except Exception as e:
        raise Exception(
            "variable file would not import. It should be the \
            vars.json file written by AutoGrow in the output folder of the run."
        ) from e
    return params


def process_inputs(inputs: Dict[str, Any]) -> Dict[str, Any]:
    """
    This will handle processing all parameters.

    inputs:
    :params dict inputs: dictionary of argparse parameters
    Returns:
    :returns: dict vars_dict: dictionary of argparse parameters
    """

    # handle input information
    inputs["infolder"] = os.path.abspath(inputs["infolder"]) + os.sep
    if os.path.exists(inputs["infolder"]) is False:
        raise Exception(
            "Input folder {} does not\
            exist.".format(
                inputs["infolder"]
            )
        )

    # get vars dict from last run
    inputs["vars_json"] = inputs["infolder"] + "vars.json"
    if os.path.exists(inputs["vars_json"]) is False:
        raise Exception(
            "Input folder {} does not contain the vars.json file \
            necessary to run script. Please make sure the vars.json is in the \
            folder.".format(
                inputs["infolder"]
            )
        )

    try:
        with open(inputs["vars_json"], "r") as f:
            vars_dict = json.load(f)
    except Exception as e:
        raise Exception(
            "variable file would not import. It should be the \
            vars.json file written by AutoGrow in the output folder of the run."
        ) from e

    if "outfile_format" in inputs:
        if inputs["outfile_format"] is None:
            inputs["outfile_format"] = "svg"
        if inputs["outfile_format"].lower() not in ["svg", "png", "jpg", "pdf"]:
            raise Exception("outfile_format not a valid format")

    if "outfile" in inputs and inputs["outfile"] is not None:
        if os.path.dirname(inputs["outfile"]) is False:
            try:
                os.mkdir(os.path.dirname(inputs["outfile"]))
            except Exception as exc:
                raise Exception("outfile directory does not exist") from exc
        if os.path.dirname(inputs["outfile"]) is False:
            raise Exception("outfile directory does not exist")
    else:
        inputs["outfile"] = (
            inputs["infolder"] + os.sep + "data_line_plot." + inputs["outfile_format"]
        )

    # update --plot_reference_lines
    if "plot_reference_lines" not in inputs.keys():
        inputs["plot_reference_lines"] = None

    if inputs["plot_reference_lines"] is not None:
        _parse_plot_reference_lines(inputs)
    # overwrite and return vars_dict with input commands
    for key in inputs:
        vars_dict[key] = inputs[key]

    return vars_dict


def _parse_plot_reference_lines(inputs: Dict[str, Any]) -> None:
    # names of all matplotlib color options
    matplot_colors = matplotlib.colors.get_named_colors_mapping().keys()

    ref_lines = inputs["plot_reference_lines"].replace("[[", "[").replace("]]", "]")
    ref_lines = ref_lines.split("],")
    ref_lines = [ref.replace("]", "").replace("[", "").split(",") for ref in ref_lines]

    new_ref_lines = []
    failed_io = False
    for ref_info in ref_lines:
        if len(ref_info) != 3:
            failed_io = True
            break

        # make a new list with 1st item the str name
        temp_ref_lines: List[Union[str, float]] = [str(ref_info[0])]

        try:
            temp_ref_lines.append(float(ref_info[1]))
        except Exception:
            failed_io = True
            break
        if str(ref_info[2]) not in matplot_colors:
            print(f"COULD NOT FIND COLOR: {str(ref_info[2])}")
            failed_io = True
            break
        temp_ref_lines.append(str(ref_info[2]))
        new_ref_lines.append(temp_ref_lines)

    if failed_io is True:
        printout = (
            "\n --plot_reference_lines must be list of lists where each "
            + "sublist has three pieces of information in this "
        )
        printout += "order:\n\t [name, value, matplotlib_color]\n"
        printout += "more details can be found using the -h option\n"
        print(printout)
        raise Exception(printout)
    inputs["plot_reference_lines"] = new_ref_lines


def main(**kwargs):
    # copying ARGSDICT so we can delete out of while iterating through the
    # original ARGSDICT
    INPUTS = copy.deepcopy(kwargs)

    for k, v in kwargs.items():
        if v is None:
            del INPUTS[k]

    USER_VARS = process_inputs(INPUTS)

    generate_figures(USER_VARS, True)

    print(f'FINISHED {USER_VARS["outfile"]}')

    print("finished")


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()

    # Get needed info
    PARSER.add_argument(
        "--outfile",
        "-o",
        metavar="param.outfile",
        required=False,
        default=None,
        help="Path to folder to output files. It will be created if does not exist. \
            If not provide it will be placed in the infolder/data_line_plot.svg",
    )
    PARSER.add_argument(
        "--outfile_format",
        metavar="param.outfile_format",
        type=str,
        default="svg",
        choices=["svg", "png", "jpg", "pdf"],
        help="The type of file for figure to be exported as default is .svg file.",
    )
    PARSER.add_argument(
        "--infolder",
        "-i",
        metavar="param.infolder",
        required=True,
        help="Path to input folder containing the AutoGrow run. This should be the \
                top folder which contains the vars.json file.",
    )
    PARSER.add_argument(
        "--plot_reference_lines",
        default=None,
        help="This will be a list of lists, with each sublist being a different \
                dotted-line reference to plot. For each sublist the order of \
                information should be: [name, value, matplotlib_color] \
                For example a [['Olaparib score',-12.8,'y'],['Niraparib score',-10.7,'k']] \
                will add horizontal dotted lines at -12.8 (yellow) and -10.7 (black) \
                with Olaparib and Niraparib added to the legend. \
                Spaces must be within quotes and not be between variables. \
                matplotlib colors can be found with mcolors.get_named_colors_mapping().keys()",
    )

    kwargs = vars(PARSER.parse_args())
    main(**kwargs)
