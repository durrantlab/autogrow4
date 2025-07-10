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
from autogrow.utils.rank_file import get_gen_number_from_folder_name, load_rank_file
import prolif
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE


# ============================================================================
# HELPER FUNCTIONS FOR COMMON OPERATIONS
# ============================================================================

def get_generation_folders(infolder: str) -> List[str]:
    """
    Get all generation folders from the input directory.
    
    Args:
        infolder: Path to the folder containing all generation folders.
        
    Returns:
        List of generation folder paths.
    """
    generation_folders = []
    for gen_folder in os.listdir(infolder):
        if "generation" in gen_folder and os.path.isdir(os.path.join(infolder, gen_folder)):
            generation_folders.append(infolder + os.sep + gen_folder + os.sep)
    return generation_folders


def load_molecules_from_file(source_file: str) -> List:
    """
    Load molecules from either .smi or .sdf file.
    
    Args:
        source_file: Path to the source file.
        
    Returns:
        List of RDKit molecules.
    """
    if ".smi" in source_file:
        return read_smi_file(source_file)
    elif ".sdf" in source_file:
        return read_sdf_file(source_file)
    else:
        raise ValueError(f"Unsupported file format: {source_file}")


def calculate_ligand_efficiency(docking_score: float, smiles: str) -> float:
    """
    Calculate ligand efficiency from docking score and SMILES.
    
    Args:
        docking_score: Docking score value.
        smiles: SMILES string of the molecule.
        
    Returns:
        Ligand efficiency value.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    num_heavy_atoms = Lipinski.HeavyAtomCount(mol)
    return docking_score / num_heavy_atoms


def get_morgan_fingerprint(mol, radius: int = 10, use_features: bool = True):
    """
    Get Morgan fingerprint for a molecule.
    
    Args:
        mol: RDKit molecule object.
        radius: Fingerprint radius.
        use_features: Whether to use features.
        
    Returns:
        Morgan fingerprint.
    """
    if use_features:
        return GetMorganFingerprint(mol, radius=radius, useFeatures=True)
    else:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=2048)


def filter_valid_compounds(ranked_cmpds: List, require_docking_score: bool = False) -> List:
    """
    Filter out invalid compounds from ranked compounds list.
    
    Args:
        ranked_cmpds: List of ranked compounds.
        require_docking_score: Whether to require valid docking score.
        
    Returns:
        Filtered list of valid compounds.
    """
    if require_docking_score:
        return [c for c in ranked_cmpds if c.docking_score is not None]
    else:
        return [c for c in ranked_cmpds if c.smiles is not None]


def setup_plot_styling():
    """Set up common plot styling."""
    plt.clf()


def save_plot(outfile: str, params: Dict[str, Any]):
    """
    Save plot with common parameters.
    
    Args:
        outfile: Output file path.
        params: Parameters dictionary containing format and DPI settings.
    """
    plt.savefig(outfile, bbox_inches="tight", format=params["outfile_format"], dpi=1000)


def get_plot_labels(params: Dict[str, Any], ligand_efficiency: bool = False) -> Tuple[str, str, str]:
    """
    Get standard plot labels based on parameters.
    
    Args:
        params: Parameters dictionary.
        ligand_efficiency: Whether calculating ligand efficiency.
        
    Returns:
        Tuple of (title, x_label, y_label).
    """
    receptor_name = os.path.basename(params["receptor_path"])
    
    if ligand_efficiency:
        title = f"Ligand efficiencies for {receptor_name}"
        y_label = "Ligand Efficiency"
    else:
        title = f"Docking scores for {receptor_name}"
        scoring_type = params.get("vina_like_executable", "")
        if "vina" in str(scoring_type):
            y_label = "Docking Score"
        else:
            y_label = "Fitness Score"
    
    return title, "Generation Number", y_label


# ============================================================================
# CORE ANALYSIS FUNCTIONS
# ============================================================================

def exist_generation_0(infolder: str) -> bool:
    """
    Returns if generation 0 exists.

    Args:
        infolder: Path to the folder containing all generation folders.

    Returns:
        True if generation 0 exists.
    """
    generation_folders = get_generation_folders(infolder)
    return any("generation_0_" in folder for folder in generation_folders)


def read_smi_file(source_file: str):
    """Read molecules from SMILES file."""
    rdkit_molecules = []
    with open(source_file) as smiles_file:
        for tsv_line in smiles_file:
            prts = tsv_line.replace("    ", "\t").strip().split("\t")
            rdkit_molecules.append(Chem.MolFromSmiles(prts[0]))

    if None in rdkit_molecules:
        raise Exception("An input molecule was not successfully read from " + source_file)

    return rdkit_molecules


def read_sdf_file(source_file: str):
    """Read molecules from SDF file."""
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
        List[float]: similarity values.

    Note:
        - Lower diversity scores indicate more diverse molecules.
        - Removes any None entries from the input list.
        - Uses Morgan Fingerprints with radius 10 and feature-based encoding.
    """
    similarity_values = []
    reference_comp_fp = get_morgan_fingerprint(reference_comp)
    
    for mol in new_comps:
        if mol is not None:
            mol_fp = get_morgan_fingerprint(mol)

            # if DiceSimilarity=1.0 it is a perfect match, the smaller the
            # number, the more diverse it is.
            diversity_score = DataStructs.DiceSimilarity(reference_comp_fp, mol_fp)
            similarity_values.append(diversity_score)

    return similarity_values


def calc_interaction_fp_per_generation(params: Dict[str, Any], infolder: str, 
                                     interactions: List[str], analyze_gen_0: bool):
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
    
    generation_folders = get_generation_folders(infolder)
    
    for gen_folder_name in generation_folders:
        ranked_cmpds = load_rank_file(gen_folder_name)
        gen_num = get_gen_number_from_folder_name(gen_folder_name)
        num_compounds_per_generation[gen_num] = 0

        for ranked_cmpd in ranked_cmpds:
            if ranked_cmpd.sdf_path == "None":
                continue

            num_compounds_per_generation[gen_num] += 1
            docked_comp = read_sdf_file(ranked_cmpd.sdf_path)[0]
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
                number_interactions_per_gen[interaction][gen_num] += 1

            # rank_file = glob.glob(f"{gen_folder_name}*_ranked.smi")[0]
            # gen_num = os.path.basename(rank_file).split("_")[1]
            # num_compounds_per_generation[gen_num] = 0

            # with open(rank_file, "r") as f:
            #     for line in f:
            #         line = line.replace("\n", "")

            #         # split line into parts separated by 5-spaces
            #         parts = line.split("\t")

            #         docked_comp_path = [parts[i] for i in range(len(parts))][5]
            #         if docked_comp_path != "None":
            #             num_compounds_per_generation[gen_num] += 1
            #             docked_comp = read_sdf_file(docked_comp_path)[0]
            #             prot = prolif.Molecule.from_rdkit(receptor)
            #             lig = prolif.Molecule.from_rdkit(docked_comp)
            #             ifp = fingerprints.generate(lig, prot, metadata=True)
            #             df = prolif.to_dataframe({0: ifp}, fingerprints.interactions)

            #             for column_name in df.columns:
            #                 interaction = column_name[1].split(".")[0] + "_" + column_name[2]
            #                 if interaction not in number_interactions_per_gen:
            #                     number_interactions_per_gen[interaction] = {}
            #                     all_interactions.append(interaction)
            #                 if gen_num not in number_interactions_per_gen[interaction]:
            #                     number_interactions_per_gen[interaction][gen_num] = 0
            #                 number_interactions_per_gen[interaction][gen_num] = \
            #                     number_interactions_per_gen[interaction][gen_num] + 1

    all_generations = []
    result_matrix = np.zeros((len(num_compounds_per_generation), len(all_interactions)))
    
    for gen_num_idx in range(len(num_compounds_per_generation)):
        gen_num = (gen_num_idx + 1) if not analyze_gen_0 else gen_num_idx
        all_generations.append(f"generation_{gen_num}")
        
        for interaction_num in range(len(all_interactions)):
            interaction = all_interactions[interaction_num]
            if str(gen_num) in number_interactions_per_gen[interaction]:
                percent_of_interactions = number_interactions_per_gen[interaction][str(gen_num)]
                percent_of_interactions = round(
                    float(percent_of_interactions / num_compounds_per_generation[str(gen_num)]), 2
                )
                result_matrix[gen_num_idx, interaction_num] = percent_of_interactions

    return result_matrix, all_generations, all_interactions


def generate_tSNE_scatterplot(infolder: str, params: Dict[str, Any], outfile: str, exist_gen_0: bool):
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
    
    generation_folders = get_generation_folders(infolder)
    
    for gen_folder_name in generation_folders:
        ranked_cmpds = load_rank_file(gen_folder_name)
        new_comps = [Chem.MolFromSmiles(cmpd.smiles) for cmpd in ranked_cmpds]
        gen_num = get_gen_number_from_folder_name(gen_folder_name)

        for mol in new_comps:
            if mol is None:
                continue
            mol_fp = get_morgan_fingerprint(mol, use_features=False)
            result_matrix.append(list(mol_fp))
            label_list.append(f"generation_{str(gen_num)}")

    # Generate t-SNE plot
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

    setup_plot_styling()
    fig, ax = plt.subplots()
    
    for gen_id in range(len(x)):
        gen_name = "generation_" + (str(gen_id + 1) if not exist_gen_0 else str(gen_id))
        x_gen = x.loc[gen_name]
        y_gen = y.loc[gen_name]
        ax.scatter(x_gen, y_gen, label=gen_name, c='black' if gen_id == 0 else None)

    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

    # Add titles and labels
    plt.xlabel("Dimension 1", fontweight="semibold")
    plt.ylabel("Dimension 2", fontweight="semibold")

    save_plot(outfile, params)


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
    generation_folders = get_generation_folders(infolder)
    
    for gen_folder_name in generation_folders:
        ranked_cmpds = load_rank_file(gen_folder_name)
        gen_num = get_gen_number_from_folder_name(gen_folder_name)
        new_comps = [Chem.MolFromSmiles(cmpd.smiles) for cmpd in ranked_cmpds]

        similarity_values = []
        for i in range(len(new_comps)):
            reference_comp = new_comps[i]
            if reference_comp is None:
                continue

            reference_comp_fp = get_morgan_fingerprint(reference_comp)
            for j in range(i + 1, len(new_comps)):
                mol = new_comps[j]
                if mol is None:
                    continue

                mol_fp = get_morgan_fingerprint(mol)
                diversity_score = DataStructs.DiceSimilarity(reference_comp_fp, mol_fp)
                
                if diversity_score >= 0.9:
                    print(f"Compounds {i+1} and {j+1} of generation {gen_num} are similar: {diversity_score}")
                
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
    input_rdkit_molecules = load_molecules_from_file(source_file)
    
    ranked_file = infolder + os.sep + "summary_ranked.sdf"
    ranked_molecules = read_sdf_file(ranked_file)

    for i, reference_mol in enumerate(input_rdkit_molecules):
        similarity_dict[f"compound_{i+1}"] = calc_diversity_scores(reference_mol, ranked_molecules)

    return similarity_dict


def get_ave_similarity_per_generated_comp(infolder: str, source_file: str) -> Dict[str, float]:
    """
    Calculate average similarity between generated compounds and input compounds.

    Args:
        infolder (str): Path to the folder containing all generation folders.
        source_file (str): Path to the input files containing the compounds.

    Returns:
        Dict[str, float]: Dictionary of average similarity values, with compound ID as keys.

    Notes:
        if a generated molecule is not successfully read, a similarity equal to 1 is assigned.
    """
    average_similarity_dict = {}
    input_rdkit_molecules = load_molecules_from_file(source_file)
    
    ranked_file = infolder + os.sep + "summary_ranked.sdf"
    ranked_molecules = read_sdf_file(ranked_file)

    for i, reference_mol in enumerate(ranked_molecules):
        if reference_mol is None:
            average_similarity_dict[f"compound_{i+1}"] = -0.05
        else:
            similarity_values = calc_diversity_scores(reference_mol, input_rdkit_molecules)
            average_similarity_dict[f"compound_{i+1}"] = sum(similarity_values) / len(input_rdkit_molecules)

    return average_similarity_dict


def get_efficiency_per_generated_comp(infolder: str) -> Dict[str, float]:
    """
    Calculate ligand efficiency for generated compounds.

    Args:
        infolder (str): Path to the folder containing all generation folders.

    Returns:
        Dictionary of ligand efficiency values with compound ID as keys.

    Notes:
        if a generated molecule is not successfully read, an efficiency equal to 0 is assigned.
    """
    efficiency_dict = {}
    ranked_file = infolder + os.sep + "summary_ranked.sdf"
    ranked_molecules = read_sdf_file(ranked_file)

    for i, mol in enumerate(ranked_molecules):
        if mol is None:
            efficiency_dict[f"compound_{i+1}"] = 0.0
        else:
            docking_score = float(mol.GetProp('Docking_Score'))
            num_heavy_atoms = Lipinski.HeavyAtomCount(mol)
            efficiency_dict[f"compound_{i+1}"] = docking_score / num_heavy_atoms

    return efficiency_dict


def get_scores_per_generation(infolder: str, ligand_efficiency: bool = False) -> Dict[str, List[float]]:
    """
    Get scores for each generation.

    This function reads ranked .smi files from each generation folder and
    gets the docking scores (or ligand efficiencies) for each generation.

    Args:
    :param str infolder: Path to the folder containing all generation folders.
    :param bool ligand_efficiency: If True, the ligand efficiency will be calculated instead
        of the affinity score.

    Returns:
        Dict[str, List[float]]: Dictionary of lists of docking scores (or ligand efficiencies)
        for each generation, with generation names as keys.
    """
    score_dict = {}
    generation_folders = get_generation_folders(infolder)
    
    for gen_folder_name in generation_folders:
        ranked_cmpds = load_rank_file(gen_folder_name)
        gen_num = get_gen_number_from_folder_name(gen_folder_name)
        gen_list = []

        for cmpd in ranked_cmpds:
            if ligand_efficiency:
                score = calculate_ligand_efficiency(cmpd.docking_score, cmpd.smiles)
            else:
                score = cmpd.docking_score
            gen_list.append(score)

        gen_name = f"generation_{gen_num}"
        score_dict[gen_name] = gen_list

        # ranked_file = glob.glob(f"{gen_folder_name}*_ranked.smi")

        # for rank_file in ranked_file:
        #     with open(rank_file, "r") as f:
        #         for line in f:
        #             line = line.replace("\n", "")
        #             parts = line.split(
        #                 "\t"
        #             )  # split line into parts separated by 4-spaces

        #             choice_list = [parts[i] for i in range(len(parts))]

        #             docking_score = float(choice_list[2])
        #             num_heavy_atoms = Lipinski.HeavyAtomCount(Chem.MolFromSmiles(choice_list[0], sanitize=False))
        #             ligand_efficiency_value = docking_score / num_heavy_atoms

        #             gen_list.append(docking_score if not ligand_efficiency else ligand_efficiency_value)

        #     gen_num = os.path.basename(rank_file).split("_")[1]
        #     gen_name = f"generation_{gen_num}"
        #     average_dict[gen_name] = gen_list

    return score_dict


def get_average_score_per_gen(infolder: str, ligand_efficiency: bool = False) -> Dict[str, float]:
    """
    Calculate the average docking score (or ligand efficiency) for each generation.

    This function reads ranked .smi files from each generation folder and
    calculates the average docking score (or ligand efficiency) for each generation.

    Args:
    :param str infolder: Path to the folder containing all generation folders.
    :param bool ligand_efficiency: If True, the ligand efficiency will be calculated instead
        of the average affinity score.

    Returns:
        Dict[str, List[float]]: Dictionary of average docking scores (or ligand efficiencies)
        for each generation, with generation names as keys.
    """
    average_dict = {}
    generation_folders = get_generation_folders(infolder)
    
    for gen_folder_name in generation_folders:
        ranked_cmpds = load_rank_file(gen_folder_name)
        ranked_cmpds = filter_valid_compounds(ranked_cmpds, require_docking_score=True)
        
        if not ranked_cmpds:
            continue

        gen_num = get_gen_number_from_folder_name(gen_folder_name)
        total_score = 0.0
        
        for cmpd in ranked_cmpds:
            if ligand_efficiency:
                score = calculate_ligand_efficiency(cmpd.docking_score, cmpd.smiles)
            else:
                score = cmpd.docking_score
            total_score += score
            
        average_dict[f"generation_{gen_num}"] = total_score / len(ranked_cmpds)

            # TODO: Several of these old versions account for multiple ranked
            # files per generation. Does this actually happen? Good to check
            # with Cesar.

            # for rank_file in ranked_file:
            #     # write as a tab delineated .smi file
            #     with open(rank_file, "r") as f:
            #         gen_affinity_sum = 0.0
            #         num_lines_counter = 0.0
            #         for line in f:
            #             line = line.replace("\n", "")
            #             parts = line.split(
            #                 "\t"
            #             )  # split line into parts separated by 4-spaces

            #             choice_list = [parts[i] for i in range(len(parts))]

            #             docking_score = float(choice_list[2])
            #             num_heavy_atoms = Lipinski.HeavyAtomCount(Chem.MolFromSmiles(choice_list[0], sanitize=False))
            #             ligand_efficiency_value = docking_score / num_heavy_atoms

            #             gen_affinity_sum = gen_affinity_sum + docking_score \
            #                 if not ligand_efficiency \
            #                 else gen_affinity_sum + ligand_efficiency_value
            #             num_lines_counter = num_lines_counter + 1.0

            #     gen_affinity_average = gen_affinity_sum / num_lines_counter

            #     gen_num = os.path.basename(rank_file).split("_")[1]
            #     gen_name = f"generation_{gen_num}"
            #     average_affinity_dict[gen_name] = gen_affinity_average

    print_generation_scores(average_dict)
    return average_dict


def get_average_top_score_per_gen(infolder: str, top_n: int, ligand_efficiency: bool = False) -> Dict[str, Union[float, str]]:
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
    average_dict = {}
    generation_folders = get_generation_folders(infolder)
    
    for gen_folder_name in generation_folders:
        ranked_cmpds = load_rank_file(gen_folder_name)
        ranked_cmpds = filter_valid_compounds(ranked_cmpds, require_docking_score=True)
        gen_num = get_gen_number_from_folder_name(gen_folder_name)
        
        # Sort compounds by score
        if ligand_efficiency:
            ranked_cmpds.sort(key=lambda c: calculate_ligand_efficiency(c.docking_score, c.smiles))
        else:
            ranked_cmpds.sort(key=lambda c: c.docking_score)

        if len(ranked_cmpds) < top_n:
            average_dict[f"generation_{gen_num}"] = "N/A"
            continue

        total_score = 0.0
        for cmpd in ranked_cmpds[:top_n]:
            if ligand_efficiency:
                score = calculate_ligand_efficiency(cmpd.docking_score, cmpd.smiles)
            else:
                score = cmpd.docking_score
            total_score += score

        average_dict[f"generation_{gen_num}"] = total_score / top_n

        #             gen_affinity_average = gen_affinity_sum / top_score_per_gen

        #             gen_num = os.path.basename(rank_file).split("_")[1]
        #             gen_name = f"generation_{gen_num}"
        #             average_affinity_dict[gen_name] = gen_affinity_average

        #     else:
        #         gen_num = os.path.basename(rank_file).split("_")[1]
        #         gen_name = f"generation_{gen_num}"
        #         average_affinity_dict[gen_name] = "N/A"

    print_generation_scores(average_dict)
    return average_dict


def print_generation_scores(score_dict: Dict[str, Union[float, str]]) -> None:
    """
    This prints out the average scores for each generation

    Inputs:
    :param dict average_affinity_dict: dictionary of average values
        for top_score_per_gen number of ligands
    """
    print("generation_number              average")
    sorted_keys = sorted(score_dict.keys(), key=lambda x: int(x.split("_")[1]))
    for gen in sorted_keys:
        print(f"{gen}                  {score_dict[gen]}")


def make_graph(dictionary: Dict[str, Union[float, str]], analyze_gen_0: bool, 
               exist_gen_0: bool) -> Tuple[Optional[List[int]], Optional[List[float]]]:
    """
    Because some generations may not have 50 ligands this basically checks to see if
    there are enough ligands and prepares lists to be plotted

    Inputs:
    :param dict dictionary: dictionary of average affinity scores for
        top_score_per_gen number of ligands
    :param bool analyze_gen_0: if True, the compounds of the generation 0, which are the input compounds, will be
    analyzed.
    :param bool exist_gen_0: If True, the generation exists.

    Returns:
    :returns: list list_generations: list of ints for each generation to be plotted.
        if a generation lacks ligands to generate the average it will return "N/A"
    :returns: list list_of_scores: list of averages for each generation;
        if a generation lacks ligands to generate the average it will return "N/A"
    """
    if any(val == "N/A" for val in dictionary.values()):
        return None, None

    sorted_keys = sorted(dictionary.keys(), key=lambda k: int(k.split("_")[1]))

    if not analyze_gen_0 and exist_gen_0 and "generation_0" in sorted_keys:
        sorted_keys.remove("generation_0")

    if not sorted_keys:
        return [], []

    list_generations = [int(k.split("_")[1]) for k in sorted_keys]
    list_of_scores = [dictionary[k] for k in sorted_keys]

    if "N/A" in list_of_scores:
        return None, None

    return list_generations, list_of_scores


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def run_score_plotter(params: Dict[str, Any], dict_of_averages: Dict[str, Dict[str, Union[float, str]]],
                     outfile: str, ligand_efficiency: bool, analyze_gen_0: bool, exist_gen_0: bool) -> None:
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
    :param bool exist_gen_0: If True, the generation exists.
    """
    setup_plot_styling()
    ax = plt.subplot(111)

    # Plot different top-N averages
    plot_configs = [
        ("average_dict", "Average", "b"),
        ("top_15", "Top 15", "c"),
        ("top_10", "Top 10", "m"),
        ("top_5", "Top 5", "g"),
        ("top_1", "Top 1", "r")
    ]

    for key, label, color in plot_configs:
        if key in dict_of_averages:
            generations, scores = make_graph(dict_of_averages[key], analyze_gen_0, exist_gen_0)
            if generations is not None and scores is not None:
                ax.plot(generations, scores, color=color, label=label)

    # Add reference lines if specified
    if params.get("plot_reference_lines"):
        for ref_info in params["plot_reference_lines"]:
            ax.axhline(y=ref_info[1], color=ref_info[2], linestyle=":", label=ref_info[0])

    # Configure plot
    title, x_label, y_label = get_plot_labels(params, ligand_efficiency)
    plt.title(title, fontweight="semibold")
    plt.xlabel(x_label, fontweight="semibold")
    plt.ylabel(y_label, fontweight="semibold")
    
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.274), fontsize="small")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    save_plot(outfile, params)


def run_boxplot(params: Dict[str, Any], dictionary_of_values: Dict[str, List[float]], outfile: str,
               key_start_with: str, x_label: str, y_label: str, title_of_figure: str = None,
               analyze_gen_0: bool = False, exist_gen_0: bool = False) -> None:
    """
    Create boxplot for the given data.

    Args:
        params: Parameters dictionary.
        dictionary_of_values: Dictionary of value lists.
        outfile: Output file path.
        key_start_with: Key prefix for data selection.
        x_label: X-axis label.
        y_label: Y-axis label.
        title_of_figure: Plot title.
        analyze_gen_0: Whether to analyze generation 0.
        exist_gen_0: Whether generation 0 exists.
    """
    data = []
    yticklabels = []
    
    num_generations = len(dictionary_of_values) - (1 if exist_gen_0 and not analyze_gen_0 else 0)
    
    for i in range(num_generations):
        key = f"{key_start_with}_{i + 1 if not analyze_gen_0 else i}"
        if key in dictionary_of_values:
            data.append(dictionary_of_values[key])
            yticklabels.append(key)

    setup_plot_styling()
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)

    bp = ax.boxplot(data, patch_artist=True, notch=False, vert=0, showmeans=True,
                    meanprops={"markerfacecolor": "black", "markeredgecolor": "black"})
    
    for median in bp['medians']:
        median.set_color('black')

    # Setting x-axis labels

    plt.xlabel(x_label, fontweight="semibold")
    plt.ylabel(y_label, fontweight="semibold")
    ax.set_yticklabels(yticklabels)
    
    if title_of_figure:
        plt.title(title_of_figure, fontweight="semibold")

    save_plot(outfile, params)


def run_plotter(params: Dict[str, Any], dictionary_of_values: Dict[str, Any], outfile: str,
               key_start_with: str, x_label: str, y_label: str, title_of_figure: str = None) -> None:
    """
    Create line plot for the given data.

    Args:
        params: Parameters dictionary.
        dictionary_of_values: Dictionary of values.
        outfile: Output file path.
        key_start_with: Key prefix for data selection.
        x_label: X-axis label.
        y_label: Y-axis label.
        title_of_figure: Plot title.
    """
    x = []
    y = []
    
    for i in range(len(dictionary_of_values)):
        key = f"{key_start_with}_{i + 1}"
        if key in dictionary_of_values:
            y.append(dictionary_of_values[key])
            x.append(i + 1)

    setup_plot_styling()
    plt.plot(x, y, marker='o', linestyle='-', color='b')

    if title_of_figure:
        plt.title(title_of_figure, fontweight="semibold")
    if x_label:
        plt.xlabel(x_label, fontweight="semibold")
    if y_label:
        plt.ylabel(y_label, fontweight="semibold")

    plt.legend()
    save_plot(outfile, params)


def run_heatmap(params: Dict[str, Any], outfile: str, result_matrix, x_label_list, y_label_list,
               title_of_figure: str) -> None:
    """
    Create heatmap for the given data.

    Args:
        params: Parameters dictionary.
        outfile: Output file path.
        result_matrix: Data matrix for heatmap.
        x_label_list: X-axis labels.
        y_label_list: Y-axis labels.
        title_of_figure: Plot title.
    """
    setup_plot_styling()
    fig, ax = plt.subplots()
    ax.imshow(result_matrix, cmap="Blues")

    ax.set_xticks(range(len(x_label_list)), labels=x_label_list, rotation=45, ha="right", 
                  rotation_mode="anchor", size="x-small")
    ax.set_yticks(range(len(y_label_list)), labels=y_label_list, size="x-small")

    # Add text annotations
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
    save_plot(outfile, params)


# ============================================================================
# MAIN ANALYSIS FUNCTIONS
# ============================================================================

def print_data_table(infolder: str, ligand_efficiency: bool) -> Dict[str, Any]:
    """
    Print data table and return analysis results.

    Args:
        infolder: Path to the folder containing all generation folders.
        ligand_efficiency: Whether to calculate ligand efficiency.

    Returns
    :returns: dict dict_of_averages: a dictionary of dictionaries containing
        the average of each generation for the top 50,20, 10, and 1 ligand(s) and
        the overall average for each generation.
    """
    score_type = "Ligand Efficiency" if ligand_efficiency else "Scoring"
    
    print(f"Overall {score_type} Average for all Compounds")
    average_affinity_dict = get_average_score_per_gen(infolder, ligand_efficiency)
    print("")
    
    # Calculate top N averages
    top_results = {}
    for n in [15, 10, 5, 1]:
        print(f"Average for Top {score_type} Compounds")
        print(f"Number of top compounds: {n}")
        top_results[f"top_{n}"] = get_average_top_score_per_gen(infolder, n, ligand_efficiency)
        print("")
    
    return {
        "average_dict": average_affinity_dict,
        **top_results
    }


def generate_all_plots(params: Dict[str, Any], infolder: str, outfile: str, 
                      analyze_gen_0: bool, exist_gen_0: bool) -> None:
    """
    Generate all plots for the analysis.

    Args:
        params: Parameters dictionary.
        infolder: Input folder path.
        outfile: Output folder path.
        analyze_gen_0: Whether to analyze generation 0.
        exist_gen_0: Whether generation 0 exists.
    """
    # Score plots
    dict_of_averages = print_data_table(infolder, False)
    run_score_plotter(params, dict_of_averages,
                      f"{outfile}{os.sep}plotter_by_generation_for_scores.{params['outfile_format']}",
                      ligand_efficiency=False, analyze_gen_0=analyze_gen_0, exist_gen_0=exist_gen_0)

    # Ligand efficiency plots
    dict_of_averages = print_data_table(infolder, True)
    run_score_plotter(params, dict_of_averages,
                      f"{outfile}{os.sep}plotter_by_generation_for_ligand_efficiencies.{params['outfile_format']}",
                      ligand_efficiency=True, analyze_gen_0=analyze_gen_0, exist_gen_0=exist_gen_0)

    # Boxplots for scores
    dict_of_score_lists = get_scores_per_generation(infolder, False)
    title, _, y_label = get_plot_labels(params, False)
    run_boxplot(params, dict_of_score_lists,
                f"{outfile}{os.sep}boxplot_by_generation_for_scores.{params['outfile_format']}",
                key_start_with="generation", x_label=y_label, y_label="Number of Generations",
                title_of_figure=title, analyze_gen_0=analyze_gen_0, exist_gen_0=exist_gen_0)

    # Boxplots for ligand efficiencies
    dict_of_score_lists = get_scores_per_generation(infolder, True)
    title, _, y_label = get_plot_labels(params, True)
    run_boxplot(params, dict_of_score_lists,
                f"{outfile}{os.sep}boxplot_by_generation_for_ligand_efficiencies.{params['outfile_format']}",
                key_start_with="generation", x_label=y_label, y_label="Number of Generations",
                title_of_figure=title, analyze_gen_0=analyze_gen_0, exist_gen_0=exist_gen_0)

    # Similarity analysis (only if generation 0 exists)
    source_file = str(params["source_compound_file"])
    if exist_gen_0:
        dict_of_similarity_lists = get_similarity_list_per_input_comp(infolder, source_file)
        run_boxplot(params, dict_of_similarity_lists,
                    f"{outfile}{os.sep}boxplot_of_similarity_between_input_and_new_compounds.{params['outfile_format']}",
                    key_start_with="compound", x_label="Dice Similarity Values", y_label="Input Compounds")

    # Diversity analysis
    dict_of_score_lists = calc_diversity_scores_per_generation(infolder)
    run_boxplot(params, dict_of_score_lists,
                f"{outfile}{os.sep}boxplot_by_generation_for_similarities.{params['outfile_format']}",
                key_start_with="generation", x_label="Dice Similarity Values", y_label="Number of Generations",
                analyze_gen_0=analyze_gen_0, exist_gen_0=exist_gen_0)

    # Additional plots
    if exist_gen_0:
        dict_of_similarity_lists = get_ave_similarity_per_generated_comp(infolder, source_file)
        run_plotter(params, dict_of_similarity_lists,
                    f"{outfile}{os.sep}plotter_of_ave_similarities_between_each_new_compound_and_input_compounds.{params['outfile_format']}",
                    key_start_with="compound",
                    x_label="New compounds (ID) sorted from the highest to lowest affinity",
                    y_label="Average similarity")

    dict_of_efficiency_lists = get_efficiency_per_generated_comp(infolder)
    run_plotter(params, dict_of_efficiency_lists,
                f"{outfile}{os.sep}plotter_of_ligand_efficiency_for_every_new_compound.{params['outfile_format']}",
                key_start_with="compound",
                x_label="New compounds (ID) sorted from the highest to lowest affinity",
                y_label="Ligand efficiency")

    # t-SNE plot
    try:
        generate_tSNE_scatterplot(infolder=infolder,
                                  outfile=f"{outfile}{os.sep}tsne_for_input_and_new_compounds.{params['outfile_format']}",
                                  params=params, exist_gen_0=exist_gen_0)
    except:
        pass

    # Interaction analysis
    generate_interaction_analysis(params, infolder, outfile, analyze_gen_0)


def generate_interaction_analysis(params: Dict[str, Any], infolder: str, outfile: str, analyze_gen_0: bool) -> None:
    """
    Generate interaction analysis plots.

    Args:
        params: Parameters dictionary.
        infolder: Input folder path.
        outfile: Output folder path.
        analyze_gen_0: Whether to analyze generation 0.
    """
    interactions = ["Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "Anionic", "Cationic", "CationPi", "PiCation"]
    result_matrix, y_label_list, x_label_list = calc_interaction_fp_per_generation(params, infolder, interactions, analyze_gen_0)
    
    # All interactions heatmap
    run_heatmap(params,
                f"{outfile}{os.sep}heatmap_of_all_interactions_by_generation.{params['outfile_format']}",
                result_matrix, x_label_list, y_label_list, "Analysis of Interactions")

    # Individual interaction heatmaps
    for interaction in interactions:
        pos_for_interaction = [i for i in range(len(x_label_list)) if interaction in x_label_list[i]]
        if pos_for_interaction:
            result_matrix_c = result_matrix[:, pos_for_interaction]
            x_label_list_c = [x_label_list[i] for i in pos_for_interaction]
            run_heatmap(params,
                        f"{outfile}{os.sep}heatmap_of_{interaction}_interactions_by_generation.{params['outfile_format']}",
                        result_matrix_c, x_label_list_c, y_label_list, f"Analysis of {interaction} Interactions")


def generate_figures(params: Dict[str, Any], analyze_gen_0: bool) -> None:
    """
    Main function to generate all figures for the analysis.

    Args:
        params: Parameters dictionary.
        analyze_gen_0: Whether to analyze generation 0.
    """
    infolder = params["infolder"]
    outfile = params["outfile"]
    exist_gen_0 = exist_generation_0(infolder)
    analyze_gen_0 = analyze_gen_0 if exist_gen_0 else False

    generate_all_plots(params, infolder, outfile, analyze_gen_0, exist_gen_0)


# ============================================================================
# CONFIGURATION AND UTILITY FUNCTIONS
# ============================================================================

def retrieve_vars_dict(autogrow_vars_json: str) -> Dict[str, Any]:
    """
    Retrieve variable dictionary from AutoGrow vars json file.

    Args:
        autogrow_vars_json: Path to AutoGrow json variable file.

    Returns:
        Dictionary of variables.
    """
    if not os.path.exists(autogrow_vars_json):
        raise Exception("Variable file could not be found. It should be the vars.json file written by AutoGrow.")
    
    try:
        with open(autogrow_vars_json, "r") as f:
            params = json.load(f)
    except Exception as e:
        raise Exception("Variable file would not import. It should be the vars.json file written by AutoGrow.") from e
    
    return params


def validate_input_folder(infolder: str) -> str:
    """
    Validate and process input folder path.

    Args:
        infolder: Input folder path.

    Returns:
        Absolute path to input folder.
    """
    infolder = os.path.abspath(infolder) + os.sep
    if not os.path.exists(infolder):
        raise Exception(f"Input folder {infolder} does not exist.")
    
    vars_json = infolder + "vars.json"
    if not os.path.exists(vars_json):
        raise Exception(f"Input folder {infolder} does not contain the vars.json file necessary to run script.")
    
    return infolder


def validate_output_settings(inputs: Dict[str, Any]) -> None:
    """
    Validate and process output settings.

    Args:
        inputs: Input parameters dictionary.
    """
    if "outfile_format" in inputs:
        if inputs["outfile_format"] is None:
            inputs["outfile_format"] = "svg"
        if inputs["outfile_format"].lower() not in ["svg", "png", "jpg", "pdf"]:
            raise Exception("outfile_format not a valid format")

    if "outfile" in inputs and inputs["outfile"] is not None:
        if not os.path.exists(os.path.dirname(inputs["outfile"])):
            try:
                os.makedirs(os.path.dirname(inputs["outfile"]))
            except Exception as exc:
                raise Exception("outfile directory does not exist") from exc
    else:
        inputs["outfile"] = inputs["infolder"] + os.sep + "data_line_plot." + inputs["outfile_format"]


def parse_plot_reference_lines(inputs: Dict[str, Any]) -> None:
    """
    Parse and validate plot reference lines.

    Args:
        inputs: Input parameters dictionary.
    """
    if inputs.get("plot_reference_lines") is None:
        return

    matplot_colors = matplotlib.colors.get_named_colors_mapping().keys()
    ref_lines = inputs["plot_reference_lines"].replace("[[", "[").replace("]]", "]")
    ref_lines = ref_lines.split("],")
    ref_lines = [ref.replace("]", "").replace("[", "").split(",") for ref in ref_lines]

    new_ref_lines = []
    for ref_info in ref_lines:
        if len(ref_info) != 3:
            raise Exception("Invalid reference line format")

        temp_ref_lines = [str(ref_info[0])]
        
        try:
            temp_ref_lines.append(float(ref_info[1]))
        except Exception:
            raise Exception("Invalid reference line value")
        
        if str(ref_info[2]) not in matplot_colors:
            raise Exception(f"Could not find color: {str(ref_info[2])}")
        
        temp_ref_lines.append(str(ref_info[2]))
        new_ref_lines.append(temp_ref_lines)

    inputs["plot_reference_lines"] = new_ref_lines


def process_inputs(inputs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process and validate all input parameters.

    Args:
        inputs: Dictionary of input parameters.

    Returns:
        Dictionary of processed parameters.
    """
    # Validate input folder
    inputs["infolder"] = validate_input_folder(inputs["infolder"])
    
    # Get vars dict from last run
    inputs["vars_json"] = inputs["infolder"] + "vars.json"
    vars_dict = retrieve_vars_dict(inputs["vars_json"])
    
    # Validate output settings
    validate_output_settings(inputs)
    
    # Parse reference lines
    parse_plot_reference_lines(inputs)
    
    # Update vars_dict with input commands
    for key in inputs:
        vars_dict[key] = inputs[key]

    return vars_dict


def main(**kwargs):
    """
    Main function to run the analysis.

    Args:
        **kwargs: Command line arguments.
    """
    # Process inputs
    inputs = copy.deepcopy(kwargs)
    for k, v in kwargs.items():
        if v is None:
            del inputs[k]

    user_vars = process_inputs(inputs)
    
    # Generate figures
    generate_figures(user_vars, bool(user_vars.get("process_input_compounds", False)))

    print(f'FINISHED {user_vars["outfile"]}')
    print("finished")


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()

    # Command line arguments
    PARSER.add_argument(
        "--outfile", "-o", metavar="param.outfile", required=False, default=None,
        help="Path to folder to output files. If not provided, will be placed in infolder/data_line_plot.svg"
    )
    PARSER.add_argument(
        "--outfile_format", metavar="param.outfile_format", type=str, default="svg",
        choices=["svg", "png", "jpg", "pdf"],
        help="The type of file for figure to be exported as default is .svg file."
    )
    PARSER.add_argument(
        "--infolder", "-i", metavar="param.infolder", required=True,
        help="Path to input folder containing the AutoGrow run with vars.json file."
    )
    PARSER.add_argument(
        "--plot_reference_lines", default=None,
        help="List of reference lines to plot. Format: [['name',value,'color'],...]"
    )
    PARSER.add_argument(
        "--process_input_compounds", action="store_true", default=False,
        help="Use information of reference compounds in processing results"
    )

    kwargs = vars(PARSER.parse_args())
    main(**kwargs)