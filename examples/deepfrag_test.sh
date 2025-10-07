#!/bin/bash

# Clean up previous runs
rm -rf with_deepfrag_out
rm -f log.txt

# ============================================================================
# OUTPUT AND INPUT FILES
# ============================================================================
output_directory="./with_deepfrag_out"
receptor_path="./input_files/7KK4-molmoda-protonated.pdb"
source_compound_file="./input_files/7KK4_frag_modified_final.sdf"
rxn_library_path="../autogrow/plugins/mutation/reaction_libraries/all_rxns"

# ============================================================================
# EXTERNAL SOFTWARE PATHS AND EXECUTABLES
# Paths to required external programs
# ============================================================================
obabel_path="/Users/jdurrant/miniconda3/envs/openbabel-env/bin/obabel"
docking_executable="/Applications/vina_1.2.6_mac_aarch64"

# ============================================================================
# DOCKING BOX CONFIGURATION
# Define the 3D search space for molecular docking
# ============================================================================
center_x=-6.382
center_y=4.374
center_z=26.197
size_x=20.0 
size_y=20.0
size_z=20.0

# ============================================================================
# GENETIC ALGORITHM - FIRST GENERATION PARAMETERS
# Special settings for the initial generation of molecules
# ============================================================================
number_of_mutants_first_generation=20
number_of_crossovers_first_generation=0
top_mols_to_seed_next_generation_first_generation=0
number_elitism_advance_from_previous_gen_first_generation=0
diversity_mols_to_seed_first_generation=0

# ============================================================================
# GENETIC ALGORITHM - SUBSEQUENT GENERATIONS PARAMETERS
# Settings for generations 2 and beyond
# ============================================================================
number_of_mutants=90
number_of_crossovers=0
number_elitism_advance_from_previous_gen=10
top_mols_to_seed_next_generation=10
diversity_seed_depreciation_per_gen=1
num_generations=15

# ============================================================================
# COMPUTATIONAL RESOURCES AND PERFORMANCE
# Settings for parallel processing and system resources
# ============================================================================
procs_per_node=-1
docking_exhaustiveness=1  # For debugging, set to 1
multithread_mode="serial"  # Because cuda
docking_nprocs=16
mutants_per_batch=-1
max_pass_per_batch=10

# ============================================================================
# FRAGMENT AND MUTATION OPERATIONS
# Settings for chemical modifications and constraints
# ============================================================================
max_fragment_mol_weight=100

# ============================================================================
# MOLECULAR FILTERING - INTERACTION CONSTRAINTS
# Define required molecular interactions with specific residues
# ============================================================================
hb_donor_interaction_filter="GLY863"
hb_acceptor_interaction_filter="GLY863"

# ============================================================================
# MOLECULAR FILTERING - STRUCTURAL CONSTRAINTS  
# Substructure patterns that molecules must contain
# ============================================================================
substructure_smiles="c1cc(c([nH]n1)=O)"

# ============================================================================
# DEEP LEARNING FILTERING
# AI-based molecular filtering using DeepFrag
# ============================================================================
deepfrag_cosine_similarity_cutoff=0.45
deepfrag_model="all_best"

# ============================================================================
# RUN AUTOGROW4
# ============================================================================
python ../run_autogrow.py \
    --output_directory="$output_directory" \
    --receptor_path="$receptor_path" \
    --center_x="$center_x" \
    --center_y="$center_y" \
    --center_z="$center_z" \
    --size_x="$size_x" \
    --size_y="$size_y" \
    --size_z="$size_z" \
    --source_compound_file="$source_compound_file" \
    --number_of_mutants_first_generation="$number_of_mutants_first_generation" \
    --number_of_crossovers_first_generation="$number_of_crossovers_first_generation" \
    --top_mols_to_seed_next_generation_first_generation="$top_mols_to_seed_next_generation_first_generation" \
    --number_elitism_advance_from_previous_gen_first_generation="$number_elitism_advance_from_previous_gen_first_generation" \
    --diversity_mols_to_seed_first_generation="$diversity_mols_to_seed_first_generation" \
    --number_of_mutants="$number_of_mutants" \
    --number_of_crossovers="$number_of_crossovers" \
    --number_elitism_advance_from_previous_gen="$number_elitism_advance_from_previous_gen" \
    --MergeMCS \
    --top_mols_to_seed_next_generation="$top_mols_to_seed_next_generation" \
    --diversity_seed_depreciation_per_gen="$diversity_seed_depreciation_per_gen" \
    --num_generations="$num_generations" \
    --procs_per_node="$procs_per_node" \
    --FragmentAddition \
    --fragment_addition_prevent_backtracking \
    --FragmentRemoval \
    --fragment_removal_prevent_backtracking \
    --rxn_library_path="$rxn_library_path" \
    --VinaLikeDocking \
    --RankSelector \
    --multithread_mode="$multithread_mode" \
    --ObabelSmiTo3DSDF \
    --PythonMultiprocessing \
    --slurm_template_file="$slurm_template_file" \
    --obabel_path="$obabel_path" \
    --docking_executable="$docking_executable" \
    --RDKitToolkit \
    --HBDonorInteractionFilter "$hb_donor_interaction_filter" \
    --HBAcceptorInteractionFilter "$hb_acceptor_interaction_filter" \
    --SubstructureFilter \
    --substructure_smiles "$substructure_smiles" \
    --DeepFragFilterRDKit \
    --deepfrag_cosine_similarity_cutoff "$deepfrag_cosine_similarity_cutoff" \
    --deepfrag_model "$deepfrag_model" \
    --docking_nprocs "$docking_nprocs" \
    --max_fragment_mol_weight "$max_fragment_mol_weight" \
    --mutants_per_batch "$mutants_per_batch" \
    --max_pass_per_batch "$max_pass_per_batch"
