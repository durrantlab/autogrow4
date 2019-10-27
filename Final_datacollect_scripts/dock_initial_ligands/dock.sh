#!/bin/bash
#SBATCH --job-name=exhautive_dock
#SBATCH --output=/bgfs/jdurrant/jspiegel/docked_source/exhautive_dock.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=23:59:00                                                                       
#SBATCH --cluster=invest                                                                            
#SBATCH --partition=jdurrant                                                                        

#SBATCH --mail-type BEGIN,END,FAIL,ARRAY_TASKS  
#SBATCH --mail-user jako134j@gmail.com       
## Define the environment
#py2.7 settings 
module purge
module load mgltools
module load gcc/8.2.0
# module load python/anaconda2.7-2018.12_westpa

#py3.6 settings 
# module purge
# module load mgltools
# module load gcc/8.2.0
# module load python/anaconda3.7-2018.12_westpa

highest_folder="/bgfs/jdurrant/jspiegel/docked_source/"
mkdir $highest_folder

average_time=0
for i in 1
do
    outfolder_four=$highest_folder"Run_$i/"
    mkdir $outfolder_four
    
    echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_QVINA2_3"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13

    ~/miniconda3/envs/py37/bin/python /bgfs/jdurrant/jspiegel/autogrow/RunAutogrow.py \
        --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        --center_x -70.76 --center_y  21.82 --center_z 28.33 \
        --size_x 25.0 --size_y 16.0 --size_z 25.0 \
        --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/ZINC_fragments.smi \
        --root_output_folder $outfolder_four \
        --number_of_mutants_first_generation 0 \
        --number_of_crossovers_first_generation 0 \
        --number_to_advance_from_previous_gen_first_generation 930 \
        --number_of_mutants 500 \
        --number_of_crossovers 500 \
        --number_to_advance_from_previous_gen 5000 \
        --top_mols_to_seed_next_generation_first_generation 930 \
        --top_mols_to_seed_next_generation 500 \
        --diversity_mols_to_seed_first_generation 500 \
        --diversity_seed_depreciation_per_gen 5 \
        --num_generations 1 \
        --mgltools_directory $MGLTOOLS_HOME/ \
        --number_of_processors -1 \
        --Dock_choice QuickVina2Docking \
        --Scoring_choice VINA \
        --Selector_Choice Rank_Selector \
        --No_Filters \
        --reduce_files_sizes False \
        --max_variants_per_compound 10 \
        --redock_advance_from_previous_gen False \
        --filter_source_compounds False \
        --use_docked_source_compounds True \
        --Rxn_library All_Rxns \
        --docking_exhaustiveness 50 \
        --generate_plot True \
        --multithread_mode multithreading \
        >>  $outfolder_four"test_output_$i.txt" 2>>  $outfolder_four"test_error_$i.txt"

done