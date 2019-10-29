#!/bin/bash
#SBATCH --job-name=PARPis
#SBATCH --output=/bgfs/jdurrant/jspiegel/DockPARPi/PARPis.conf.out
#SBATCH --time=70:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=20
#SBATCH --cluster=mpi
#SBATCH --partition=opa


## Define the environment
#py2.7 settings 
module purge
module load mgltools
module load gcc/8.2.0
module load python/anaconda2.7-2018.12_westpa

#py3.6 settings 
# module purge
# module load mgltools
# module load gcc/8.2.0
# module load python/anaconda3.7-2018.12_westpa


highest_folder="/bgfs/jdurrant/jspiegel/DockPARPi/"
mkdir $highest_folder

##Pre Run the process to prevent EOF Errors                                     
mpirun -n 1 --mpi=pmi2 python -m mpi4py RunAutogrow.py -c

average_time=0
for i in 1
do
    outfolder_four=$highest_folder"Run_$i/"
    mkdir $outfolder_four
    
    echo "START Autogrow 4.0 Run number $i"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13
    
    python -m mpi4py /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
        --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        --center_x -70.76 --center_y  21.82 --center_z 28.33 \
        --size_x 25.0 --size_y 16.0 --size_z 25.0 \
        --source_compound_file /ihome/jdurrant/jspiegel/to_dock/LGMNB_PAINS.smi \
        --root_output_folder $outfolder_four \
        --number_of_mutants_first_generation 0 \
        --number_of_crossovers_first_generation 0 \
        --number_of_mutants 0 \
        --number_of_crossovers 0 \
        --number_to_advance_from_previous_gen 6312 \
        --top_mols_to_seed_next_generation 6312 \
        --diversity_mols_to_seed_first_generation 0 \
        --diversity_seed_depreciation_per_gen 0 \
        --num_generations 1 \
        --mgltools_directory $MGLTOOLS_HOME/ \
        --number_of_processors -1 \
        --Dock_choice VinaDocking \
        --Scoring_choice NN1 \
        --Selector_Choice Rank_Selector \
        --docking_exhaustiveness 100 \
        --docking_num_modes 25 \
        --max_variants_per_compound 10 \
        --gypsum_thoroughness 5 \
        --start_a_new_run \
        --reduce_files_sizes True \
        --multithread_mode mpi \
        --generate_plot False \
        --redock_advance_from_previous_gen True \
        --filter_source_compounds False \
        >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

done    