#!/bin/bash
#SBATCH --job-name=PARP_Run_0
#SBATCH --output=/bgfs/jdurrant/jspiegel/test_mpi/PARP_Run/PARP_Run_0.conf.out
#SBATCH --time=00:45:00
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=12
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

highest_folder="/bgfs/jdurrant/jspiegel/test_mpi/PARP_Run/"
mkdir $highest_folder


average_time=0
for i in 1
do
    outfolder_four=$highest_folder"Run_$i/"
    mkdir $outfolder_four
    
    echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_QVINA2_3"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13

    srun -n $SLURM_NTASKS --mpi=pmi2 python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
        --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        --center_x -70.76 --center_y  21.82 --center_z 28.33 \
        --size_x 25.0 --size_y 16.0 --size_z 25.0 \
        --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/ZINC_fragments.smi \
        --root_output_folder $outfolder_four \
        --number_of_mutants_first_generation 500 \
        --number_of_crossovers_first_generation 500 \
        --number_of_mutants 500 \
        --number_of_crossovers 500 \
        --number_to_advance_from_previous_gen 500 \
        --top_mols_to_seed_next_generation 500 \
        --diversity_mols_to_seed_first_generation 500 \
        --diversity_seed_depreciation_per_gen 5 \
        --num_generations 40 \
        --mgltools_directory $MGLTOOLS_HOME/ \
        --number_of_processors -1 \
        --scoring_choice VINA \
        --Lipinski_Strict \
        --Ghose \
        --PAINS_Filter \
        --start_a_new_run \
        --filter_source_compounds False \
        --rxn_library Robust_Rxns \
        --selector_choice Rank_Selector \
        --dock_choice QuickVina2Docking \
        --max_variants_per_compound 3 \
        --multithread_mode mpi
        >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

done