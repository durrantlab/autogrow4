#!/bin/bash
#SBATCH --job-name=exhautive_dock_MW_150
#SBATCH --output=/bgfs/jdurrant/jspiegel/docked_source/exhautive_dock_MW_150.txt
#SBATCH --time=23:55:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=28
#SBATCH --cluster=mpi
#SBATCH --partition=opa
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
source activate py37

# run precache
~/miniconda3/envs/py37/bin/python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py -c

highest_folder="/bgfs/jdurrant/jspiegel/docked_source/"
# mkdir $highest_folder

average_time=0
for i in 1
do
    outfolder_four=$highest_folder"MW_150_to_200/"
    mkdir $outfolder_four
    
    echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_QVINA2_3"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13

    mpirun -n $SLURM_NTASKS \
    ~/miniconda3/envs/py37/bin/python -m mpi4py /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
        --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        --center_x -70.76 --center_y  21.82 --center_z 28.33 \
        --size_x 25.0 --size_y 16.0 --size_z 25.0 \
        --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/Fragment_MW_150_to_200.smi \
        --root_output_folder $outfolder_four \
        --num_generations 0 \
        --mgltools_directory $MGLTOOLS_HOME/ \
        --number_of_processors -1 \
        --dock_choice QuickVina2Docking \
        --scoring_choice VINA \
        --selector_choice Rank_Selector \
        --No_Filters \
        --reduce_files_sizes True \
        --max_variants_per_compound 25 \
        --redock_advance_from_previous_gen False \
        --filter_source_compounds False \
        --use_docked_source_compounds True \
        --docking_exhaustiveness 50 \
        --multithread_mode multithreading \
        >>  $outfolder_four"test_output_$i.txt" 2>>  $outfolder_four"test_error_$i.txt"

done