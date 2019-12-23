#!/bin/bash
#SBATCH --job-name=SMP_PARP_Run_0
#SBATCH --output=/bgfs/jdurrant/jspiegel/test_autogrow/SMP_PARP_Run_0.conf.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=10:00:00
#SBATCH --cluster=smp
#SBATCH --partition=high-mem

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

highest_folder="/bgfs/jdurrant/jspiegel/test_autogrow/"
mkdir $highest_folder


average_time=0
for i in 1
do
    outfolder_four=$highest_folder"Run_$i/"
    mkdir $outfolder_four
    
    echo "START Autogrow 4.0 Run number $i"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13
    
    python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py -j /bgfs/jdurrant/jspiegel/autogrow4/sample_sub_scripts/PARP_SMP.json \
        >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

done    