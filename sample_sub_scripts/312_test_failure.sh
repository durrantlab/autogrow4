#!/bin/bash
#SBATCH --job-name=312_failrate
#SBATCH --output=/bgfs/jdurrant/jspiegel/output_crc/312_failrate/312_failrate.conf.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
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

highest_folder="/bgfs/jdurrant/jspiegel/output_crc/312_failrate/"


higher_outdir=$highest_folder"Run_0/"
mkdir $higher_outdir


average_time=0
for i in 1 2 3 4 5 6 7 8 9 10
do
    outfolder=$higher_outdir"Run_$i/"
    mkdir $outfolder
    output_txt=$outfolder"test_output.txt"
    error_txt=$outfolder"test_error.txt"
    echo "START Autogrow 3.1.2 Run number $i"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13

    time timeout 6000 python /bgfs/jdurrant/jspiegel/test_old_autogrow/autogrow/3.1.2/autogrow/autogrow_3_1_2.py \
        --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        -center_x -70.76 -center_y  21.82 -center_z 28.33 \
        -size_x 25.0 -size_y 16.0 -size_z 25.0 \
        -additional_autoclickchem_parameters "-all_reactions +azide_and_alkyne_to_azole" \
        -allow_modification_without_frag_addition FALSE \
        -directory_of_source_compounds /bgfs/jdurrant/jspiegel/test_old_autogrow/autogrow/3.1.2/autogrow/tutorial/starting_compounds/ \
        -directory_of_fragments /bgfs/jdurrant/jspiegel/test_old_autogrow/autogrow/3.1.2/autogrow/fragments/MW_250/ \
        -number_of_mutants_first_generation 100 \
        -number_of_crossovers_first_generation 100 \
        -number_of_mutants 100 -number_of_crossovers 100 \
        -top_ones_to_advance_to_next_generation 20 \
        -num_generations 5 -max_seconds_per_generation 6000 \
        -use_lipinski_filter FALSE -use_strict_lipinski_filter TRUE -use_ghose_filter FALSE \
        -scoring_function NN1 -score_by_ligand_efficiency FALSE -maintain_core FALSE \
        -minimum_core_atoms_required 4 \
        -vina_executable /bgfs/jdurrant/jspiegel/autogrow/autogrow/Docking/Docking_Executables/Vina/autodock_vina_1_1_2_linux_x86/bin/vina \
        -openbabel_bin_directory /ihome/crc/install/python/anaconda2.7-4.2.0/bin/ \
        -mgltools_directory $MGLTOOLS_HOME/ \
        -num_processors -1 -output_dir $outfolder \
        >  $outfolder"test_output.txt" 2>  $outfolder"test_error.txt"
        exit_status=$?
    echo ""
    echo "##############################"
    echo "EXIT STATUS!!!!!!: "$exit_status
    echo "##############################"
    echo ""
    end_time="$(date +%s%N | cut -b1-13)"
    echo $end_time "   end"
    tot_time=$(($end_time-$start_time))
    echo $tot_time " Total time"
    tot_time_all=$(($tot_time_all + $tot_time))
    echo $tot_time_all "Total time ALL"

done