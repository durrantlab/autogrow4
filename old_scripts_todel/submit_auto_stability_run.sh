#!/bin/bash
#SBATCH --job-name=stability_run_PARP
#SBATCH --output=stability_run_PARP.conf.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=48:05:00
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

highest_folder="/bgfs/jdurrant/jspiegel/output_crc/BenchMarking/Run_0/"
mkdir $highest_folder




higher_outdir=$highest_folder"312/"
mkdir $higher_outdir

average_time=0
for i in 1
do
    outfolder_four=$higher_outdir"Run_$i/"
    mkdir $outfolder_four
    output_txt=$outfolder"test_output.txt"
    error_txt=$outfolder"test_error.txt"
    echo "START Autogrow 3.1.2 Run number $i"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13

    time python /bgfs/jdurrant/jspiegel/old_autogrow/autogrow/autogrow/autogrow_3_1_2.py \
            -filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/4.0.0/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
            -center_x -70.76 -center_y  21.82 -center_z 28.33 \
            -size_x 25.0 -size_y 16.0 -size_z 25.0 \
            -additional_autoclickchem_parameters "-all_reactions +azide_and_alkyne_to_azole" \
            -allow_modification_without_frag_addition FALSE \
            -directory_of_source_compounds /bgfs/jdurrant/jspiegel/old_autogrow/autogrow/autogrow/tutorial/starting_compounds/ \
            -directory_of_fragments /bgfs/jdurrant/jspiegel/old_autogrow/autogrow/autogrow/fragments/MW_250/ \
            -number_of_mutants_first_generation 100 \
            -number_of_crossovers_first_generation 100 \
            -number_of_mutants 100 -number_of_crossovers 100 \
            -top_ones_to_advance_to_next_generation 20 \
            -num_generations 5 -max_seconds_per_generation 6000 \
            -use_lipinski_filter FALSE -use_strict_lipinski_filter TRUE -use_ghose_filter FALSE \
            -scoring_function VINA -score_by_ligand_efficiency FALSE -maintain_core FALSE \
            -minimum_core_atoms_required 4 \
            -vina_executable /bgfs/jdurrant/jspiegel/autogrow4/autogrow/Docking/Docking_Executables/Vina/autodock_vina_1_1_2_linux_x86/bin/vina \
            -openbabel_bin_directory /ihome/crc/install/python/anaconda2.7-4.2.0/bin/ \
            -mgltools_directory $MGLTOOLS_HOME/ \
            -num_processors -1 -output_dir $outfolder \
            >  $outfolder"test_output.txt" 2>  $outfolder"test_error.txt"

    end_time="$(date +%s%N | cut -b1-13)"
    echo $end_time "   end"
    tot_time=$(($end_time-$start_time))
    echo $tot_time " Total time"
    tot_time_all=$(($tot_time_all + $tot_time))
    echo $tot_time_all "Total time ALL"

done




echo ""
echo "" 
echo ""
date +%s%N | cut -b1-13
echo "Start Autogrow 4.0 STABILITY RUN USING with Rank_VINA_1"
date +%s%N | cut -b1-13



higher_outdir=$highest_folder"4/"
mkdir $higher_outdir
runtypedir=$higher_outdir"Rank_VINA_1/"
mkdir $runtypedir

average_time=0
for i in 1 2 3 4 5
do
    outfolder_four=$runtypedir"Run_$i/"
    mkdir $outfolder_four
    output_txt=$outfolder_four"test_output.txt"
    error_txt=$outfolder_four"test_error.txt"
    echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_VINA_1"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13

    time python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
        --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        --center_x -70.76 --center_y  21.82 --center_z 28.33 \
        --size_x 25.0 --size_y 16.0 --size_z 25.0 \
        --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/Fragment_MW_100_to_150_docked.smi \
        --root_output_folder $outfolder_four \
        --number_of_mutants_first_generation 100 \
        --number_of_crossovers_first_generation 100 \
        --number_of_mutants 100 \
        --number_of_crossovers 100 \
        --number_elitism_advance_from_previous_gen 0 \
        --top_mols_to_seed_next_generation 100 \
        --diversity_mols_to_seed_first_generation 100 \
        --diversity_seed_depreciation_per_gen 0 \
        --num_generations 5 \
        --mgltools_directory $MGLTOOLS_HOME/ \
        --number_of_processors -1 \
        --scoring_choice VINA \
        --Lipinski_Strict \
        --start_a_new_run \
        --filter_source_compounds False \
        --rxn_library ClickChem \
        --selector_choice Rank_Selector \
        --dock_choice VinaDocking \
        --max_variants_per_compound 1 \
        >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

    end_time="$(date +%s%N | cut -b1-13)"
    echo $end_time "   end"
    tot_time=$(($end_time-$start_time))
    echo $tot_time " Total time"
    tot_time_all=$(($tot_time_all + $tot_time))
    echo $tot_time_all "Total time ALL"
done
echo "finished Autogrow 4.0 STABILITY RUN USING with Rank_VINA_1"
echo "Total time for $i runs Autogrow 4.0 STABILITY RUN USING with Rank_VINA_1"
echo $tot_time_all 
average_time=$(($tot_time_all/$i))
echo $average_time "AVERAGE TIME FOR $i RUNS Autogrow 4.0 STABILITY RUN USING with Rank_VINA_1"
date +%s%N | cut -b1-13
echo "#######################"
