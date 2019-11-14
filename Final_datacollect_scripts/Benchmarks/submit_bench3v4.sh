#!/bin/bash
#SBATCH --job-name=bench_3v4_smallSMP_Run_0
#SBATCH --output=/bgfs/jdurrant/jspiegel/Benchmarks/bench_3v4_Run_0.conf.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=23:59:00
#SBATCH --cluster=smp
#SBATCH --partition=smp

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

highest_folder="/bgfs/jdurrant/jspiegel/Benchmarks/Run_0/"
mkdir $highest_folder


##Pre Run the process to prevent EOF Errors                                     
python RunAutogrow.py -c

higher_outdir=$highest_folder"313/"
mkdir $higher_outdir

# 3.1.3 needs the top_ones_to_advance_to_next_generation to seed
# but for every 1 from advanced then you create and dock 1 less.
# so 3.1.3 needs to be have a higher mutant/crossover number
#     so for 3.1.3 there will be 150mut and 150cross with an advance of 100
#       This will result in 200 docking events per generation
#       for 4.0.0 the equivalent is 100mut and 100cross with an advance of 0
#           -this results in 200 docking events per generations.
# Also generation_num needs to be 1 higher since the last gen is a copy over
#
# with timeout of 100seconds ~80% of all dockings complete for Vina Docking
#       Default setting (which is used in this benchmark) is 200 seconds
#   
average_time=0
for i in 1
do
    outfolder=$higher_outdir"Run_$i/"
    mkdir $outfolder
    output_txt=$outfolder"test_output.txt"
    error_txt=$outfolder"test_error.txt"
    echo "START Autogrow 3.1.3 Run number $i"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13

    time timeout 86400 python /bgfs/jdurrant/jspiegel/old_autogrow/autogrow/autogrow/autogrow_3_1_3.py \
        --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        -center_x -70.76 -center_y  21.82 -center_z 28.33 \
        -size_x 25.0 -size_y 16.0 -size_z 25.0 \
        -additional_autoclickchem_parameters "-all_reactions +azide_and_alkyne_to_azole" \
        -allow_modification_without_frag_addition FALSE \
        -directory_of_source_compounds /bgfs/jdurrant/jspiegel/old_autogrow/autogrow/autogrow/tutorial/starting_compounds/ \
        -directory_of_fragments /bgfs/jdurrant/jspiegel/old_autogrow/autogrow/autogrow/fragments/MW_150/ \
        -number_of_mutants_first_generation 50 \
        -number_of_crossovers_first_generation 50 \
        -number_of_mutants 85 \
        -number_of_crossovers 85 \
        -top_ones_to_advance_to_next_generation 70 \
        -num_generations 6 -max_seconds_per_generation 21600 \
        -use_lipinski_filter TRUE -use_strict_lipinski_filter TRUE -use_ghose_filter TRUE \
        -scoring_function VINA -score_by_ligand_efficiency FALSE -maintain_core FALSE \
        -minimum_core_atoms_required 4 \
        -vina_executable /bgfs/jdurrant/jspiegel/autogrow4/autogrow/docking/docking_executables/vina/autodock_vina_1_1_2_linux_x86/bin/vina \
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


if [ $exit_status = "124" ]
then
    echo "Fail to Finish AUTOGROW 3.1.3"
    echo "TERMINATE this Run!!!!!!1"
fi

#CHECK IT FINISHED GEN 6 in 3.1.3 run
if [ -d $higher_outdir"/Run_1/generation6/" ]
then
    echo "AUTOGROW 3.1.3 completed"
else
    echo "AUTOGROW 3.1.3 failed"
    echo "TERMINATE this Run!!!!!!"
    exit 0
fi

if  [ $exit_status = "0" ]
then

    # RANKING VINA MaxVar 1  

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
    for i in 1 
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
            --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/naphthalene_smiles.smi \
            --root_output_folder $outfolder_four \
            --number_of_mutants_first_generation 50 \
            --number_of_crossovers_first_generation 50 \
            --number_of_mutants 50 \
            --number_of_crossovers 50 \
            --top_mols_to_seed_next_generation 70 \
            --number_elitism_advance_from_previous_gen 70 \
            --number_elitism_advance_from_previous_gen_first_generation 0 \
            --diversity_mols_to_seed_first_generation 0 \
            --diversity_seed_depreciation_per_gen 0 \
            --num_generations 5 \
            --mgltools_directory $MGLTOOLS_HOME/ \
            --number_of_processors -1 \
            --scoring_choice VINA \
            --LipinskiStrictFilter \
            --GhoseFilter \
            --start_a_new_run \
            --rxn_library click_chem_rxns \
            --selector_choice Rank_Selector \
            --dock_choice VinaDocking \
            --max_variants_per_compound 1 \
            --redock_elite_from_previous_gen False \
            --generate_plot False \
            --debug_mode \
            --reduce_files_sizes False \
            --use_docked_source_compounds False \
            --gypsum_timeout_limit 60 \
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


    # RANKING VINA MaxVar 3

    echo ""
    echo "" 
    echo ""
    date +%s%N | cut -b1-13
    echo "Start Autogrow 4.0 STABILITY RUN USING with Rank_VINA_3"
    date +%s%N | cut -b1-13



    higher_outdir=$highest_folder"4/"
    mkdir $higher_outdir
    runtypedir=$higher_outdir"Rank_VINA_3/"
    mkdir $runtypedir

    average_time=0
    for i in 1 
    do
        outfolder_four=$runtypedir"Run_$i/"
        mkdir $outfolder_four
        output_txt=$outfolder_four"test_output.txt"
        error_txt=$outfolder_four"test_error.txt"
        echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_VINA_3"
        start_time="$(date +%s%N | cut -b1-13)"
        date +%s%N | cut -b1-13

        time python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
            --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
            --center_x -70.76 --center_y  21.82 --center_z 28.33 \
            --size_x 25.0 --size_y 16.0 --size_z 25.0 \
            --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/naphthalene_smiles.smi \
            --root_output_folder $outfolder_four \
            --number_of_mutants_first_generation 50 \
            --number_of_crossovers_first_generation 50 \
            --number_of_mutants 50 \
            --number_of_crossovers 50 \
            --top_mols_to_seed_next_generation 70 \
            --number_elitism_advance_from_previous_gen 70 \
            --number_elitism_advance_from_previous_gen_first_generation 0 \
            --diversity_mols_to_seed_first_generation 0 \
            --diversity_seed_depreciation_per_gen 0 \
            --num_generations 5 \
            --mgltools_directory $MGLTOOLS_HOME/ \
            --number_of_processors -1 \
            --scoring_choice VINA \
            --LipinskiStrictFilter \
            --GhoseFilter \
            --start_a_new_run \
            --rxn_library click_chem_rxns \
            --selector_choice Rank_Selector \
            --dock_choice VinaDocking \
            --max_variants_per_compound 3 \
            --redock_elite_from_previous_gen False \
            --generate_plot False \
            --debug_mode \
            --reduce_files_sizes False \
            --use_docked_source_compounds False \
            --gypsum_timeout_limit 60 \
            >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

        end_time="$(date +%s%N | cut -b1-13)"
        echo $end_time "   end"
        tot_time=$(($end_time-$start_time))
        echo $tot_time " Total time"
        tot_time_all=$(($tot_time_all + $tot_time))
        echo $tot_time_all "Total time ALL"
    done
    echo "finished Autogrow 4.0 STABILITY RUN USING with Rank_VINA_3"
    echo "Total time for $i runs Autogrow 4.0 STABILITY RUN USING with Rank_VINA_3"
    echo $tot_time_all 
    average_time=$(($tot_time_all/$i))
    echo $average_time "AVERAGE TIME FOR $i RUNS Autogrow 4.0 STABILITY RUN USING with Rank_VINA_3"
    date +%s%N | cut -b1-13
    echo "#######################"





    # RANKING VINA MaxVar 5 

    echo ""
    echo "" 
    echo ""
    date +%s%N | cut -b1-13
    echo "Start Autogrow 4.0 STABILITY RUN USING with Rank_VINA_5"
    date +%s%N | cut -b1-13



    higher_outdir=$highest_folder"4/"
    mkdir $higher_outdir
    runtypedir=$higher_outdir"Rank_VINA_5/"
    mkdir $runtypedir

    average_time=0
    for i in 1 
    do
        outfolder_four=$runtypedir"Run_$i/"
        mkdir $outfolder_four
        output_txt=$outfolder_four"test_output.txt"
        error_txt=$outfolder_four"test_error.txt"
        echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_VINA_5"
        start_time="$(date +%s%N | cut -b1-13)"
        date +%s%N | cut -b1-13

        time python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
            --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
            --center_x -70.76 --center_y  21.82 --center_z 28.33 \
            --size_x 25.0 --size_y 16.0 --size_z 25.0 \
            --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/naphthalene_smiles.smi \
            --root_output_folder $outfolder_four \
            --number_of_mutants_first_generation 50 \
            --number_of_crossovers_first_generation 50 \
            --number_of_mutants 50 \
            --number_of_crossovers 50 \
            --top_mols_to_seed_next_generation 70 \
            --number_elitism_advance_from_previous_gen 70 \
            --number_elitism_advance_from_previous_gen_first_generation 0 \
            --diversity_mols_to_seed_first_generation 0 \
            --diversity_seed_depreciation_per_gen 0 \
            --num_generations 5 \
            --mgltools_directory $MGLTOOLS_HOME/ \
            --number_of_processors -1 \
            --scoring_choice VINA \
            --LipinskiStrictFilter \
            --GhoseFilter \
            --start_a_new_run \
            --rxn_library click_chem_rxns \
            --selector_choice Rank_Selector \
            --dock_choice VinaDocking \
            --max_variants_per_compound 5 \
            --redock_elite_from_previous_gen False \
            --generate_plot False \
            --debug_mode \
            --reduce_files_sizes False \
            --use_docked_source_compounds False \
            --gypsum_timeout_limit 60 \
            >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

        end_time="$(date +%s%N | cut -b1-13)"
        echo $end_time "   end"
        tot_time=$(($end_time-$start_time))
        echo $tot_time " Total time"
        tot_time_all=$(($tot_time_all + $tot_time))
        echo $tot_time_all "Total time ALL"
    done
    echo "finished Autogrow 4.0 STABILITY RUN USING with Rank_VINA_5"
    echo "Total time for $i runs Autogrow 4.0 STABILITY RUN USING with Rank_VINA_5"
    echo $tot_time_all 
    average_time=$(($tot_time_all/$i))
    echo $average_time "AVERAGE TIME FOR $i RUNS Autogrow 4.0 STABILITY RUN USING with Rank_VINA_5"
    date +%s%N | cut -b1-13
    echo "#######################"



    # USING QUICKVINA



    # RANKING QVINA2 MaxVar 1  

    echo ""
    echo "" 
    echo ""
    date +%s%N | cut -b1-13
    echo "Start Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_1"
    date +%s%N | cut -b1-13



    higher_outdir=$highest_folder"4/"
    mkdir $higher_outdir
    runtypedir=$higher_outdir"Rank_QVINA2_1/"
    mkdir $runtypedir

    average_time=0
    for i in 1 
    do
        outfolder_four=$runtypedir"Run_$i/"
        mkdir $outfolder_four
        output_txt=$outfolder_four"test_output.txt"
        error_txt=$outfolder_four"test_error.txt"
        echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_QVINA2_1"
        start_time="$(date +%s%N | cut -b1-13)"
        date +%s%N | cut -b1-13

        time python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
            --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
            --center_x -70.76 --center_y  21.82 --center_z 28.33 \
            --size_x 25.0 --size_y 16.0 --size_z 25.0 \
            --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/naphthalene_smiles.smi \
            --root_output_folder $outfolder_four \
            --number_of_mutants_first_generation 50 \
            --number_of_crossovers_first_generation 50 \
            --number_of_mutants 50 \
            --number_of_crossovers 50 \
            --top_mols_to_seed_next_generation 70 \
            --number_elitism_advance_from_previous_gen 70 \
            --number_elitism_advance_from_previous_gen_first_generation 0 \
            --diversity_mols_to_seed_first_generation 0 \
            --diversity_seed_depreciation_per_gen 0 \
            --num_generations 5 \
            --mgltools_directory $MGLTOOLS_HOME/ \
            --number_of_processors -1 \
            --scoring_choice VINA \
            --LipinskiStrictFilter \
            --GhoseFilter \
            --start_a_new_run \
            --rxn_library click_chem_rxns \
            --selector_choice Rank_Selector \
            --dock_choice QuickVina2Docking \
            --max_variants_per_compound 1 \
            --redock_elite_from_previous_gen False \
            --generate_plot False \
            --debug_mode \
            --reduce_files_sizes False \
            --use_docked_source_compounds False \
            --gypsum_timeout_limit 60 \
            >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

        end_time="$(date +%s%N | cut -b1-13)"
        echo $end_time "   end"
        tot_time=$(($end_time-$start_time))
        echo $tot_time " Total time"
        tot_time_all=$(($tot_time_all + $tot_time))
        echo $tot_time_all "Total time ALL"
    done
    echo "finished Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_1"
    echo "Total time for $i runs Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_1"
    echo $tot_time_all 
    average_time=$(($tot_time_all/$i))
    echo $average_time "AVERAGE TIME FOR $i RUNS Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_1"
    date +%s%N | cut -b1-13
    echo "#######################"


    # RANKING QVINA2 MaxVar 3

    echo ""
    echo "" 
    echo ""
    date +%s%N | cut -b1-13
    echo "Start Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_3"
    date +%s%N | cut -b1-13



    higher_outdir=$highest_folder"4/"
    mkdir $higher_outdir
    runtypedir=$higher_outdir"Rank_QVINA2_3/"
    mkdir $runtypedir

    average_time=0
    for i in 1 
    do
        outfolder_four=$runtypedir"Run_$i/"
        mkdir $outfolder_four
        output_txt=$outfolder_four"test_output.txt"
        error_txt=$outfolder_four"test_error.txt"
        echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_QVINA2_3"
        start_time="$(date +%s%N | cut -b1-13)"
        date +%s%N | cut -b1-13

        time python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
            --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
            --center_x -70.76 --center_y  21.82 --center_z 28.33 \
            --size_x 25.0 --size_y 16.0 --size_z 25.0 \
            --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/naphthalene_smiles.smi \
            --root_output_folder $outfolder_four \
            --number_of_mutants_first_generation 50 \
            --number_of_crossovers_first_generation 50 \
            --number_of_mutants 50 \
            --number_of_crossovers 50 \
            --top_mols_to_seed_next_generation 70 \
            --number_elitism_advance_from_previous_gen 70 \
            --number_elitism_advance_from_previous_gen_first_generation 0 \
            --diversity_mols_to_seed_first_generation 0 \
            --diversity_seed_depreciation_per_gen 0 \
            --num_generations 5 \
            --mgltools_directory $MGLTOOLS_HOME/ \
            --number_of_processors -1 \
            --scoring_choice VINA \
            --LipinskiStrictFilter \
            --GhoseFilter \
            --start_a_new_run \
            --rxn_library click_chem_rxns \
            --selector_choice Rank_Selector \
            --dock_choice QuickVina2Docking \
            --max_variants_per_compound 3 \
            --redock_elite_from_previous_gen False \
            --generate_plot False \
            --debug_mode \
            --reduce_files_sizes False \
            --use_docked_source_compounds False \
            --gypsum_timeout_limit 60 \
            >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

        end_time="$(date +%s%N | cut -b1-13)"
        echo $end_time "   end"
        tot_time=$(($end_time-$start_time))
        echo $tot_time " Total time"
        tot_time_all=$(($tot_time_all + $tot_time))
        echo $tot_time_all "Total time ALL"
    done
    echo "finished Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_3"
    echo "Total time for $i runs Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_3"
    echo $tot_time_all 
    average_time=$(($tot_time_all/$i))
    echo $average_time "AVERAGE TIME FOR $i RUNS Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_3"
    date +%s%N | cut -b1-13
    echo "#######################"





    # RANKING QVINA2 MaxVar 5 

    echo ""
    echo "" 
    echo ""
    date +%s%N | cut -b1-13
    echo "Start Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_5"
    date +%s%N | cut -b1-13



    higher_outdir=$highest_folder"4/"
    mkdir $higher_outdir
    runtypedir=$higher_outdir"Rank_QVINA2_5/"
    mkdir $runtypedir

    average_time=0
    for i in 1 
    do
        outfolder_four=$runtypedir"Run_$i/"
        mkdir $outfolder_four
        output_txt=$outfolder_four"test_output.txt"
        error_txt=$outfolder_four"test_error.txt"
        echo "START Autogrow 4.0 Run number $i  STABILITY RUN USING with Rank_QVINA2_5"
        start_time="$(date +%s%N | cut -b1-13)"
        date +%s%N | cut -b1-13

        time python /bgfs/jdurrant/jspiegel/autogrow4/RunAutogrow.py \
            --filename_of_receptor /bgfs/jdurrant/jspiegel/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
            --center_x -70.76 --center_y  21.82 --center_z 28.33 \
            --size_x 25.0 --size_y 16.0 --size_z 25.0 \
            --source_compound_file /bgfs/jdurrant/jspiegel/autogrow4/source_compounds/naphthalene_smiles.smi \
            --root_output_folder $outfolder_four \
            --number_of_mutants_first_generation 50 \
            --number_of_crossovers_first_generation 50 \
            --number_of_mutants 50 \
            --number_of_crossovers 50 \
            --top_mols_to_seed_next_generation 70 \
            --number_elitism_advance_from_previous_gen 70 \
            --number_elitism_advance_from_previous_gen_first_generation 0 \
            --diversity_mols_to_seed_first_generation 0 \
            --diversity_seed_depreciation_per_gen 0 \
            --num_generations 5 \
            --mgltools_directory $MGLTOOLS_HOME/ \
            --number_of_processors -1 \
            --scoring_choice VINA \
            --LipinskiStrictFilter \
            --GhoseFilter \
            --start_a_new_run \
            --rxn_library click_chem_rxns \
            --selector_choice Rank_Selector \
            --dock_choice QuickVina2Docking \
            --max_variants_per_compound 5 \
            --redock_elite_from_previous_gen False \
            --generate_plot False \
            --debug_mode \
            --reduce_files_sizes False \
            --use_docked_source_compounds False \
            --gypsum_timeout_limit 60 \
            >  $outfolder_four"test_output.txt" 2>  $outfolder_four"test_error.txt"

        end_time="$(date +%s%N | cut -b1-13)"
        echo $end_time "   end"
        tot_time=$(($end_time-$start_time))
        echo $tot_time " Total time"
        tot_time_all=$(($tot_time_all + $tot_time))
        echo $tot_time_all "Total time ALL"
    done
    echo "finished Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_5"
    echo "Total time for $i runs Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_5"
    echo $tot_time_all 
    average_time=$(($tot_time_all/$i))
    echo $average_time "AVERAGE TIME FOR $i RUNS Autogrow 4.0 STABILITY RUN USING with Rank_QVINA2_5"
    date +%s%N | cut -b1-13
    echo "#######################"

fi