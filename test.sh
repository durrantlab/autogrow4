rm -f ./tutorial/PARP/4r6eA_PARP1_prepared.pdbqt
#rm -rf output_directory/*

python run_autogrow.py \
    --receptor_path ./tutorial/PARP/4r6eA_PARP1_prepared.pdb \
    --center_x -70.76 --center_y  21.82 --center_z 28.33 \
    --size_x 25.0 --size_y 16.0 --size_z 25.0 \
    --obabel_path /Users/jdurrant/opt/anaconda3/bin/obabel \
    --source_compound_file tmp_20.smi \
    --root_output_folder ./output_directory/ \
    --vina_like_executable /Applications/vina_1.2.5_mac_x86_64 \
    --number_of_mutants_first_generation 5 \
    --number_of_crossovers_first_generation 5 \
    --number_elitism_advance_from_previous_gen_first_generation 5 \
    --diversity_mols_to_seed_first_generation 5 \
    --number_of_mutants 5 \
    --number_of_crossovers 5 \
    --number_elitism_advance_from_previous_gen 5 \
    --MergeMCS \
    --top_mols_to_seed_next_generation 5 \
    --diversity_seed_depreciation_per_gen 1 \
    --num_generations 5 \
    --number_of_processors -1 \
    --FragmentAddition \
    --rxn_library_path ./autogrow/plugins/mutation/reaction_libraries/all_rxns \
    --VinaLikeDocking \
    --docking_exhaustiveness 1 \
    --redock_elite_from_previous_gen False \
    --generate_plot True \
    --RankSelector \
    --multithread_mode serial \
    --ObabelSmiTo3DSDF \
    --ParallelExecPlugin \
    --parallel_exec_path /opt/homebrew/bin/parallel

    # --PythonMultiprocessing
    
    #\
    #>  ./output_directory/text_file.txt 2> ./output_directory/text_errormessage_file.txt

