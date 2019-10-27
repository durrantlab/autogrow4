nohup time python RunAutogrow.py \
    --filename_of_receptor /home/jspiegel/projects/autogrow/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
    --center_x -70.76 --center_y  21.82 --center_z 28.33 \
    --size_x 25.0 --size_y 16.0 --size_z 25.0 \
    --source_compound_file /home/jspiegel/projects/parpi_test/parpi_frags.smi \
    --root_output_folder /home/jspiegel/DataB/jspiegel/projects/Autogrow_output/ \
    --number_of_mutants_first_generation 50 \
    --number_of_crossovers_first_generation 50 \
    --number_of_mutants 50 \
    --number_of_crossovers 50 \
    --number_to_advance_from_previous_gen 50 \
    --top_mols_to_seed_next_generation 50 \
    --diversity_mols_to_seed_first_generation 50 \
    --diversity_seed_depreciation_per_gen 10 \
    --num_generations 5 \
    --mgltools_directory /home/jspiegel/DataB/spinel/programs/MGLTools-1.5.6/ \
    --number_of_processors 40 \
    --Dock_choice QuickVina2Docking \
    --Scoring_choice VINA \
    --Lipinski_Strict \
    --Ghose \
    --filter_source_compounds True \
    --Rxn_library Robust_Rxns \
    --Selector_Choice Rank_Selector \
    --max_variants_per_compound 3 \
    --reduce_files_sizes True \
    --generate_plot True \
    --use_docked_source_compounds True \
    >  /home/jspiegel/projects/Autogrow_output/test_output.txt 2>  /home/jspiegel/projects/Autogrow_output/test_error.txt
    