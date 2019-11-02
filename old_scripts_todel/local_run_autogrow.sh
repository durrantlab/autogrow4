highest_folder = /home/jacob/Desktop/Outputfolder/trial/

##Pre Run the process to prevent EOF Errors                                     
python -m mpi4py RunAutogrow.py -c
    
average_time=0
for i in 1
do
    outfolder_four=$highest_folder"Run_$i/"
    mkdir $outfolder_four
    
    echo "START Autogrow 4.0 Run number $i"
    start_time="$(date +%s%N | cut -b1-13)"
    date +%s%N | cut -b1-13
    
    time python /home/jacob/Documents/autogrow4/RunAutogrow.py \
        --filename_of_receptor /home/jacob/Documents/autogrow4/tutorial/PARP/4r6e_removed_smallmol_aligned_Hs.pdb \
        --center_x -70.76 --center_y  21.82 --center_z 28.33 \
        --size_x 25.0 --size_y 16.0 --size_z 25.0 \
        --source_compound_file /home/jacob/Documents/autogrow4/source_compounds/Fragment_MW_100_to_150_docked.smi \
        --root_output_folder $outfolder_four \
        --number_of_mutants_first_generation 60 \
        --number_of_crossovers_first_generation 60 \
        --number_of_mutants 60 \
        --number_of_crossovers 60 \
        --number_to_advance_from_previous_gen 60 \
        --top_mols_to_seed_next_generation 30 \
        --diversity_mols_to_seed_first_generation 30 \
        --diversity_seed_depreciation_per_gen 0 \
        --num_generations 5 \
        --mgltools_directory /home/jacob/MGLTools-1.5.6/ \
        --number_of_processors -1 \
        --dock_choice QuickVina2Docking \
        --scoring_choice VINA \
        --Lipinski_Lenient \
        --start_a_new_run \
        --rxn_library ClickChem \
        --selector_choice Rank_Selector \
        --max_variants_per_compound 3 \
        --filter_source_compounds False \
        --redock_advance_from_previous_gen False \
        --generate_plot 'False'
        
done    
