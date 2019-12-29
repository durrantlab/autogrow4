#!/bin/bash
python ./autogrow_in_docker.py
# python ./autogrow_in_docker.py -filename_of_receptor ./example_files/1xdn_receptor.pdb \
#      -center_x 39.5 -center_y 23.7 -center_z 15.0 \
#      -size_x 30.0 -size_y 30.0 -size_z 20.0 \
#      -additional_autoclickchem_parameters "-all_reactions +azide_and_alkyne_to_azole" \
#      -allow_modification_without_frag_addition FALSE \
#      -directory_of_source_compounds ./example_files/starting_compounds/ \
#      -directory_of_fragments MW_250 \
#      -number_of_mutants_first_generation 3 -number_of_crossovers_first_generation 3 \
#      -number_of_mutants 3 -number_of_crossovers 3 \
#      -top_ones_to_advance_to_next_generation 3 -num_generations 1 \
#      -max_seconds_per_generation 1000000 \
#      -use_lipinski_filter TRUE -use_strict_lipinski_filter TRUE -use_ghose_filter TRUE \
#      -scoring_function VINA -score_by_ligand_efficiency FALSE \
#      -maintain_core FALSE -minimum_core_atoms_required 0 \
#      -num_processors 4 #\
#     #  -output_dir ./autogrow_output/ #\
#      #-vina_executable /PATH/TO/VINA/EXECUTABLE/vina 
#     #  -openbabel_bin_directory /PATH/TO/OPENBABEL/BIN/DIR/bin/ \
#     #  -mgltools_directory /PATH/TO/MGLTOOLS/DIRECTORY/MGLTools-1.5.4/ \
pwd