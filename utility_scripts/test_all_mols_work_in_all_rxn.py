"""
This script is used to test all compounds within a complimentary molecule library 
work with all reactions that may call them.
"""

import __future__
import glob
import os
import sys
import pickle
import copy
import argparse

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import support_scripts.Multiprocess as mp


num_processors = -1

def dothing(file):
    zinc_list = []
    smile_list = []


    with open(file, "r") as f:
        for line in f.readlines():


            line = line.replace("\n","")
            parts = line.split("\t")
            if len(parts) != 2:
                parts = line.split("    ")
                if len(parts) != 2:
                    parts = line.split(" ")
                    if len(parts) != 2:
                        print("Fail because line could split in {}".format(file))
                        print(line)
                        print("")
                        raise Exception("Fail because line could split in {}".format(file))
            smile = parts[0]
            zinc_id = parts[1]

            zinc_list.append(zinc_id)
            smile_list.append(smile)


    if len(smile_list)<50 or len(zinc_list)<50:
        print("FAIL NO LIGANDS!!!!!!!")
        raise Exception("FAIL NO LIGANDS!!!!!!!")
    job_input = [[smile_list[i], zinc_list[i]]  for i in range(0, len(smile_list))]

    output = mp.multi_threading(job_input, num_processors,  dothing_to_mol)

    output = [x for x in output if x is not None]
    if len(output)==0:
        return None

    else:
        for i in output:
            if i is None:
                continue
            else:
                print(i)

    return file
    #

def dothing_to_mol(smile,zinc_id):
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)

    return [zinc_id, [smile,zinc_id, mol]]


def get_mols(file_name):

    job_input = []
    # counter = 0
    with open(file_name, "r") as f:

        for line in f.readlines():

            # # FOR DEBUGGING
            # counter= counter + 1
            # if counter >100:
            #     break

            line = line.replace("\n","")
            parts = line.split("\t")
            if len(parts) != 2:
                parts = line.split("    ")
                if len(parts) != 2:
                    parts = line.split(" ")
                if len(parts) != 2:
                    print("Fail because line could split in {}".format(file_name))
                    print(line)
                    print("")
                    raise Exception("Fail because line could split in {}".format(file))

            # counter = counter +1
            job_input.append((parts[0], parts[1]))
            continue
    print("Multithread testing: ", file_name)
    print("      len multithread:", len(job_input))
    # job_input = [[smile_list[i], zinc_list[i]] for i in range(0, len(smile_list))]
    output = mp.multi_threading(job_input, num_processors,  dothing_to_mol)

    mol_dic = {}
    for x in output:
        mol_dic[x[0]] = x[1]

    return mol_dic

##############################
def conduct_reaction_one(mol_set, rxn, sub, substructure_from_rxn):
    mol = mol_set[2]

    if mol.HasSubstructMatch(sub) is False:
        mol = Chem.RemoveHs(mol)
        if mol.HasSubstructMatch(sub) is False:
            return ["Missing", mol_set]
    if mol.HasSubstructMatch(substructure_from_rxn) is False:
        mol = Chem.RemoveHs(mol)
        if mol.HasSubstructMatch(substructure_from_rxn) is False:
            return ["Missing", mol_set]

    try:
        rxn.RunReactants((mol,))[0][0]
        return  [True, mol_set]
    except:
        return [False, mol_set]

def conduct_reaction_two_mol2control(mol_set, mol_2, rxn, sub_1, substructure_from_rxn):

    mol_1 = mol_set[2]
    if mol_1.HasSubstructMatch(sub_1) is False:
        mol_1 = Chem.RemoveHs(mol_1)
        if mol_1.HasSubstructMatch(sub_1) is False:
            mol_1 = Chem.AddHs(mol_1)
            if mol_1.HasSubstructMatch(sub_1) is False:
                return ["Missing", mol_set]
    if mol_1.HasSubstructMatch(substructure_from_rxn) is False:
        mol_1 = Chem.RemoveHs(mol_1)
        if mol_1.HasSubstructMatch(substructure_from_rxn) is False:
            mol_1 = Chem.AddHs(mol_1)
            if mol_1.HasSubstructMatch(substructure_from_rxn) is False:
                return ["Missing", mol_set]
    try:
        rxn.RunReactants((mol_1, mol_2))[0][0]
        return [True, mol_set]
    except:
        try:
            mol_2 = Chem.AddHs(mol_2)
            rxn.RunReactants((mol_1,mol_2))[0][0]
            return [True, mol_set]
        except:
            return [False, mol_set]

def conduct_reaction_two_mol1control(mol_set, mol_1, rxn, sub_2, substructure_from_rxn):
    mol_2 = mol_set[2]

    if mol_2.HasSubstructMatch(sub_2) is False:
        mol_2 = Chem.RemoveHs(mol_2)
        if mol_2.HasSubstructMatch(sub_2) is False:
            return ["Missing", mol_set]
    if mol_2.HasSubstructMatch(substructure_from_rxn) is False:
        mol_2 = Chem.RemoveHs(mol_2)
        if mol_2.HasSubstructMatch(substructure_from_rxn) is False:
            return ["Missing", mol_set]

    try:
        rxn.RunReactants((mol_1, mol_2))[0][0]
        return [True, mol_set]
    except:
        try:
            mol_2 = Chem.AddHs(mol_2)
            rxn.RunReactants((mol_1,mol_2))[0][0]
            return [True, mol_set]
        except:
            return [False, mol_set]

def write_to_badmol_library(missing_substructure,failed_reaction, output_folder, functional_group_name):


    # missing_substructure
    missing_sub_folder = output_folder + "missing_sub/"
    if os.path.isdir(missing_sub_folder) is False:
        os.makedirs(missing_sub_folder)
    missing_sub_file = missing_sub_folder + functional_group_name + ".smi"
    printout = ""
    for missing_mol in missing_substructure:
        temp_list = [missing_mol[1][0],missing_mol[1][1]]
        mol_info = "\t".join(temp_list)
        printout = printout + mol_info + "\n"
    if os.path.exists(missing_sub_file) is True:
        with open(missing_sub_file, "a") as f:
            f.write(printout)
    else:
        with open(missing_sub_file, "w") as f:
            f.write(printout)

    # failed_reaction
    printout = ""
    failed_reaction_folder = output_folder + "failed_reaction/"
    if os.path.isdir(failed_reaction_folder) is False:
        os.makedirs(failed_reaction_folder)
    failed_reaction_file = failed_reaction_folder + functional_group_name + ".smi"

    for failed_mol in failed_reaction:
        temp_list = [failed_mol[1][0],failed_mol[1][1]]
        mol_info = "\t".join(temp_list)
        printout = printout + mol_info + "\n"

    if os.path.exists(failed_reaction_file) is True:
        with open(failed_reaction_file, "a") as f:
            f.write(printout)
    else:
        with open(failed_reaction_file, "w") as f:
            f.write(printout)


def write_to_goodmol_library(mols_dict, final_output_folder, functional_group_name):


    # missing_substructure
    final_name = final_output_folder +functional_group_name + ".smi"

    printout = ""
    for key in list(mols_dict.keys()):
        temp = [mols_dict[key][0],mols_dict[key][1]]
        printout = printout + "\t".join(temp) + "\n"

    with open(missing_sub_file, "w") as f:
        f.write(printout)

def get_mols_dict_from_pickle(file_path):

    with open(file_path, 'rb') as handle:
        mols_dict = pickle.load(handle)
    return mols_dict

def write_pickle_to_file(file_path, obj):
    with open(file_path, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)


def write_final_outputs(out_file, pickle_file):
    mols_dict = get_mols_dict_from_pickle(pickle_file)

    printout = ""
    with open(out_file, "w") as f:
        f.write(printout)

    with open(out_file, "a") as f:
        for zinc_id in mols_dict:
            printout = ""
            temp = [mols_dict[zinc_id][0], mols_dict[zinc_id][0]]
            printout = printout + "\t".join(temp) + "\n"
            f.write(printout)
            printout = ""
    mols_dict = None

def check_id_in_list(zinc_id, list_to_remove):
    """
    if in list return None (will be excluded)
    if NOT IN LIST RETURN zinc_id
    """
    if zinc_id in list_to_remove:
        return None
    else:
        return zinc_id

def run_post_multithread(output, functional_group, output_folder, modified_pickle_folder):
    """
    if in list return None (will be excluded)
    if NOT IN LIST RETURN zinc_id
    """
    missing_substructure = [x for x in output if x[0]=="Missing"]
    failed_reaction = [x for x in output if x[0] is False]
    passed_reaction = [x[1] for x in output if x[0] is True]
    output = []
    print("            {} ligands are missing_substructure.".format(len(missing_substructure)))
    print("            {} ligands failed reaction.".format(len(failed_reaction)))
    print("            {} ligands which PASSED Everything!.".format(len(passed_reaction)))

    print("     Writing bad mols to file")

    # write the failures to a bad mol library
    write_to_badmol_library(missing_substructure,failed_reaction, output_folder, functional_group)
    missing_substructure = []
    failed_reaction = []


    print("     Making new mol_dict without bad mols")
    if len(passed_reaction)!=0:
        print("example of passed_reaction:  ", passed_reaction[0])
    mols_dict = {}
    for x in passed_reaction:
        mols_dict[x[1]] = x

    print("     Finished making modified mols_dict")
    # Repickle the mols_dict folder to modified.
    temp_pickle_name = modified_pickle_folder + functional_group + "_pickle"
    print("     Pickle modified mols_dict")
    write_pickle_to_file(temp_pickle_name, mols_dict)
    # Modify the pickled_file_dictionary to use the modified pickle dictionary
    mols_dict = None

    print("     Finished with run_post_multithread")
    return temp_pickle_name




if __name__ == "__main__":

    rxn_set = "Robust"
    rxn_set = "all_rxns"
    if rxn_set=="Robust": reactions=Robust_reactions
    elif rxn_set=="AUTOCLICKCHEM": reactions=AUTOCLICKCHEM_reactions
    elif rxn_set=="all_rxns": reactions=All_Rxns_reactions
    else:
        raise Exception("WHICH REACTION SET?")


    folder = "/home/jacob/Desktop/test_all_rxns/complimentary_mol_dir/"

    output_folder = "/home/jacob/Desktop/test_all_rxns/Robust_reactions_FILTERED/"
    source_pickle_folder = "/home/jacob/Desktop/test_all_rxns/Robust_reactions_FILTERED/source_pickled_lib/"
    modified_pickle_folder = "/home/jacob/Desktop/test_all_rxns/Robust_reactions_FILTERED/modified_pickled_lib/"
    good_mols_folder = "/home/jacob/Desktop/test_all_rxns/Robust_reactions_FILTERED/Final/"

    if os.path.exists(folder) is False:
        raise Exception("folder HERE!!")

    if os.path.exists(output_folder) is False:
        os.mkdir(output_folder)

    if os.path.exists(source_pickle_folder) is False:
        os.mkdir(source_pickle_folder)
    if os.path.exists(modified_pickle_folder) is False:
        os.mkdir(modified_pickle_folder)
    if os.path.exists(good_mols_folder) is False:
        os.mkdir(good_mols_folder)


    list_file_basename = [os.path.basename(x).replace(".smi","") for x in glob.glob(folder+"*.smi")]


    print("Geting molecules from library")
    pickled_file_dictionary = {}

    # Make a dictionary of pickle file locations

    for name in list_file_basename:
        # 1st see if we've made a modified pickle file before
        # This would be a file which has already been filtered down.
        # # If so lets use this

        modified_pickle_file = modified_pickle_folder + name+"_pickle"
        if os.path.exists(modified_pickle_file) is True:
            pickled_file_dictionary[name] = modified_pickle_file
            modified_pickle_file = None
        # If no pre-modified pickle file exists lets check if there is a pickled version of the source file.
        # Presumably this is bigger than the modified would've been
        else:
            temp_pickle_name = source_pickle_folder + name+"_pickle"
            # print(temp_pickle_name)
            pickled_file_dictionary[name] = temp_pickle_name

            if os.path.exists(temp_pickle_name) is True:
                # print("pickled version Already Exists use this")
                temp_pickle_name = None
                continue

            else:
                # Import all molecules from a large list. This often takes a while.
                file_name = folder + name + ".smi"
                temp_mol = get_mols(file_name)
                write_pickle_to_file(temp_pickle_name, temp_mol)
                temp_mol = None
                temp_pickle_name = None

    #Close out excess variables
    list_file_basename = None

    for top_level_key in list(reactions.keys()):
        sub = reactions[top_level_key]

        reaction_string =sub["reaction_string"]
        functional_groups = sub["functional_groups"]
        group_smarts = sub["group_smarts"]
        example_rxn_reactants = sub["example_rxn_reactants"]
        num_reactants = sub["num_reactants"]
        reaction_string_split = reaction_string.split(">>")[0]

        rxn_num = sub["RXN_NUM"]


        print("")
        print("")
        print("Running Reaction: ", top_level_key)
        print("")

        rxn = AllChem.ReactionFromSmarts(reaction_string)
        rxn.Initialize()


        list_mols_to_react_1 = []
        mol_standardized = []
        mols = []


        if num_reactants == 1:
            print("     Running onestep Reactions for: ", top_level_key)

            # Get mols_dict from pickled file
            mols_dict = get_mols_dict_from_pickle(pickled_file_dictionary[functional_groups[0]])
            for key in list(mols_dict.keys()):
                mols.append(mols_dict[key])
            mols_dict = {}

            substructure = Chem.MolFromSmarts(group_smarts[0])
            substructure_from_rxn = Chem.MolFromSmarts(reaction_string_split)
            job_input = tuple([tuple([mols[i], rxn, substructure, substructure_from_rxn]) for i in range(0, len(mols))])
            mols = []
            substructure = None
            rxn = None
            print("     multi_threading onestep")
            print("          {} Number of reactions to perform.".format(len(job_input)))
            output = mp.multi_threading(job_input, num_processors,  conduct_reaction_one)
            job_input = None
            output = [x for x in output if x is not None]

            pickled_file_dictionary[functional_groups[0]] = run_post_multithread(output, functional_groups[0], output_folder, modified_pickle_folder)

            output = None


        if num_reactants == 2:

            print("     Running Two Step Reactions for: ", top_level_key)

            reaction_string_split = reaction_string_split.split(".")


            # Handle mol_1 as variable and mol_2 as a control.
            # Get mols_dict from pickled file
            mols_dict = get_mols_dict_from_pickle(pickled_file_dictionary[functional_groups[0]])

            control_mol = Chem.MolFromSmiles(example_rxn_reactants[1])
            sub_control = Chem.MolFromSmarts(group_smarts[1])
            if control_mol.HasSubstructMatch(sub_control) is False:
                control_mol = Chem.AddHs(control_mol)
                if control_mol is None:
                    print("THIS SHOULDNT HAPPEN {}".format(reaction_string))

                    raise Exception("THIS SHOULDNT HAPPEN")
                if control_mol.HasSubstructMatch(sub_control) is False:
                    print("THIS SHOULDNT HAPPEN 2 {}".format(reaction_string))
                    raise Exception("THIS SHOULDNT HAPPEN 2")
                else:
                    control_mol = Chem.AddHs(control_mol)
            sub_control = None


            for key in list(mols_dict.keys()):
                mols.append(mols_dict[key])

            substructure = Chem.MolFromSmarts(group_smarts[0])
            substructure_from_rxn = Chem.MolFromSmarts(reaction_string_split[0])
            job_input = tuple([tuple([mols[i], control_mol, rxn, substructure,substructure_from_rxn]) for i in range(0, len(mols))])
            mols_dict = {}
            mols = []
            substructure = None
            print("     multi_threading mol_1 as variable, mol_2 as control")
            print("          {} Number of reactions to perform.".format(len(job_input)))

            output = mp.multi_threading(job_input, num_processors,  conduct_reaction_two_mol2control)
            job_input = None
            print("     Removing bad mols")

            output = [x for x in output if x is not None]

            pickled_file_dictionary[functional_groups[0]] = run_post_multithread(output, functional_groups[0], output_folder, modified_pickle_folder)
            output = None
            ###################################################################
            ###################################################################
            ###################################################################
            ###################################################################

            # Handle mol_2 as variable and mol_1 as a control.
            # Get mols_dict from pickled file
            mols_dict = get_mols_dict_from_pickle(pickled_file_dictionary[functional_groups[1]])

            control_mol = Chem.MolFromSmiles(example_rxn_reactants[0])
            sub_control = Chem.MolFromSmarts(group_smarts[0])
            if control_mol.HasSubstructMatch(sub_control) is False:
                control_mol = Chem.AddHs(control_mol)
                if control_mol is None:
                    print("THIS SHOULDNT HAPPEN {}".format(reaction_string))
                    raise Exception("THIS SHOULDNT HAPPEN")
                if control_mol.HasSubstructMatch(sub_control) is False:
                    print("THIS SHOULDNT HAPPEN 2 {}".format(reaction_string))
                    raise Exception("THIS SHOULDNT HAPPEN 2")
                else:
                    control_mol = Chem.AddHs(control_mol)
            sub_control = None

            mols = []
            for key in list(mols_dict.keys()):
                mols.append(mols_dict[key])
            mols_dict = {}

            substructure = Chem.MolFromSmarts(group_smarts[1])
            substructure_from_rxn = Chem.MolFromSmarts(reaction_string_split[1])

            # mol_set, mol_1, rxn, sub_2, substructure_from_rxn
            job_input = tuple([tuple([mols[i], control_mol, rxn, substructure, substructure_from_rxn]) for i in range(0, len(mols))])
            mols = None
            substructure = None
            rxn = None
            print("     multi_threading mol_2 as variable, mol_1 as control")
            print("          {} Number of reactions to perform.".format(len(job_input)))
            output = mp.multi_threading(job_input, num_processors,  conduct_reaction_two_mol1control)
            job_input = None
            print("     Removing bad mols")

            pickled_file_dictionary[functional_groups[1]] = run_post_multithread(output, functional_groups[1], output_folder, modified_pickle_folder)

            output = None



    print("")
    print("")
    print("")
    print("#################")

    print("Finished Testing Reactions!!!")

    # Save Final sets
    if os.path.isdir(good_mols_folder) is False:
        os.makedirs(good_mols_folder)


    job_input = []
    for groups in list(pickled_file_dictionary.keys()):
        out_file = good_mols_folder + groups + ".smi"
        pickle_file = pickled_file_dictionary[groups]
        temp = (out_file, pickle_file)
        job_input.append(temp)


    output = mp.multi_threading(job_input, -1,  write_final_outputs)


    print("finished!")
