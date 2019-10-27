import __future__


def Run_Rank_selector(usable_list_of_smiles, number_to_chose, column_idx_to_select, reverse_sort=False):
    """
    Given a data set and an idx number to select based on it will select the top rank scores for that critera.
    The number is choses is defined by number_to_chose. 

    This is an alternative to the weight roulette style selectors. 

    Input:
    :param list usable_list_of_smiles: a list with all the information of all the mols in the previous generation
    :param int number_to_chose: the number of molecules to chose based on diversity score
    :param int column_idx_to_select: the idx to use as the criteria for selection.
                                    In the case of docking affinity score column_idx_to_select is -2
                                    For diversity score column_idx_to_select is -1
    :param bol reverse_sort:    Set to True if you want to select the most positive number is the best choice
                                Set to False if you want to select the most negative number

    Returns:
    :returns: list top_choice_smile_order: list of ligands chosen by a elitism selection, without replacement, 
    """
    if type(usable_list_of_smiles) is not type([]):
        raise Exception("usable_list_of_smiles Must be a list, wrong data type")
        
    num_ligands = len(usable_list_of_smiles)  
    if num_ligands == 0:
        raise Exception("usable_list_of_smiles is an empty list. There is nothing to chose from.")

    if number_to_chose <= 0:
        top_choice_smile_order = []
        return top_choice_smile_order

    # Sort by chosen idx property
    sorted_list = sorted(usable_list_of_smiles, key=lambda x: float(x[column_idx_to_select]), reverse=reverse_sort)

    if len(sorted_list) < number_to_chose:

        raise Exception("There are more ligands to chose to seed the list than ligands to select from. \
                        Please lower the top_mols_to_seed_next_generation and/or diversity_mols_to_seed_first_generation")

    top_choice_smile_order = []
    for i in range(0, number_to_chose):
        smile = sorted_list[i]
        top_choice_smile_order.append(smile[0])
        
    return top_choice_smile_order

#
