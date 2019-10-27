import __future__

import random
import math
import copy


def Run_Tournement_Selector(list_of_ligands, num_to_chose, tourn_size, idx_to_sel, favor_most_negative=True):
    """
    This runs a tournement style selector given a list of ligands and specified metric. It will randomly 
    select ligands for tournements. The best scoring ligand for each of these groups will end up in the
    chosen_ligands list. 

    This is done WITHOUT REPLACEMENT. This does provide an opportunity for any ligand to make 
    it into the chosen_ligands list even if it doesn't have a high score,
    but that chance of random incorporation decreases as tourn_size increases.
        -ie tourn_size=1.0 will randomly pick N number of ligands equal to the total number of ligands in the list
            this means theres a high chance that the top ligand will be chosen enter every tournement 
            and will win everytime. This could result in a very homogenous choice. 

    :param list list_of_ligands: The list of lists containing info about ligands with scores to select from.
    :param int num_to_chose: the number of ligands to be chosen total 
                            this also is the number of tournements that will be conducted.
    :param float tourn_size: percentage of the total pool of ligands to be tested in each tournement.
    :param int idx_to_sel: the idx within each sublist which will serve as the metric for each tournement.
    :param bol favor_most_negative: True if the most negative number is the best solution. 
                                    False if the most positive number is the best solution
                                    default to True.
    Returns:
    :returns: list chosen_ligands: a list of chosen ligands containing all the info for each ligand
                                    with potential for redundancy

    """
    
        
    if type(list_of_ligands) is not type([]):
        raise Exception("list_of_ligands Must be a list, wrong data type")
        
    num_ligands = len(list_of_ligands)  
    if num_ligands == 0:
        raise Exception("list_of_ligands is an empty list. There is nothing to chose from.")

    if idx_to_sel != -1:
        if len(list_of_ligands[0]) < idx_to_sel:
            raise Exception("The idx to select by does not exist in the provided list_of_ligand.")

    num_per_tourn = int(math.ceil(num_ligands * tourn_size))

    chosen_ligands = []
    list_of_ligands_reduced = copy.deepcopy(list_of_ligands)
    for i in range(0, num_to_chose):
        chosen_ligand = run_one_tournement(list_of_ligands, num_per_tourn, idx_to_sel, favor_most_negative)
        list_of_ligands_reduced = [x for x in list_of_ligands_reduced if x!=chosen_ligand]
        chosen_ligands.append(chosen_ligand)

    return chosen_ligands   
#

def run_one_tournement(list_of_ligands, num_per_tourn, idx_to_sel, favor_most_negative=True):
    """
    This runs a single tournement style selection given a list of ligands and specified metric. It will randomly 
    select ligands for the tournement. The best scoring ligand from the tournement will be returned.

    This is done WITHOUT REPLACEMENT. This does provide an opportunity for any ligand to make it \into the chosen_ligands list even if it doesn't have a high score,
    but that chance of random incorporation decreases as tourn_size increases.
        -ie tourn_size=1.0 will randomly pick N number of ligands equal to the total number of ligands in the list
            this means theres a high chance that the top ligand will be chosen enter every tournement 
            and will win everytime. This could result in a very homogenous choice. 

        -num_per_tourn is the int(math.ceil(num_ligands * tourn_size)) so that it rounds to the nearest int
            with a minimum values of 1


    :param list list_of_ligands: The list of lists containing info about ligands with scores to select from.
    :param int num_per_tourn: the number of ligands to be tested in each tournement.
    :param int idx_to_sel: the idx within each sublist which will serve as the metric for each tournement.
    :param bol favor_most_negative: True if the most negative number is the best solution. 
                                    False if the most positive number is the best solution
                                    default to True.
    Returns:
    :returns: list chosen_option: a list with a single ligand chosen from a single tournement
    """
    num_ligands = len(list_of_ligands)
    
    chosen_option = []
    temp = []
    for i in range(0,num_per_tourn):
        temp.append(i)
        if i == 0:
            chosen_option = list_of_ligands[random.randint(0, num_ligands-1)]
        else:
            choice = list_of_ligands[random.randint(0, num_ligands-1)]
            if favor_most_negative == True:
                if float(chosen_option[idx_to_sel]) > float(choice[idx_to_sel]):
                    chosen_option = choice    
                else:
                    continue        
            elif favor_most_negative == False:
                if float(chosen_option[idx_to_sel]) < float(choice[idx_to_sel]):
                    chosen_option = choice
                else:
                    continue

    return chosen_option
#
