"""
Dictionary and dictionary handling functions for molecular operations.

This module provides utility functions for handling dictionaries related to
molecular structures, R-groups, and B-groups in the context of molecular
crossover operations.
"""
import __future__

import copy
from typing import Dict, List

import rdkit  # type: ignore
from rdkit import Chem  # type: ignore

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.plugins.crossover.merge_functions.mapping_class as mapping_class


def handle_dicts_and_select_b_groups(
    mol_1: Chem.rdchem.Mol, mol_2: Chem.rdchem.Mol, mcs_mol: Chem.rdchem.Mol
):
    """
    Create necessary dictionaries, mappings, and select ligands for the final
    molecule.

    Args:
        mol_1 (Chem.rdchem.Mol): RDKit mol for ligand 1.
        mol_2 (Chem.rdchem.Mol): RDKit mol for ligand 2.
        mcs_mol (Chem.rdchem.Mol): RDKit mol for shared common core between
            mol_1 and mol_2.

    Returns:
        List[str]: SMILES strings for the R groups corresponding to the chosen
            B's. Returns None if it fails.
    """
    # Confirm that mcs_mol can be replaced in mol_1 and mol_2 around 0.8% of
    # the time this function fails so we will filter this 1st
    they_pass = _check_replace_mol(mol_1, mol_2, mcs_mol)
    if they_pass is False:
        return None

    (
        r_smiles_dict_1,
        b_to_r_master_dict_1,
        b_to_anchor_master_dict_1,
    ) = _mol_handling_of_fragmenting_labeling_and_indexing(mol_1, mcs_mol, 1)

    # check that this worked (ie if it failed they will return None)
    if r_smiles_dict_1 is None:
        return None
    if b_to_r_master_dict_1 is None:
        return None
    if b_to_anchor_master_dict_1 is None:
        return None

    (
        r_smiles_dict_2,
        b_to_r_master_dict_2,
        b_to_anchor_master_dict_2,
    ) = _mol_handling_of_fragmenting_labeling_and_indexing(mol_2, mcs_mol, 2)

    # check that this worked (ie if it failed they will return None)
    if r_smiles_dict_2 is None:
        return None
    if b_to_r_master_dict_2 is None:
        return None
    if b_to_anchor_master_dict_2 is None:
        return None

    # Merge b_to_anchor_master_dict into 1 master dictionary of B_to_anchors.
    # the keys will be all be unique so we can add these dictionaries together
    # without worry of overrighting an entry. We will invert the dict after to
    # get anchors as the keys and the B's as the items. example
    # b_to_anchor_master {'1B1':[10008,10007],'1B2':[10000],'1B3':[10006],
    # '2B3':[10006,10007],'2B2':[10000],'2B1':[10008]}
    b_to_anchor_master = b_to_anchor_master_dict_1
    for i in list(b_to_anchor_master_dict_2.keys()):
        b_to_anchor_master[i] = b_to_anchor_master_dict_2[i]

    # Invert b_dictionary to produce a master I dictionary. example
    # anchor_to_b_master = {10008:['1B1','2B1'],10000:['1B2','2B2']}
    anchor_to_b_master = _invert_dictionary(b_to_anchor_master)

    bs_chosen = mapping_class.run_mapping(b_to_anchor_master, anchor_to_b_master)

    # Get the R groups which correspond to the chosen B's
    # ['1R1', '1R5', '2R2']
    rs_chosen = _get_rs_chosen_from_bs(
        bs_chosen, b_to_r_master_dict_1, b_to_r_master_dict_2
    )

    # Get the smiles strings for the the R groups which correspond to the
    # chosen B's
    return _get_rs_chosen_smiles(rs_chosen, r_smiles_dict_1, r_smiles_dict_2)


def _mol_handling_of_fragmenting_labeling_and_indexing(mol, mcs_mol, lig_number):
    """
    This takes an rdkit mol for a ligand and 1 for the mcs_mol. It fragments
    the ligand by replacing the MCS. and it determines which anchors are in
    each fragment. These fragments are our R-groups and the assignment of
    anchors. is how we determine which R-group goes where relative to the MCS.

    lig_number  int    is the number of the ligand that is mol
                       ie if mol is mol_1 lig_number = 1
                       ie if mol is mol_2 lig_number = 2


    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit mol (either mol_1 or mol_2)
    :param rdkit.Chem.rdchem.Mol mcs_mol: rdkit mol for shared common core
        between mol_1 and mol_2
    :param int lig_number: an int either 1 or 2 for (mol_1 or mol_2
        respectively)

    Returns:
    :returns: dict r_smiles_dictionary: a dictionary of the R-groups which
        branch off the common core keys are the R-groups; items are the SMILES
        strings of that R-groups returns None if it fails. Example: {'1R1':
        '[10003*][1007N]=[1013O]', '1R2': '[10000*][1011CH2]=[1008O]'}
    :returns: dict b_to_r_master_dict: A dictionary which tracks the R groups
        which belong to a B-group keys are the B-groups; items are the R-groups
        which belong to the B-group. returns None if it fails. Example: {'1B1':
        ['1R2'], '1B2': ['1R1']}
    :returns: dict b_to_anchor_master_dict: A dictionary which tracks the iso
        label of the anchor atoms for B-group. keys are the B-groups; items are
        the iso label of the anchor atoms for B-group returns None if it fails.
        Example:{'1B1': [10000], '1B2': [10003]}
    """
    # Find which MCS Atoms Rs branch from Function to find all neighbors for a
    # set of molecules touching an Isolabeled core
    mcs_touches = _get_atoms_touch_mcs(mol)

    # invert dictionary
    lig_r_atoms_touch_mcs = _invert_dictionary(mcs_touches)

    # remove the Core atoms from each ligand this gives us the R-groups
    replace_core = _r_group_list(mol, mcs_mol)
    if replace_core is None:
        # replace_core failed to handle fragments"
        return None, None, None

    replace_core = _replace_core_mol_dummy_atoms(mol, mcs_mol, replace_core)
    if replace_core is None:
        # replace_core failed to handle fragments"
        return None, None, None

    # A single anchor (isotope label) may now be present multiple times in the
    # replace_core_mols as they are fragmented replace_frag_w_anchor_isolabels
    # can return a None if failed so lets check for None before we move on
    if replace_core is None:
        # replace_core failed to handle fragments"
        return None, None, None

    # MAKE NEW MOL FRAGS FROM LABELED replace_core
    mol_frags = Chem.GetMolFrags(replace_core, asMols=True, sanitizeFrags=False)
    list_r_groups = []
    i = 0
    while i < len(mol_frags):
        val = Chem.MolToSmiles(mol_frags[i], isomericSmiles=True)
        list_r_groups.append(val)
        i += 1

    # Generate all the R-libraries with full R-groups using the index of its
    # respective Lig r_chain_dictionary is the master dictionary for R-groups
    r_chain_dictionary, r_smiles_dictionary = _r_groups_dict(mol_frags, lig_number)

    # r_dict is a secondary dictionary for searching I's in R's. this
    # dictionary is limited to only the R-group and anchor(I).
    r_dict = _get_r_dict(r_chain_dictionary, lig_r_atoms_touch_mcs)

    # make inversion of r_dict. keys are the Anchor atom iso_labels while the
    # items are the R-group numbers which are attached to that anchor atom.
    # Example: {10008: ['2R3'], 10000: ['2R2'], 10006: ['2R1'], 10007:
    # ['2R1']}
    i_dict = _invert_dictionary(r_dict)

    """
    B-dictionaries:
    Ligmerge will randomly select R-groups to append to a shared common core
    from 2 separate ligands but so anchor atoms in the common core may have
    more than 1 R-group attached to it.

    ie. if an anchor carbon has di-methyls attached to it (which are not part
        of the shared core) are these di-methyls 2 separate R-groups or is the
        contextual chemical environment created by having both methyls
        different from having one alone and thus should be treated as a single
        R-group which just happen to branch. This author would argue that
        context is important here and so this version of Ligmerge treats
        anything attached to an anchor atom in the common core as a singular
        contextual functional group which shall be referred to as a B-groups.
    ie. a B-group consists of 1 or more R-groups which are attached to an
        anchor atom in the shared common core. This makes a significant
        difference in how we select for which pieces are added to build our
        child molecule. Additionally this has significance in the decision
        tree use to build a child molecule. An R/B group can be connected to
        multiple anchor atoms so once we chose a B group we will need to know
        which anchor atoms are affected by that decision. This is something
        handled more in the Mapping class, but this is why the nomenclature
        change from R-groups to B-groups and why the next several steps are
        important.
    make_b_dictionaries (B is the name we gave to R-groups sets)
    """
    b_to_r_master_dict, b_to_anchor_master_dict = _make_b_dic(
        i_dict, r_dict, lig_number
    )

    return r_smiles_dictionary, b_to_r_master_dict, b_to_anchor_master_dict


def _check_replace_mol(mol_1, mol_2, mcs_mol):
    """
    Confirm that mcs_mol can be replaced in mol_1 and mol_2.

    Args:
        mol_1 (Chem.rdchem.Mol): An RDKit mol.
        mol_2 (Chem.rdchem.Mol): An RDKit mol.
        mcs_mol (Chem.rdchem.Mol): RDKit mol for shared common core between
            mol_1 and mol_2.

    Returns:
        bool: True if it passes for both mol_1 and mol_2, False if either fails.
    """
    temp = _r_group_list(mol_1, mcs_mol)
    if temp is None:
        return False
    temp = _r_group_list(mol_2, mcs_mol)
    return temp is not None


###########################################################
###########################################################
# HANDLE THE OBTAINING THE R-Groups for a given mol


def _r_group_list(mol, core_mol):
    """
    Find all R-groups by replacing the atoms in the ligand that make up the
    common core with nothing.

    This fragments the ligand and from those fragments we are able to determine
    what our R-groups are. For any common core atom which touched the fragment,
    a * will replace that atom in the fragments.

    Args:
        mol (Chem.rdchem.Mol): An RDKit molecule.
        core_mol (Chem.rdchem.Mol): An RDKit molecule for the shared common
            core.

    Returns:
        Chem.rdchem.Mol: An RDKit molecule (replace_core_mol) with the common
            core removed from a ligand. This fragments the mol which can be
            used to make lists of R-groups. Returns None if the mol either did
            not contain the core_mol or the core_mol is the same mol as the
            mol. This is rare but does occur when the only difference is H's,
            which means it can't be replaced within because it's the same mol.
    """
    # This returns all the mol frags for a particular compound against the
    # core molecule
    replace_core_mol = Chem.ReplaceCore(
        mol, core_mol, labelByIndex=True, replaceDummies=True, requireDummyMatch=False
    )

    if len(replace_core_mol.GetAtoms()) == 0:
        # This means that the mol either did not contain the core_mol or the
        # core_mol is the same mol as the mol. ie) if mol_string
        # ="[10000N-]=[10001N+]=[10002N][10003CH]1[10004O][10005CH]([10006CH2][10007OH])[10008CH]([10013OH])[10009CH]([10012OH])[10010CH]1[10011OH]"
        # and core_string
        # ="[10000NH]=[10001N+]=[10002N][10003CH]1[10004O][10005CH]([10006CH2][10007OH])[10008CH]([10013OH])[10009CH]([10012OH])[10010CH]1[10011OH]"
        # the only difference is the H's which means it can be replaced within
        # because its the same mol This is rare but does occur.
        return None

    return replace_core_mol


def _replace_core_mol_dummy_atoms(mol, mcs, replace_core_mol):
    """
    Replace the dummy atoms (*) with the isotope label from the core atoms.

    This function will replace the dummy atoms (*) with the isotope label from
    the core atoms.

    Args:
        mol (Chem.rdchem.Mol): An RDKit molecule.
        mcs (Chem.rdchem.Mol): An RDKit molecule for the shared common core.
        replace_core_mol (Chem.rdchem.Mol): The mol with the MCS anchors
            labeled with * and an isotope label of the idx of the core anchor
            atom.

    Returns:
        Chem.rdchem.Mol: An RDKit molecule with the common core removed from a
            ligand fragments the mol which can be used to make lists of
            R-groups. The * atoms will be isotope labeled with the isotope
            label from the core.

    Example:
        mol = Chem.MolFromSmiles("[10000N-]=[10001N+]=[10002N][10003CH2][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        mcs = Chem.MolFromSmiles("[10003CH3][10002N]=[10001N+]=[10000NH]")
        replace_core = Chem.MolFromSmiles("[3*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")

        resulting replace_core = '[10003*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]'
    """
    replace_core_mol_original = copy.deepcopy(replace_core_mol)
    anchor_dict = {}
    anchor_to_set_dict = {}
    for atom in replace_core_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            anchor_iso = atom.GetIsotope() + 10000
            neighbors = atom.GetNeighbors()
            tmp = [n_atom.GetIsotope() for n_atom in neighbors]
            anchor_dict[anchor_iso] = tmp

            anchor_to_set_dict[atom.GetIdx()] = anchor_iso

    for idx in list(anchor_to_set_dict.keys()):

        atom = replace_core_mol.GetAtomWithIdx(idx)
        anchor_iso = anchor_to_set_dict[idx]
        atom.SetIsotope(anchor_iso)

    return replace_core_mol


def _r_groups_dict(mol_frags, lig_number_for_multiplier):
    """
    Create dictionaries of R-groups and their SMILES for given molecular fragments.

    Args:
        mol_frags (rdkit.Chem.rdchem.Mol): RDKit molecule containing fragments.
        lig_number_for_multiplier (int): Ligand number (1 for mol_1, 2 for 
            mol_2) used for labeling.

    Returns:
        tuple: Two dictionaries:
            - r_chain_dictionary (dict): R-groups and their anchor atoms.
                e.g., {'1R1':[13,14], '1R2':[21,22], '1R3':[25]}
            - r_smiles_dictionary (dict): R-groups and their SMILES strings.
                e.g., {'1R1':'[1*]:[1013c]([1020H])[1014c]([1019H])[1015c]
                ([1018H])[1016c](:[2*])[1017H]', '1R2':'[3*][1024C]([1026H])
                ([1027H])[1023N]=[1022N+]=[1021N-]', '1R3':'[4*][1025O][1029H]'}
    """
    num_frags = len(mol_frags)
    r_chain_dictionary = {}
    r_smiles_dictionary = {}
    k = int(lig_number_for_multiplier)
    for i in range(num_frags):
        frag = mol_frags[i]
        r_list_temp = []
        r_list_smiles = Chem.MolToSmiles(frag, isomericSmiles=True)
        lig_num_r_r_num = f"{k}R{i + 1}"
        for atoms in frag.GetAtoms():
            iso = atoms.GetIsotope()
            if 3000 > iso > 100:
                r_list_temp.append(iso - (1000 * k))
                atoms.SetIsotope(0)
            if iso > 3000:
                name = f"I{iso - 10000}"
                r_list_temp.append(iso)
            r_chain_dictionary[lig_num_r_r_num] = r_list_temp
            r_smiles_dictionary[lig_num_r_r_num] = r_list_smiles
    return r_chain_dictionary, r_smiles_dictionary


def _get_r_dict(r_chain_dict, lig_r_atom_touch_mcs):
    """
    Create a dictionary of R-groups and their connected anchor atoms.

    Args:
        r_chain_dict (dict): Dictionary of R-groups and their atom isolabels.
            e.g., {'1R1': [3, 4, 5, 6, 7, 8, 9, 10, 11, 10000]}
        lig_r_atom_touch_mcs (dict): Dictionary of atoms touching the core and 
            their anchors. e.g., {3: [10000]}

    Returns:
        dict: Dictionary of R-groups and their connected anchor atoms.
            e.g., {'1R1': [10000]}
    """
    r_s_dict = {}
    for key in list(r_chain_dict.keys()):
        temp_r_list = r_chain_dict[key]
        node_list = []
        for atom in temp_r_list:
            for key_id in list(lig_r_atom_touch_mcs.keys()):
                if atom == key_id:
                    node_list.extend(iter(lig_r_atom_touch_mcs[key_id]))
                    r_s_dict[key] = node_list

    return r_s_dict


##########
# Mapping functions and finding neighbors
#########
def get_idx_using_unique_iso(mol, iso_val):
    """
    Find the atom index in a molecule with a specific isotope label.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule with uniquely labeled atoms.
        iso_val (int): Isotope value to search for.

    Returns:
        int: Index of the atom with the specified isotope label, or None if not 
            found.
    """
    # TODO: Never used?
    return next(
        (atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsotope() == iso_val),
        None,
    )


def _make_b_dic(i_dictionary, r_dict_num, lig_number):
    """
    Generate dictionaries for B-groups, tracking R-groups and anchor atoms.

    Args:
        i_dictionary (dict): Dictionary of R-groups bound to nodes (I's).
            e.g., {'10008':[1R1,1R2],'10009':[1R2,1R3]}
        r_dict_num (dict): Dictionary of anchors attached to R-groups.
            e.g., {'1R1':[10008],'1R2':[10008,10009],'1R3':[10009]}
        lig_number (int): Ligand number (1 or 2).

    Returns:
        tuple: Two dictionaries:
            - b_to_r_master_dict (dict): B-groups and their represented R-groups.
                e.g., {'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3':['1R5']}
            - b_to_anchor_master_dict (dict): B-groups and their connected anchors.
                e.g., {'1B1':[10008,10007],'1B2':[10000],'1B3':[10006]}
    """
    k = lig_number
    b_to_r_master_dict = {}
    b_to_anchor_master_dict = {}
    counter = 1
    anchor_list = list(i_dictionary.keys())
    # anchor_list = [10008, 10000, 10006, 10007]

    while anchor_list:
        anchor = anchor_list[0]
        B_key = f"{k}B{counter}"
        temp_r_list = []
        temp_anchor_list = []

        for Rs in i_dictionary[anchor]:
            # example Rs in i_dictionary[anchor]: '1R1')
            temp_r_list.append(Rs)
            r_dict_i = r_dict_num[Rs]
            for I in r_dict_i:
                # example Rs in i_dictionary[anchor]: '1R1')
                temp_anchor_list.append(I)
        # remove any redundancies in the list by list(set(list_of_things))
        temp_anchor_list = list(set(temp_anchor_list))
        temp_r_list = list(set(temp_r_list))

        # make new B-group entry in the dictionaries
        b_to_r_master_dict[B_key] = temp_r_list  # This B-represents these R-groups
        b_to_anchor_master_dict[
            B_key
        ] = temp_anchor_list  # This B connects to these anchor atoms

        counter = counter + 1

        # make a list of atoms to remove in the next iteration if they are in
        # both the temp_anchor_list and anchor_list
        for i in temp_anchor_list:
            if i in anchor_list:
                anchor_list.remove(i)

    # example
    # b_to_r_master_dict:{'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3':
    # ['1R5']}. example
    # b_to_anchor_master_dict:{'1B1':[10008,10007],'1B2':[10000],'1B3':[10006]}.
    return b_to_r_master_dict, b_to_anchor_master_dict


def _invert_dictionary(old_dic):
    """
    Invert a dictionary, making values the keys and keys the values.

    Args:
        old_dic (dict): Dictionary to invert.

    Returns:
        dict: Inverted dictionary where old values are keys and old keys are 
            values.
    """
    # inverted_dic = {}
    # for k, v in old_dic.iteritems():
    # keys = inverted_dic.setdefault(v, [])
    # keys.append(k)
    values = {a for b in list(old_dic.values()) for a in b}
    values = list(values)
    return {
        new_key: [key for key, value in list(old_dic.items()) if new_key in value]
        for new_key in values
    }


def _get_atoms_touch_mcs(mol):
    """
    Find all non-core atoms touching core atoms in a molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule with isotope-labeled atoms.
            Core atoms have labels idx + 10000, non-core atoms have 
            idx + 1000 (lig_1) or idx + 2000 (lig_2).

    Returns:
        dict: Dictionary with core atom isotope labels as keys and lists of 
            touching non-core atom indices as values.
    """
    mcs_touches = {}
    all_atoms = mol.GetAtoms()

    for atom in all_atoms:
        # find all atoms in the mol which are also in the Core using iso
        # labels

        iso = atom.GetIsotope()
        if iso > 9999:
            # then its a core atom
            neighbors = atom.GetNeighbors()
            values = []

            for neighbor_atom in neighbors:
                # compile list of the Indexes of all neighbor of the atom
                Idx_neighbor = neighbor_atom.GetIdx()

                # Select for only neighbors which are not in the core using
                # iso
                iso_neighbor_x = neighbor_atom.GetIsotope()
                if iso_neighbor_x < 9999:
                    # Then this is an atom which is not in the core but
                    # touches the core
                    idx_of_neighbor = neighbor_atom.GetIdx()
                    values.append(idx_of_neighbor)
                    mcs_touches[iso] = values

    return mcs_touches


##########
# Handling after B-groups are chosen
##########
def _get_rs_chosen_from_bs(bs_chosen, b_to_r_master_dict_1, b_to_r_master_dict_2):
    """
    Get a list of R-groups based on chosen B-groups.

    Args:
        bs_chosen (list): List of chosen B-groups. e.g., ['1B1', '1B2', '2B3']
        b_to_r_master_dict_1 (dict): B to R-group dictionary for mol_1.
            e.g., {'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3':['1R5']}
        b_to_r_master_dict_2 (dict): B to R-group dictionary for mol_2.
            e.g., {'2B1':['2R1'],'2B2':['2R2','2R3','2R4'],'2B3':['2R5','2R6']}

    Returns:
        list: List of R-groups represented by the chosen B-groups.
            e.g., ['1R1', '1R2', '1R3','1R4', '2R5', '2R6']
    """
    rs_chosen = []
    for B in bs_chosen:
        Rs_for_the_B = []
        lig_number = B[0]
        B_number = B[2]
        if lig_number == str(1):
            Rs_for_the_B.extend(iter(b_to_r_master_dict_1[B]))
        elif lig_number == str(2):
            Rs_for_the_B.extend(iter(b_to_r_master_dict_2[B]))
        rs_chosen.extend(iter(Rs_for_the_B))
    # rs_chosen looks like ['1R1', '1R5', '2R2']
    return rs_chosen


def _get_rs_chosen_smiles(
    rs_chosen: List, r_smiles_dict_1: Dict, r_smiles_dict_2: Dict
):
    """
    Get SMILES strings for chosen R-groups.

    Args:
        rs_chosen (List): List of chosen R-groups. e.g., ['2R2', '1R1']
        r_smiles_dict_1 (Dict): R-group SMILES dictionary for Ligand 1.
            e.g., {'1R1': '[10006*][1009N]=[1008N+]=[1007N-]'}
        r_smiles_dict_2 (Dict): R-group SMILES dictionary for Ligand 2.
            e.g., {'2R2': '[10006*][2009OH]', '2R1': '[10003*][2007CH2][2008OH]'}

    Returns:
        List: List of SMILES strings for chosen R-groups, each as a sublist.
            e.g., [['[10006*][1009N]=[1008N+]=[1007N-]'],['[10006*][2009OH]']]
    """
    rs_chosen_smiles = []
    for R in rs_chosen:
        Rs_for_the_R = []
        lig_number = R[0]
        R_number = R[2]
        if lig_number == str(1):
            Rs_for_the_R.append(r_smiles_dict_1[R])
        elif lig_number == str(2):
            Rs_for_the_R.append(r_smiles_dict_2[R])

        rs_chosen_smiles.append(Rs_for_the_R)

    return rs_chosen_smiles
