from __future__ import absolute_import
import unittest
import os
import glob
import copy

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.rdFMCS as rdFMCS


from autogrow.Operators.Crossover.smiles_merge.MergeFunctions import Alignment_and_Breaks as ANB
import autogrow.Operators.Crossover.smiles_merge.MergeFunctions.Dict_and_R_Groups as DnR
from autogrow.Operators.Crossover.smiles_merge.MergeFunctions import MappingClass as MC
from autogrow.Operators.Crossover.smiles_merge.MergeFunctions import Merge_w_core as MWC

from autogrow.Operators.Crossover import Execute_Crossover as EC


class CrossoverTests(unittest.TestCase):
    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        
        #Variables for ANB
        frag_string = "[10000CH3][10001CH2][10002CH2][10003CH2][10004CH2][10005CH2][10006CH2][10007CH2][10008CH2][10009CH2][10010CH2][10011CH2][10012CH2][10013CH2][10014CH2][10015CH3].[10016CH3][10017CH2][10018CH3].[10019OH2]"
        self.frag = Chem.MolFromSmiles(frag_string, sanitize=False)
        self.smallest_frag_list =[19]
        self.middle_frag_list = [16, 17, 18]
        self.largest_frag_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        for atom in self.frag.GetAtoms():
            if atom.GetIsotope() != atom.GetIdx()+10000:
                printout = "This failed because frag mol has diff iso label than Idx \n"
                printout = printout + "This is used later for verifying ANB.ringbreak_frag_handling"
                raise Exception(printout)
        cyclic_break_frag_str = "[10000CH3][10001CH2][10002CH2][10003CH2][10004CH2][10005CH2][10006CH2][10007CH2][10008CH2][10009CH2][10010CH2][10011CH2][10012CH2][10013CH]([10014CH2][10015CH2][10016OH])[10017c]([10018CH3])[10019cH2].[10020CH3][10021CH2][10022CH2][10023OH].[10024OH2]"
        self.CB_frag = Chem.MolFromSmiles(cyclic_break_frag_str, sanitize=False)
        self.CB_smallest_frag_list = [10024]
        self.CB_middle_frag_list = [10020, 10021, 10022, 10023]
        self.CB_largest_frag_list = [10000, 10001, 10002, 10003, 10004, 10005, 10006, 10007, 10008, 10009, 10010, 10011, 10012, 10013, 10014, 10015, 10016, 10017, 10018, 10019]
        for atom in self.CB_frag.GetAtoms():
            if atom.GetIsotope() != atom.GetIdx()+10000:
                printout = "This failed because frag mol has diff iso label than Idx \n"
                printout = printout + "This is used later for verifying ANB.ringbreak_frag_handling"
                raise Exception(printout)
        self.mol1 = Chem.MolFromSmiles("SCCCCCCCc1ccccc1CCCCCCCC")
        self.mol2 = Chem.MolFromSmiles("NCCCCCCCc1ccccc1CCCCCCCCN")
        self.mcs_results = rdFMCS.FindMCS([self.mol1, self.mol2], matchValences = False, ringMatchesRingOnly = True, completeRingsOnly = True)
        self.mcs_mol = Chem.MolFromSmarts(self.mcs_results.smartsString)
        self.index_tuple = ((10, 9, 8, 7, 5, 6), (0, 1, 2, 3, 4, 5), (0, 1, 2, 3, 4, 5))
        self.mol_1 = Chem.MolFromSmiles('CC1=C(OC[10000CH]([10001OH])[10002CH2][10003N]=[10004N+]=[10005N-])C(C)SC1=O')
        self.mol_2 = Chem.MolFromSmiles('[10005N-]=[10004N+]=[10003N][10002CH2][10000C]1([10001OH])OC(CO)C(O)C(O)C1O')
        self.mcs_mol_2 = Chem.MolFromSmiles('[10005NH]=[10004N]=[10003N][10002CH2][10000CH2][10001OH]',sanitize=False)
        
        
        #Variables for MC
        self.Bs_to_Is = {'1B1': [10000], '1B2': [10004], '2B3': [10003, 10002], '2B2': [10002], '2B1': [10001,10003]}
        self.Is_to_Bs = {10000: ['1B1'], 10001: ['2B1'], 10002: ['2B2','2B3'], 10003: ['2B3','2B1'], 10004: ['1B2']}
        self.aMapping = MC.Mapping(self.Bs_to_Is, self.Is_to_Bs)

        
        #Variables for MWC
        # bad_SMILES_list are the wrong data type for the MOH functions. 
        # MOH functions mostly take rdkit.Chem.rdchem.Mol instead of strings
        # None of these strings can be sanitize in rdkit
        # c1ccccc1(C)(C) is invalid because it isnt viable to be an aromatic ring with a dimethlys on the same carbon.
        self.bad_SMILES_list = ["None","c1ccccc1(C)(C)", "CCC[N]=[N]=[N]", 123,"123",1.2, "1.2",["CCC"],("CCC"),{"Mol":["CCC"]}]

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        #Variables for ANB
        self.mol1 = None
        self.mol2 = None
        self.mcs_results =  None
        self.index_tuple = None
        self.mol_1 = None
        self.mol_2 = None
        self.mcs_mol_2 = None
        self.frag = None
        self.smallest_frag_list = None
        self.middle_frag_list = None
        self.largest_frag_list = None
        self.CB_frag = None
        self.CB_smallest_frag_list = None
        self.CB_middle_frag_list = None
        self.CB_largest_frag_list = None

        #Variables for MC
        self.Bs_to_Is = None
        self.Is_to_Bs = None
        self.aMapping = None

        #Variables for MWC
        self.bad_SMILES_list = None

    # # Alignment_and_Breaks
    def test_handle_mcs_align_labeling_and_cyclicbreaks(self):
        
        results = ANB.handle_mcs_align_labeling_and_cyclicbreaks(self.mol1, self.mol2, self.mcs_mol)

        self.assertEqual(type(results[0]), type(Chem.MolFromSmiles("C")))
        self.assertEqual(type(results[1]), type(Chem.MolFromSmiles("C")))
        self.assertEqual(type(results[2]), type(Chem.MolFromSmiles("C")))
        
        for atom in results[0].GetAtoms():
            if atom.GetIsotope() < 1000:
                self.assertEqual("A","B")
        for atom in results[1].GetAtoms():
            if atom.GetIsotope() < 2000:
                self.assertEqual("A","B")
    def test_check_cyclic_breaks(self):
        
        results = ANB.check_cyclic_breaks(self.index_tuple, self.mol_1, self.mol_2, self.mcs_mol_2)

        self.assertEqual(type(results[0]), type(Chem.MolFromSmiles("C")))
        
        # Asserts that the mcs has been changed. Should be 1 atom less than prior.
        self.assertNotEqual(results[0].GetNumAtoms(), self.mcs_mol_2.GetNumAtoms())
        self.assertNotEqual(Chem.MolToSmiles(results[0]), Chem.MolToSmiles(self.mcs_mol_2))
        # This is the modified mcs SMILES string
        self.assertEqual(Chem.MolToSmiles(results[0]), '[10000CH2][10002CH2][10003N]=[10004N+]=[10005NH]')

        self.assertEqual(results[1], ([10, 9, 8, 7, 5], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4]))
        self.assertEqual(results[2], True)  
    def test_ringbreak_frag_handling(self):
        """
        This will take a mol and a list of frags and either expand the list to include additional atoms
            from preexisting frags or frags which would be created if we deleted an atom in the original list
        """
        small_and_mid_frag_list = copy.deepcopy(self.middle_frag_list)
        tmp = copy.deepcopy(self.smallest_frag_list)
        small_and_mid_frag_list.extend(tmp)
        small_and_mid_frag_list.sort()
        
        copy_frag = copy.deepcopy(self.frag)

        test_1 = ANB.ringbreak_frag_handling(copy_frag, small_and_mid_frag_list)
        test_1.sort()

        self.assertEqual(test_1, small_and_mid_frag_list)
        self.assertEqual(Chem.MolToSmiles(copy_frag), Chem.MolToSmiles(self.frag))
    def test_ringbreak_frag_handling_2(self):
        """
        This will take a mol and a list of frags and either expand the list to include additional atoms
            from preexisting frags or frags which would be created if we deleted an atom in the original list
        """
        small_and_mid_frag_list = copy.deepcopy(self.middle_frag_list)
        tmp = copy.deepcopy(self.smallest_frag_list)
        small_and_mid_frag_list.extend(tmp)
        small_and_mid_frag_list.sort()
        
        copy_frag = copy.deepcopy(self.frag)
        middle_frag_list = copy.deepcopy(self.middle_frag_list)
             


        test_2 = ANB.ringbreak_frag_handling(copy_frag, middle_frag_list)
        test_2.sort()


        self.assertEqual(test_2, small_and_mid_frag_list)
        self.assertEqual(Chem.MolToSmiles(copy_frag), Chem.MolToSmiles(self.frag))    
    def test_ringbreak_frag_handling_3(self):
        """
        This will take a mol and a list of frags and either expand the list to include additional atoms
            from preexisting frags or frags which would be created if we deleted an atom in the original list
        """
        
        small_and_mid_frag_list = copy.deepcopy(self.middle_frag_list)
        tmp = copy.deepcopy(self.smallest_frag_list)
        small_and_mid_frag_list.extend(tmp)
        small_and_mid_frag_list.sort()
        
        copy_frag = copy.deepcopy(self.frag)
        small_frag_list = copy.deepcopy(self.smallest_frag_list)
             


        test_3 = ANB.ringbreak_frag_handling(copy_frag, [18,19])
        test_3.sort()
        
        self.assertEqual(test_3, small_and_mid_frag_list)
        self.assertEqual(Chem.MolToSmiles(copy_frag), Chem.MolToSmiles(self.frag))
    def test_ringbreak_frag_handling_4(self):
        """
        This will take a mol and a list of frags and either expand the list to include additional atoms
            from preexisting frags or frags which would be created if we deleted an atom in the original list
        """
        # Test if theres a new cyclic break that would need to be fixed
        CB_small_and_mid_frag_list = self.CB_middle_frag_list
        CB_small_and_mid_frag_list.extend(self.CB_smallest_frag_list)
        CB_small_and_mid_frag_list.sort()
        CB_copy_frag = self.CB_frag

        test_4 = ANB.ringbreak_frag_handling(self.CB_frag, CB_small_and_mid_frag_list)
        test_4.sort()
        self.assertNotEqual(test_4, CB_small_and_mid_frag_list)
        self.assertNotEqual(len(test_4), len(CB_small_and_mid_frag_list))
        self.assertEqual(Chem.MolToSmiles(CB_copy_frag), Chem.MolToSmiles(self.CB_frag))
    def test_find_biggest_frag(self):
        """
        This will take a mol and a list of frags and either expand the list to include additional atoms
            from preexisting frags or frags which would be created if we deleted an atom in the original list
        """
        
        frag = Chem.GetMolFrags(self.CB_frag, asMols = True, sanitizeFrags = False)
        biggest_frag, idx_big_frag = ANB.find_biggest_frag(frag)
        self.assertEqual(idx_big_frag, 0)
        self.assertEqual(Chem.MolToSmiles(biggest_frag[0]), '[10000CH3][10001CH2][10002CH2][10003CH2][10004CH2][10005CH2][10006CH2][10007CH2][10008CH2][10009CH2][10010CH2][10011CH2][10012CH2][10013CH]([10014CH2][10015CH2][10016OH])[10017c]([10018CH3])[10019cH2]') 
        frag = Chem.GetMolFrags(self.CB_frag, asMols = True, sanitizeFrags = False)
    
        #test with a mol with no fragments
        mol = Chem.GetMolFrags(self.mol1, asMols = True, sanitizeFrags = False)
        biggest_frag, idx_big_frag = ANB.find_biggest_frag(mol)
        self.assertEqual(idx_big_frag, 0)
        self.assertEqual(Chem.MolToSmiles(biggest_frag[0]), 'CCCCCCCCc1ccccc1CCCCCCCS')     
    def test_remove_iso_labels(self):
        """
        this will remove all of the idx's, the list can be reduced to select atoms
        but for testing the list is the entire idx for the mol

        This doesn't return a value but instead just adjusts the mol
        """
        list_of_idx_to_remove = []
        for atom in self.mol_1.GetAtoms():
            list_of_idx_to_remove.append(atom.GetIdx())

        ANB.remove_iso_labels(self.mol_1,list_of_idx_to_remove)
        self.assertEqual(Chem.MolToSmiles(self.mol_1),"CC1=C(OCC(O)CN=[N+]=[N-])C(C)SC1=O")
    def test_remove_iso_labels_empty(self):
        """
        same as above test but w an empty list
        """
        list_of_idx_to_remove = []
        mol_str_before = Chem.MolToSmiles(self.mol_1)
        ANB.remove_iso_labels(self.mol_1,list_of_idx_to_remove)
        self.assertEqual(Chem.MolToSmiles(self.mol_1),mol_str_before)
    def test_add_r_atom_isolabels(self):
        """
        for any atom which doesn't have an isolabel already add 1000+idx to mol_1 or 2000+idx for mol_2
        

        This doesn't return a value but instead just adjusts the mol
        """    

        # This commented out section fails before add_r_atom_isolabels as it should
        # for atom in self.mol_1.GetAtoms():
        #     if atom.GetIsotope() < 1000:
        #         self.assertEqual("A","B")
        # for atom in self.mol_2.GetAtoms():
        #     if atom.GetIsotope() < 1000:
        #         self.assertEqual("A","B")


        mol_1_str_before = Chem.MolToSmiles(self.mol_1)
        mol_2_str_before = Chem.MolToSmiles(self.mol_2)     
        ANB.add_r_atom_isolabels(self.mol_1, self.mol_2)
        mol_1_str_after = Chem.MolToSmiles(self.mol_1)
        mol_2_str_after = Chem.MolToSmiles(self.mol_2)

        self.assertNotEqual(mol_1_str_before, mol_1_str_after)
        self.assertNotEqual(mol_2_str_before, mol_2_str_after)

        for atom in self.mol_1.GetAtoms():
            if atom.GetIsotope() < 1000:
                self.assertEqual("A","B")
        for atom in self.mol_2.GetAtoms():
            if atom.GetIsotope() < 1000:
                self.assertEqual("A","B")
    def test_pick_mcs_alignment(self):

        self.assertEqual(ANB.pick_mcs_alignment(self.mol_1, self.mol_2, self.mcs_mol_2),self.index_tuple)     
    def test_add_mcs_isolabels(self):

        mol_1_str_before = Chem.MolToSmiles(self.mol_1)
        mol_2_str_before = Chem.MolToSmiles(self.mol_2)     
        mcs_str_before = Chem.MolToSmiles(self.mcs_mol_2)     

        results = ANB.add_mcs_isolabels(self.mol_1, self.mol_2, self.mcs_mol_2,self.index_tuple)
               

        mol_1_str_after = Chem.MolToSmiles(self.mol_1)
        mol_2_str_after = Chem.MolToSmiles(self.mol_2) 
        mcs_str_after = Chem.MolToSmiles(self.mcs_mol_2)     

        # Shows that the mols have been renumbered
        self.assertNotEqual(mol_1_str_before, mol_1_str_after)
        self.assertNotEqual(mol_2_str_before, mol_2_str_after)
        self.assertNotEqual(mcs_str_before, mcs_str_after)

        # prove that all the isotopes are iso labeled consistently
        for lig1, lig2, c1 in zip(self.index_tuple[0],self.index_tuple[1],self.index_tuple[2]):
            atom1 = self.mol_1.GetAtomWithIdx(lig1)
            atom2 = self.mol_2.GetAtomWithIdx(lig2)        
            atomC = self.mcs_mol_2.GetAtomWithIdx(c1)
    
            self.assertEqual(atom1.GetIsotope(),atom2.GetIsotope())
            self.assertEqual(atomC.GetIsotope(),atom1.GetIsotope())
            self.assertEqual(atomC.GetIsotope(),atom2.GetIsotope())
    def make_list_for_test_renumber_to_mcs(self,mol):

        ordered_str_atomNum = ""
        atoms = mol.GetAtoms()
        num_atom = len(atoms)
        for i in range(0,num_atom):
            atom = mol.GetAtomWithIdx(i)
            at_num = atom.GetAtomicNum()
            ordered_str_atomNum = ordered_str_atomNum + str(at_num)
        return ordered_str_atomNum
    def test_renumber_to_mcs(self):

        mol_1_order_str_before = self.make_list_for_test_renumber_to_mcs(self.mol_1)
        mol_2_order_str_before =  self.make_list_for_test_renumber_to_mcs(self.mol_2)     

        mol_1 = ANB.renumber_to_mcs(self.mol_1 ,self.index_tuple[0])
        # Normally this would be using self.index_tuple[1] but this mol is already
        #     numbered to mcs so it doesn't change after ANB.renumber_to_mcs
        #     but it can still be randomly renumbered to a different 
        mol_2 = ANB.renumber_to_mcs(self.mol_2, self.index_tuple[0])
               
        mol_1_order_str_after =  self.make_list_for_test_renumber_to_mcs(mol_1)
        mol_2_order_str_after = self.make_list_for_test_renumber_to_mcs(mol_2)     

        # The before and after changes
        self.assertNotEqual(mol_1_order_str_before, mol_1_order_str_after)
        self.assertNotEqual(mol_2_order_str_before, mol_2_order_str_after)


    # # MappingClass.py tests
    def test_mapping_locate_bs(self):
        aMapping = copy.deepcopy(self.aMapping) 
        
        self.assertEqual(aMapping.locate_b(10000),['1B1'])

        for i in self.Is_to_Bs.keys():
            items = self.Is_to_Bs[i]
             
            self.assertEqual(aMapping.locate_b(i),items)
    def test_mapping_locate_is(self):
        aMapping = copy.deepcopy(self.aMapping) 
        
        self.assertEqual(aMapping.locate_i('1B1'),[10000])
        
        for i in self.Bs_to_Is.keys():
            items = self.Bs_to_Is[i]
            self.assertEqual(aMapping.locate_i(i),items)
    def test_mapping_delete_b(self):
        aMapping = copy.deepcopy(self.aMapping) 
        self.assertEqual(self.aMapping.locate_b(10002), ['2B2','2B3'])
        self.assertEqual(aMapping.locate_b(10002),['2B2','2B3'])

        aMapping.delete_b('2B2')

        self.assertNotEqual(aMapping.locate_b(10002),['2B2','2B3'])
        self.assertNotEqual(aMapping.locate_b(10002), self.aMapping.locate_b(10002))
        self.assertEqual(aMapping.locate_b(10002),['2B3'])
    def test_mapping_delete_i(self):
        aMapping = copy.deepcopy(self.aMapping) 
        self.assertEqual(self.aMapping.locate_i('2B3'), [10003, 10002])
        self.assertEqual(aMapping.locate_i('2B3'), [10003, 10002])

        aMapping.delete_i(10002)

        self.assertNotEqual(aMapping.locate_i('2B3'), [10003, 10002])
        self.assertNotEqual(aMapping.locate_i('2B3'), self.aMapping.locate_i('2B3'))
        self.assertEqual(aMapping.locate_i('2B3'), [10003])
    def test_chose_b_from_i(self):
        aMapping = copy.deepcopy(self.aMapping) 

        original_Bs_to_Is, original_Is_to_Bs = self.aMapping.testing_function_return_self_dicts()
        self.assertEqual(aMapping.locate_i('2B3'), [10003, 10002])
        self.assertEqual(aMapping.locate_b(10002),['2B2','2B3'])

        result = aMapping.chose_b_from_i(10002)
        current_Bs_to_Is, current_Is_to_Bs = aMapping.testing_function_return_self_dicts()
 
 
        self.assertNotEqual(original_Bs_to_Is, current_Bs_to_Is)
        self.assertNotEqual(original_Is_to_Bs, current_Is_to_Bs)     
    def test_run_mapping_class(self):
        
        B_chosen = MC.run_mapping(self.Bs_to_Is, self.Is_to_Bs)

        self.assertEqual(B_chosen, ['1B1', '1B2', '2B2', '2B1'])


    # # Merge_w_core.py
    def test_remove_all_isolabels(self):

        mol = self.mcs_mol_2
        mol_no_iso = MWC.remove_all_isolabels(mol)
        self.assertEqual(Chem.MolToSmiles(mol_no_iso), Chem.MolToSmiles(mol))

        
        for i in self.bad_SMILES_list:
            with self.assertRaises(Exception) as cm:
                MWC.remove_all_isolabels(i)
    def test_make_dict_all_atoms_iso_to_idx_dict(self):


        mol = self.mcs_mol_2
        temp_dic = {}
        for atom in mol.GetAtoms():
            x = atom.GetIdx()
            atom.SetIsotope(x +10000)
            temp_dic[x+10000] = x

        mol_iso_to_Idx_dict = MWC.make_dict_all_atoms_iso_to_idx_dict(mol)
        self.assertEqual(type(mol_iso_to_Idx_dict),dict)
        self.assertEqual(mol_iso_to_Idx_dict,temp_dic)
    def test_make_anchor_to_bonds_and_type_for_frag(self):
        
        """
        In the 1st test case, there is only 1 connecting atom ([10003*]) which is only bound to 
        1 atom ([1015CH]) by a single bond.
        This should result in a dict with a key of 100003 and an item of [1015, rdkit.Chem.rdchem.BondType.SINGLE]

        """
        mol = Chem.MolFromSmiles("[10003*][1015CH]([1016OH])[1013C](=[1014O])[1012O][1011CH2][1010CH3]",sanitize=False)
        anchor_to_connection_dict = MWC.make_anchor_to_bonds_and_type_for_frag(mol)
        known_answer = {10003: [[1015, rdkit.Chem.rdchem.BondType.SINGLE]]}

        self.assertEqual(type(anchor_to_connection_dict),dict)
        self.assertEqual(anchor_to_connection_dict,known_answer)

        mol2 = Chem.MolFromSmiles("[1000CH2]1[1001CH2][1002CH2][10003*]2[1004CH2][1005CH2][1006CH2][1007CH2][1008CH]2[1009CH2]1",sanitize=False)
        anchor_to_connection_dict = MWC.make_anchor_to_bonds_and_type_for_frag(mol2)

        known_answer = {10003: [[1002, rdkit.Chem.rdchem.BondType.SINGLE], [1004, rdkit.Chem.rdchem.BondType.SINGLE], [1008, rdkit.Chem.rdchem.BondType.SINGLE]]}
      
        self.assertEqual(type(anchor_to_connection_dict),dict)
        self.assertEqual(anchor_to_connection_dict,known_answer)
    def test_unpack_lists_of_atoms_and_bond_type(self):

        #Test with 1 connection
        anchor_to_connection_dict  = {10003: [2004, rdkit.Chem.rdchem.BondType.SINGLE]}
        anchor_atom_iso = 10003
        core_merg_iso_to_idx_dict = {10000: 0, 10001: 1, 10002: 2, 10003: 3, 2004: 4, 2005: 5, 2006: 6, 2007: 7, 2008: 8, 2009: 9, 2010: 10, 2011: 11, 2012: 12, 2013: 13}

        list_of_atom_idx, list_of_bond_types = MWC.unpack_lists_of_atoms_and_bond_type(anchor_to_connection_dict, anchor_atom_iso, core_merg_iso_to_idx_dict)

        known_list_of_atom_idx = [4]
        known_list_of_bonds = [rdkit.Chem.rdchem.BondType.SINGLE]
        
        self.assertEqual(list_of_atom_idx,known_list_of_atom_idx)
        self.assertEqual(list_of_bond_types,known_list_of_bonds)

        #Test with 1 connection but in list of list
        anchor_to_connection_dict  = {10003: [[2004, rdkit.Chem.rdchem.BondType.SINGLE]]}
        anchor_atom_iso = 10003
        core_merg_iso_to_idx_dict = {10000: 0, 10001: 1, 10002: 2, 10003: 3, 2004: 4, 2005: 5, 2006: 6, 2007: 7, 2008: 8, 2009: 9, 2010: 10, 2011: 11, 2012: 12, 2013: 13}

        list_of_atom_idx, list_of_bond_types = MWC.unpack_lists_of_atoms_and_bond_type(anchor_to_connection_dict, anchor_atom_iso, core_merg_iso_to_idx_dict)

        known_list_of_atom_idx = [4]
        known_list_of_bonds = [rdkit.Chem.rdchem.BondType.SINGLE]
        
        self.assertEqual(list_of_atom_idx,known_list_of_atom_idx)
        self.assertEqual(list_of_bond_types,known_list_of_bonds)


        #Test with 2 connections
        anchor_to_connection_dict2  ={10003: [[1002,rdkit.Chem.rdchem.BondType.SINGLE],[1004, rdkit.Chem.rdchem.BondType.SINGLE],[1008, rdkit.Chem.rdchem.BondType.SINGLE]]}
        
        anchor_atom_iso2 = 10003
        core_merg_iso_to_idx_dict2 ={1000: 0, 1001: 1, 1002: 2, 1004: 4, 1005: 5, 1006: 6, 1007: 7, 1008: 8, 1009: 9, 10003: 3}
        list_of_atom_idx2, list_of_bond_types2 = MWC.unpack_lists_of_atoms_and_bond_type(anchor_to_connection_dict2, anchor_atom_iso2, core_merg_iso_to_idx_dict2)

        known_list_of_atom_idx2 = [2, 4, 8] 
        known_list_of_bonds2 = [rdkit.Chem.rdchem.BondType.SINGLE, rdkit.Chem.rdchem.BondType.SINGLE, rdkit.Chem.rdchem.BondType.SINGLE]

        self.assertEqual(list_of_atom_idx2, known_list_of_atom_idx2)
        self.assertEqual(list_of_bond_types2, known_list_of_bonds2)
    def test_unpack_lists_of_atoms_and_bond_type(self):
        
        ###########anchor_to_connection_dict, anchor_atom_iso, core_merg_iso_to_idx_dict)
        ###########
        rs_chosen_smiles = [['[10000*][1006CH]([1007CH3])[1008O][1009CH2][1010CH]1[1011CH2][1012O]1']]
        mcs_mol= Chem.MolFromSmiles("[10005NH]=[10004N+]=[10003N][10000SH](=[10001O])=[10002O]")

        rw_core_merg = MWC.merge_smiles_with_core(rs_chosen_smiles, mcs_mol)

        known_mol_sting = "[1007CH3][1006CH]([1008O][1009CH2][1010CH]1[1011CH2][1012O]1)[10000SH](=[10001O])(=[10002O])[10003N]=[10004N+]=[10005NH]"
        self.assertEqual(Chem.MolToSmiles(rw_core_merg),known_mol_sting)
        ###########


    # # Dict_and_r_groups.py  DRG
    def test_get_rs_chosen_smiles(self):

        """returns a list containing the smiles strings for each R group in rs_chosen"""

        rs_chosen=['1R1']
        r_smiles_dict_1={'1R1': '[10003*][1004CH]'}
        r_smiles_dict_2={'2R2': '[10003*][2006CH]([2007OH])[2008CH2][2009OH]', '2R1': '[10003*][2004CH2][2005OH]'}
        known_rs_chosen_smiles=[['[10003*][1004CH]']]

        rs_chosen_smiles= DRG.get_rs_chosen_smiles(rs_chosen,r_smiles_dict_1,r_smiles_dict_2)
        
        self.assertEqual(rs_chosen_smiles,known_rs_chosen_smiles)

        # Test multiple R groups
        rs_chosen2=['1R1','2R1']
        known_rs_chosen_smiles2=[['[10003*][1004CH]'], ['[10003*][2004CH2][2005OH]']]

        rs_chosen_smiles= DRG.get_rs_chosen_smiles(rs_chosen2,r_smiles_dict_1,r_smiles_dict_2)
        self.assertEqual(rs_chosen_smiles,known_rs_chosen_smiles2)
    def test_get_rs_chosen_from_bs(self):


        bs_chosen = ['1B1', '1B2', '2B1']
        b_to_r_master_dict_1 = {'1B1': ['1R1'], '1B2': ['1R2']}
        b_to_r_master_dict_2 = {'2B1': ['2R1']}
        known_rs_chosen = ['1R1', '1R2', '2R1']
        rs_chosen = DRG.get_rs_chosen_from_bs(bs_chosen,b_to_r_master_dict_1, b_to_r_master_dict_2)
        self.assertEqual(rs_chosen,known_rs_chosen)
    def test_replace_core_mol_dummy_atoms(self):
        mol = Chem.MolFromSmiles("[10000N-]=[10001N+]=[10002N][10003CH2][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        mcs = Chem.MolFromSmiles("[10003CH3][10002N]=[10001N+]=[10000NH]")
        Replace_core = Chem.MolFromSmiles("[3*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        knownmol_mod = '[10003*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]'

        mol_mod = DRG.replace_core_mol_dummy_atoms(mol, mcs, Replace_core)
        self.assertEqual(Chem.MolToSmiles(mol_mod),knownmol_mod)

    def test_get_idx_using_unique_iso(self):
             
        # this should be a dic with 8 entries, only 10000 and 100007 should have 1 entry, all else have 2
        smile = "[1000CH3][1001CH2][1002CH2][1003CH2][1004CH2][1005CH2][1006CH2][10007CH3]"
        mol = Chem.MolFromSmiles(smile)
       
        for atom in mol.GetAtoms():
            idx_known = atom.GetIdx()
            iso_val = atom.GetIsotope()

            idx = DRG.get_idx_using_unique_iso(mol, iso_val)

            self.assertEqual(idx_known, idx)
    def test_get_atoms_touch_mcs(self):
        """
        Function to find all neighbors for a set of molecules touching 
        Isolabeled core atoms"""

        
        mol = Chem.MolFromSmiles("[10000N-]=[10001N+]=[10002N][10003CH2][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        mcs_touches = DRG.get_atoms_touch_mcs(mol)

        known_dict = {10003: [4]}
        self.assertEqual(mcs_touches,known_dict)
    def test_invert_dictionary(self):
        

        orig_dict = {1000: [1,2,3,4],1001:[1,4],1002:[5],1003:[3]}
        new_dict = DRG.invert_dictionary(orig_dict)

        known_dict  = {1:[1000,1001],2:[1000],3:[1000,1003],4:[1000,1001],5:[1002]}
        self.assertEqual(new_dict,known_dict)
    def test_make_b_dic(self):
        i_dictionary = {}
        r_dict_num = {}
        lig_number=0

        i_dictionary = {10008:['1R1','1R2'],10009:['1R2','1R3']}
        r_dict_num = {'1R1':[10008],'1R2':[10008,10009],'1R3':[10009]}
        lig_number = 1
        B_to_R_master_dict, B_to_Anchor_master_dict = DRG.make_b_dic(i_dictionary, r_dict_num, lig_number)
        # print("")
        # print("B_to_R_master_dict")
        # print(B_to_R_master_dict)
        # print("")
        # print("B_to_Anchor_master_dict")
        # print(B_to_Anchor_master_dict)
        # print("")
        known_B_to_R_master_dict  = {'1B1': ['1R1', '1R2']}
        known_B_to_Anchor_master_dict  = {'1B1': [10008, 10009]}

        self.assertEqual(B_to_R_master_dict,known_B_to_R_master_dict)
        self.assertEqual(B_to_Anchor_master_dict,known_B_to_Anchor_master_dict)
    def test_get_idx_using_unique_iso(self):

        mol = Chem.MolFromSmiles("CCCCC")
        for atom in mol.GetAtoms():
            atom.SetIsotope(1000+atom.GetIdx())

        for iso in [1000,1001,1002,1003,1004]:
            idx = DRG.get_idx_using_unique_iso(mol, iso)
            self.assertEqual(idx, iso-1000)

        # floats should also work
        for iso in [1000.0,1001.0,1002.0,1003.0,1004.0]:
            idx = DRG.get_idx_using_unique_iso(mol, iso)
            self.assertEqual(idx, iso-1000)

        # Test things that shouldn't work and should return None
        for iso in [None, 100000, "1000", 1000.1]:
            idx = DRG.get_idx_using_unique_iso(mol, iso)
            self.assertEqual(None, DRG.get_idx_using_unique_iso(mol, iso))
    def test_get_r_dict(self):

        r_chain_dict = {'1R1': [3, 4, 5, 6, 7, 8, 9, 10, 11, 10000]}
        lig_r_atom_touch_mcs = {3: [10000]}


        test_Rs_dict = DRG.get_r_dict(r_chain_dict, lig_r_atom_touch_mcs)
        Known_Rs_dict = {'1R1': [10000]}
        
        self.assertEqual(Known_Rs_dict, test_Rs_dict)
    def test_r_groups_dict(self):

        mol = Chem.MolFromSmiles("[10005*]=[1007O].[10003*][1008c]1[1009cH][1010cH][1011c]([1012OH])[1013cH][1014cH]1")
        mol_frags = Chem.GetMolFrags(mol, asMols = True, sanitizeFrags = False)

        lig_number_for_multiplier = 1

        r_chain_dictionary, R_Smiles_dictionary = DRG.r_groups_dict(mol_frags, lig_number_for_multiplier)

        known_r_chain_dictionary = {'1R1': [10005, 7], '1R2': [10003, 8, 9, 10, 11, 12, 13, 14]}
        Known_R_Smiles_dictionary = {'1R1': '[10005*]=[1007O]', '1R2': '[10003*][1008c]1[1009cH][1010cH][1011c]([1012OH])[1013cH][1014cH]1'}

        self.assertEqual(r_chain_dictionary, known_r_chain_dictionary)
        self.assertEqual(R_Smiles_dictionary, Known_R_Smiles_dictionary)
    def test_iso_to_IDX_replace_dictionary(self):

        mol = Chem.MolFromSmiles("[10000N-]=[10001N+]=[10002N][10003CH2][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        mcs = Chem.MolFromSmiles("[10003CH3][10002N]=[10001N+]=[10000NH]")
        Replace_core = Chem.MolFromSmiles("[3*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
                    
        replace_core_mol = DRG.replace_core_mol_dummy_atoms(mol, mcs, Replace_core)


        known_Replace_core = '[10003*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]'

        self.assertEqual(Chem.MolToSmiles(replace_core_mol), known_Replace_core)
    def test_r_group_list(self):

        mol = Chem.MolFromSmiles("[10003N-]=[10002N+]=[10001N][10000CH2][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        mcs = Chem.MolFromSmiles("[10000CH3][10001N]=[10002N+]=[10003NH]")

        Replace_core = DRG.r_group_list(mol, mcs)

        known_Replace_core = "*[2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]"
                            
        self.assertEqual(Chem.MolToSmiles(Replace_core), known_Replace_core)
    def test_mol_handling_of_fragmenting_labeling_and_indexing(self):

        #######
        mol_1 = Chem.MolFromSmiles('[10000C][10001C][10002C][1003C]([1004C])[1005C]')
        mcs = Chem.MolFromSmiles('[10000CH3][10001CH2][10002CH3]')

        known_r_smiles_dict_1 ={'1R1': '[10002*][1003C]([1004C])[1005C]'}
        known_b_to_r_master_dict_1 = {'1B1': ['1R1']}
        known_B_to_Anchor_master_dict_1 = {'1B1': [10002]}
  

        ########
        r_smiles_dict_1, b_to_r_master_dict_1 ,B_to_Anchor_master_dict_1 = DRG.mol_handling_of_fragmenting_labeling_and_indexing(mol_1, mcs, 1)
    

        self.assertEqual(r_smiles_dict_1, known_r_smiles_dict_1)
        self.assertEqual(b_to_r_master_dict_1, known_b_to_r_master_dict_1)
        self.assertEqual(B_to_Anchor_master_dict_1, known_B_to_Anchor_master_dict_1)
    def test_handle_dicts_and_select_b_groups(self):

        #######
        mol_1 = Chem.MolFromSmiles('[10000C][10001C][10002C][1003C]([1004C])[1005C]')
        mol_2 = Chem.MolFromSmiles('[10000C][10001C][10002C][1003O][1004C]=[1005O]')
        mcs = Chem.MolFromSmiles('[10000CH3][10001CH2][10002CH3]')

        rs_chosen_smiles = DRG.handle_dicts_and_select_b_groups(mol_1, mol_2, mcs)
    
        known_rs_chosen_smiles = [['[10002*][1003C]([1004C])[1005C]']]
        self.assertEqual(rs_chosen_smiles, known_rs_chosen_smiles)




    # # Execute_Crossover.py  EC
    
    def test_test_for_mcs(self):


        vars = {}
        vars["max_time_mcs_prescreen"] = 1
        vars["min_atom_match_mcs"] = 3

        # Test normally
        mol1 = Chem.MolFromSmiles("CCC")
        mol2 = Chem.MolFromSmiles("CCCCCCCCC")
        result = EC.test_for_mcs(vars, mol1, mol2)
        self.assertEqual(result.smartsString, "[#6]-[#6]-[#6]")
        self.assertEqual(result.numAtoms, 3)
        self.assertEqual(result.numBonds, 2)


        # test the min_atom_match_mcs function (Set it 1 higher than the actual match count)
        # Result should be None
        vars = {}
        vars["max_time_mcs_prescreen"] = 1
        vars["min_atom_match_mcs"] = 4
        
        # Test normally
        mol1 = Chem.MolFromSmiles("CCC")
        mol2 = Chem.MolFromSmiles("CCCCCCCCC")
        result = EC.test_for_mcs(vars, mol1, mol2)
        self.assertEqual(result, None)

        # test the min_atom_match_mcs function (Set it 1 higher than the actual match count)
        # Result should be None
        vars = {}
        vars["max_time_mcs_prescreen"] = 1
        vars["min_atom_match_mcs"] = 3
        
        # Test with AddH's 
        mol1 = Chem.MolFromSmiles("CCC")
        mol2 = Chem.MolFromSmiles("CCCCCCCCC")
        mol1 = Chem.AddHs(mol1)
        mol2 = Chem.AddHs(mol2)
        result = EC.test_for_mcs(vars, mol1, mol2)
        self.assertEqual(result.smartsString, "[#6](-[#6](-[#6](-[#1])-[#1])(-[#1])-[#1])(-[#1])(-[#1])-[#1]")
        self.assertEqual(result.numAtoms, 10)
        self.assertEqual(result.numBonds, 9)

        # Test Timeout This should time out as its way too big, especially with H's
        # result should be None
        vars = {}
        vars["max_time_mcs_prescreen"] = 0
        vars["min_atom_match_mcs"] = 3
        mol1 = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(CCCCCCCCCCCC)CCCCCCCCCCC(N=[N+]=[NH])CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
        mol2 = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCC(CCCCCCCCCCCC)CCCCCCCCCC(CN=[N+]=[N-])CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")

        result = EC.test_for_mcs(vars, mol1, mol2)

        self.assertEqual(result, None)

        # Test Timeout This should time out as its way too big, especially with H's
        # result should be None
        vars = {}
        vars["max_time_mcs_prescreen"] = 0
        vars["min_atom_match_mcs"] = 3
        mol1 = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(CCCCCCCCCCCC)CCCCCCCCCCC(N=[N+]=[NH])CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
        mol2 = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCC(CCCCCCCCCCCC)CCCCCCCCCC(CN=[N+]=[N-])CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")

        mol1 = Chem.AddHs(mol1)
        mol2 = Chem.AddHs(mol2)
        result = EC.test_for_mcs(vars, mol1, mol2)

        self.assertEqual(result, None)
        # Test Timeout This should time out as its way too big, especially with H's
        # result should be None
        vars = {}
        vars["max_time_mcs_prescreen"] = 1
        vars["min_atom_match_mcs"] = 3
        mol1 = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(CCCCCCCCCCCC)CCCCCCCCCCC(N=[N+]=[NH])CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
        mol2 = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(CCCCCCCCCCCC)CCCCCCCCCC(CN=[N+]=[N-])CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")

        mol1 = Chem.AddHs(mol1)
        mol2 = Chem.AddHs(mol2)
        result = EC.test_for_mcs(vars, mol1, mol2)

        self.assertEqual(result, None)

        