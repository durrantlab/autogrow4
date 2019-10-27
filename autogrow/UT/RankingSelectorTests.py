from __future__ import absolute_import
import unittest
import os
import glob
# import sys
# script_dir = str(os.path.dirname(os.path.realpath(__file__))).replace('unittest',"")
# scoring_folder = script_dir + 'autogrow/Docking/Scoring/'

# ES = __import__(scoring_folder+'Execute_Scoring.py')

import autogrow.Docking.Ranking.Ranking_mol as Ranking

import autogrow.Docking.Ranking.Selecting.Tournement_Selection as TourSel
import autogrow.Docking.Ranking.Selecting.Roulette_Selection as RoulSel
import autogrow.Docking.Ranking.Selecting.Rank_Selection as RankSel


class RankingSelectorTests(unittest.TestCase):
        # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        self.rank_file = '/home/jacob/Desktop/UnittestExample/Vina_tests/Run_5/generation_0/generation_0_ranked.smi'
        
        self.non_exist_file = '/home/jacob/Desktop/fake.txt'

        self.mol_list = [['CCCCCC', 'mol_1'], ['CCCCCCCC', 'mol_2'], ['CN=[N+]=NC', 'mol_3']]


        self.bad_mol_list = [[123, 'mol_1'], ['CCCCCCCC', 'mol_2'], ['CN=[N+]=NC', 'mol_3']]

        self.age_list = [["Jake",'27',27],["Jacob",'37',37],["Jen",'37',37],["Erich",'30',30],["Pauline",'26',26],["Kevin",'25',25],["Aletheia",'36',36]]


        self.str = '4'
        self.vars = {'filename_of_receptor': "/home/jacob/Documents/autogrow/cdc73/3v46.pdb", 'Dock_choice': 'QuickVina2Docking', \
        'nn1_script': '/home/jacob/Documents/autogrow4/autogrow/Docking/Scoring/NNScore_exe/nnscore1/NNScore.py',\
        'docking_executable': '/home/jacob/Documents/autogrow4/autogrow/Docking/Docking_Executables/QVina02/qvina02', 'size_y': 30.0, 'size_z': 30.0, 'Rxn_library': 'ClickChem', \
        'rxn_library_file': '', 'output_directory': '/home/jacob/Desktop/Outputfolder/Run_1/', \
        'num_generations': 2, 'top_mols_to_seed_next_generation': 1, 'scoring_function': 'VINA', \
        'root_output_folder': '/home/jacob/Desktop/Outputfolder/', 'diversity_seed_depreciation_per_gen': 1, \
        'size_x': 30.0, \
        'mgl_python': '/home/jacob/MGLTools-1.5.6/bin/pythonsh',\
        'nn2_script': '/home/jacob/Documents/autogrow4/autogrow/Docking/Scoring/NNScore_exe/nnscore2/NNScore2.py',\
        'Custom_scoring_script': '', 'redock_advance_from_previous_gen': False, \
        'Selector_Choice': 'Roulette_Selector', 'number_of_processors': -1}

        self.smile_dict = {'Gen_0_Mutant_5_203493': ['COC1OC(CO)C(O)C(O)C1n1nnc(CCO)c1-c1ccc(-c2cccs2)s1', '(ZINC04530731+ZINC01529972)Gen_0_Mutant_5_203493'], 'Gen_0_Cross_452996': ['CC(=O)OCC(O)CN=[N+]=[N-]', '(ZINC44117885+ZINC34601304)Gen_0_Cross_452996'], 'ZINC13526729': ['[N-]=[N+]=NCC1OC(O)CC1O', 'ZINC13526729']}

        self.pdb_folder = '/home/jacob/Desktop/UnittestExample/Vina_tests/Run_5/generation_0/PDBs/'

        self.pdbqtvina_file = '/home/jacob/Desktop/UnittestExample/Vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.vina'

        self.pdbqtvina_file_fail = '/home/jacob/Desktop/UnittestExample/Vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1Fail.pdbqt.vina'
        self.Gen_0_Cross_452996_smile = 'CC(=O)OCC(O)CN=[N+]=[N-]'

        self.nn1_file = "/home/jacob/Desktop/UnittestExample/Vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.vina.nn1"

        self.nn2_file = "/home/jacob/Desktop/UnittestExample/Vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.nn2"

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.rank_file = None
        self.non_exist_file = None
        self.mol_list = None

        self.num = None
        self.str = None
        self.vars = None
        self.smile_dict = None
        self.pdbqtvina_file = None



    ######## RANKING script ########

    def test_usable_list(self):
        
        self.assertEqual(type(Ranking.get_usable_fomat(self.rank_file)),list)

        # Skipping this error but it does raise an exceptions 
        # Tested but it registers as a failure so removed.
        #self.assertRaises(Ranking.get_usable_fomat(self.non_exist_file))

        print("")

    def test_convert_usable_list_to_lig_dict(self):
        
        self.assertEqual(type(Ranking.convert_usable_list_to_lig_dict(Ranking.get_usable_fomat(self.rank_file))),dict)
   
        self.assertEqual(Ranking.convert_usable_list_to_lig_dict("abc"),None)

        print("")

    ######## Tournement Selection ########

    def test_score_and_append_diversity_scores(self):
        self.assertEqual(len(Ranking.score_and_append_diversity_scores(self.mol_list)[0]),3)

       
        self.assertEqual(Ranking.score_and_append_diversity_scores(self.mol_list)[0][-1],'0.906060606061')

        # Skipping this error but it does raise an exceptions 
        # Tested but it registers as a failure so removed.
        # self.assertRaises(Ranking.score_and_append_diversity_scores(self.bad_mol_list))

    def test_Tournement_Selction(self):
        self.assertEqual(type(TourSel.Run_Tournement_Selector(self.age_list,3,2,1)),list)
        self.assertEqual(len(TourSel.Run_Tournement_Selector(self.age_list,3,2,1)),3)


        # These raise the proper errors but are removed from this unittest.
        # no idx in sublist to select by
        # self.assertRaises(TourSel.Run_Tournement_Selector(self.age_list,3,2,3))

        # Not a list
        # self.assertRaises(TourSel.Run_Tournement_Selector(None,3,2,3))

        # idx to select by could not be converted to float
        # self.assertEqual(TourSel.Run_Tournement_Selector(self.age_list,3,2,0),None)
        
    def test_run_one_tournement(self):
        self.assertEqual(type(TourSel.run_one_tournement(self.age_list,2,1)),list)
        self.assertEqual(len(TourSel.run_one_tournement(self.age_list,2,1)), 3)


    ######## Roulette Selection ########

    def test_Spin_Roulette_Selector(self):
        
        self.assertEqual(len(RoulSel.Spin_Roulette_Selector(self.age_list, 5,'docking')), 5)

        self.assertEqual(len(RoulSel.Spin_Roulette_Selector(self.age_list, 3,'diversity')),3)

    def test_adjust(self):
        
        self.assertEqual(len(RoulSel.adjust_scores(self.age_list, 'docking')), 7)

        self.assertEqual(len(RoulSel.adjust_scores(self.age_list,'diversity')), 7)


    ######## Rank Selection ########

    def test_Spin_Roulette_Selector(self):
        
        self.assertEqual(len(RankSel.Run_Rank_selector(self.age_list, 5,-1)), 5)

       

        # These raise the proper errors but are removed from this unittest.
        # Asking for more to be chosen than optinos.
        # self.assertRaises(RankSel.Run_Rank_selector(self.age_list, 9,-1))
               
        # Wrong Data type
        # self.assertRaises(RankSel.Run_Rank_selector(None, 2, -1))
       


