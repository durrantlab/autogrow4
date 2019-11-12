from __future__ import absolute_import
import unittest
import os
import glob
# import sys
# script_dir = str(os.path.dirname(os.path.realpath(__file__))).replace('unittest',"")
# scoring_folder = script_dir + 'autogrow/docking/scoring/'

# ES = __import__(scoring_folder+'Execute_Scoring.py')
from autogrow.docking.scoring import Execute_Scoring  as ES
from autogrow.docking.scoring.scoring_classes.scoring_functions import VINA as VINA
from autogrow.docking.scoring.scoring_classes.scoring_functions import LigEfficiency  as LigEfficiency
from autogrow.docking.scoring.scoring_classes.scoring_functions import NN1 as NN1
from autogrow.docking.scoring.scoring_classes.scoring_functions import NN2 as NN2

def fun(x):
    y = x + 1
    return y


def funq(x):
    return x + 1

class ScoringTests(unittest.TestCase):
        # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        self.num = 3
        self.other_num = 4
        self.str = '4'
        self.vars = {'filename_of_receptor': "/home/jacob/Documents/autogrow4/cdc73/3v46.pdb", 'dock_choice': 'QuickVina2Docking', \
        'nn1_script': '/home/jacob/Documents/autogrow4/autogrow/docking/scoring/nn_score_exe/nnscore1/NNScore.py',\
        'docking_executable': '/home/jacob/Documents/autogrow4/autogrow/docking/docking_executables/q_vina_02/q_vina_02', 'size_y': 30.0, 'size_z': 30.0, 'rxn_library': 'ClickChem', \
        'rxn_library_file': '', 'output_directory': '/home/jacob/Desktop/Outputfolder/Run_1/', \
        'num_generations': 2, 'top_mols_to_seed_next_generation': 1, 'scoring_function': 'VINA', \
        'root_output_folder': '/home/jacob/Desktop/Outputfolder/', 'diversity_seed_depreciation_per_gen': 1, \
        'size_x': 30.0, \
        'mgl_python': '/home/jacob/MGLTools-1.5.6/bin/pythonsh',\
        'nn2_script': '/home/jacob/Documents/autogrow4/autogrow/docking/scoring/nn_score_exe/nnscore2/NNScore2.py',\
        'custom_scoring_script': '', 'redock_elite_from_previous_gen': False, \
        'selector_choice': 'Roulette_Selector', 'number_of_processors': -1}

        #self.smile_dict = {'Gen_0_Mutant_1_696701': ['CC(=O)c1cc(OC(C)(C(N)=O)C2(O)CCCCC2)c([O-])c([N+](=O)[O-])c1', '(ZINC00002778+ZINC05004233)Gen_0_Mutant_1_696701'], 'Gen_0_Cross_472868': ['O=C(O)C1C=CNNC1C(=O)[O-]', '(ZINC60190052+ZINC35881877)Gen_0_Cross_472868'], 'ZINC04727099': ['COC1OC(CO)C(N=[N+]=[N-])C(O)C1O', 'ZINC04727099']}
        self.smile_dict = {'Gen_0_Mutant_5_203493': ['COC1OC(CO)C(O)C(O)C1n1nnc(CCO)c1-c1ccc(-c2cccs2)s1', '(ZINC04530731+ZINC01529972)Gen_0_Mutant_5_203493'], 'Gen_0_Cross_452996': ['CC(=O)OCC(O)CN=[N+]=[N-]', '(ZINC44117885+ZINC34601304)Gen_0_Cross_452996'], 'ZINC13526729': ['[N-]=[N+]=NCC1OC(O)CC1O', 'ZINC13526729']}

        self.pdb_folder = '/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/'

        self.pdbqtvina_file = '/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.vina'

        self.pdbqtvina_file_fail = '/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1Fail.pdbqt.vina'
        self.Gen_0_Cross_452996_smile = 'CC(=O)OCC(O)CN=[N+]=[N-]'

        self.nn1_file = "/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.vina.nn1"

        self.nn2_file = "/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.nn2"

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.num = None
        self.other_num = None
        self.str = None
        self.vars = None
        self.smile_dict = None
        self.pdbqtvina_file = None

    def test_fun(self):
        
        self.assertEqual(fun(self.num),self.other_num)
        print("")
        self.assertEqual(fun(self.num), 4)

    def test(self):
        self.assertEqual(funq(self.num),4)

    # Test Vina scoring function
    def test_VINA_class(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        vina_obj = VINA.VINA(self.vars, self.smile_dict)
        
        # 2nd item in list should be -5.4 from demo
        self.assertEqual(vina_obj.run_scoring(self.pdbqtvina_file)[-1], "-5.4")

        # 1st item in list should be Gen_0_Cross_452996 from demo
        self.assertEqual(vina_obj.run_scoring(self.pdbqtvina_file)[2], "Gen_0_Cross_452996")
        
        # return type should be list
        self.assertEqual(type(vina_obj.run_scoring(self.pdbqtvina_file)), type([]))

    def test_VINA_class_fails(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        vina_obj = VINA.VINA(self.vars, self.smile_dict)

        # return type should be list
        self.assertEqual(vina_obj.run_scoring(self.pdbqtvina_file_fail), None)
        
    def test_get_all_VINA_files(self):
        
        all_pdbqt_vina_files = glob.glob(self.pdb_folder + "*.pdbqt.vina")

        vina_obj = VINA.VINA(self.vars, self.smile_dict)

        # return type should be list
        self.assertEqual(vina_obj.find_files_to_score(self.pdb_folder), all_pdbqt_vina_files)
        

    # Test VinaLigEFF scoring function
    def test_VINA_lig_EFF_outside_class_functions(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        vina_obj = VINA_lig_Eff.LigEfficiency(self.vars, self.smile_dict)


        # The smile string is 'CC(=O)OCC(O)CN=[N+]=[N-]' which has 3N 3O and 5C 
        # this is a total of 11 non-Hydrogens
        self.assertEqual(VINA_lig_Eff.get_number_heavy_atoms(self.Gen_0_Cross_452996_smile),11)


        # Same SMILES with Hydrogens explicit [CH3][C](=[O])[O][CH2][CH]([OH])[CH2][N]=[N+]=[N-]
        # this is a total of 11 non-Hydrogens
        self.assertEqual(VINA_lig_Eff.get_number_heavy_atoms('[CH3][C](=[O])[O][CH2][CH]([OH])[CH2][N]=[N+]=[N-]'),11)
        


        # A viable but Bad Smile string is 'CC(=O)OCC(O)CN=[N]=[N]' which has 3N 3O and 5C 
        # this is a total of 11 non-Hydrogens
        self.assertEqual(VINA_lig_Eff.get_number_heavy_atoms('CC(=O)OCC(O)CN=[N]=[N]'),11)


        # A NON Viable Smile string is 123456 Should return None
        self.assertEqual(VINA_lig_Eff.get_number_heavy_atoms(123455), None)

        # A NON Viable Smile string is None Should return None
        self.assertEqual(VINA_lig_Eff.get_number_heavy_atoms(None), None)

        
        self.assertEqual(VINA_lig_Eff.append_lig_effeciency(['CC(=O)OCC(O)CN=[N+]=[N-]', '(ZINC44117885+ZINC34601304)Gen_0_Cross_452996',str(-5.4)]), ['CC(=O)OCC(O)CN=[N+]=[N-]', '(ZINC44117885+ZINC34601304)Gen_0_Cross_452996', str(-5.4), str(-0.490909090909)])


    def test_VINA_lig_EFF_class(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        vina_obj = VINA_lig_Eff.LigEfficiency(self.vars, self.smile_dict)
        
        # 2nd item in list should be -0.490909090909 from demo
        self.assertEqual(float(vina_obj.run_scoring(self.pdbqtvina_file)[-1]), float(-0.490909090909))

        # 1st item in list should be Gen_0_Cross_452996 from demo
        self.assertEqual(vina_obj.run_scoring(self.pdbqtvina_file)[2], "Gen_0_Cross_452996")
        
        # return type should be list
        self.assertEqual(type(vina_obj.run_scoring(self.pdbqtvina_file)), type([]))

    def test_VINA_lig_EFF_class_fails(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        vina_obj = VINA_lig_Eff.LigEfficiency(self.vars, self.smile_dict)

        # return type should be list
        self.assertEqual(vina_obj.run_scoring(self.pdbqtvina_file_fail), None)
        
    def test_get_all_VINA_lig_EFF_files(self):
        
        all_pdbqt_vina_files = glob.glob(self.pdb_folder + "*.pdbqt.vina")

        vina_obj = VINA_lig_Eff.LigEfficiency(self.vars, self.smile_dict)



        # return type should be list
        self.assertEqual(vina_obj.find_files_to_score(self.pdb_folder), all_pdbqt_vina_files)
        

    ### TEST NN1
    def test_NN1_outside_class_functions(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        NN1_obj = NN1.NN1(self.vars, self.smile_dict)

        # The smile string is 'CC(=O)OCC(O)CN=[N+]=[N-]' which has 3N 3O and 5C 
        # this is a total of 11 non-Hydrogens
        self.assertEqual(NN1.run_nn_rescoring(self.vars,['/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.vina']),['/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.vina.nn1'])

        #
      
    ### NN1
    def test_NN1_class(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        NN1_obj = NN1.NN1(self.vars, self.smile_dict)
        

        # -0.5288992907179128
        
        # 2nd item in list should be -0.490909090909 from demo
        self.assertEqual(float(NN1_obj.run_scoring(self.nn1_file)[-1]), float(0.528899290718))

        # return type should be list
        self.assertEqual(type(NN1_obj.run_scoring(self.nn1_file)), type([]))

    def test_NN1_class_fails(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        NN1_obj = NN1.NN1(self.vars, self.smile_dict)

        # return type should be list
        self.assertEqual(NN1_obj.run_scoring(self.pdbqtvina_file_fail), None)
        
    def test_get_all_NN1_files(self):
        
        all_pdbqt_vina_files = glob.glob(self.pdb_folder + "*.pdbqt.vina.nn1")

        NN1_obj = NN1.NN1(self.vars, self.smile_dict)



        # return type should be list
        self.assertEqual(NN1_obj.find_files_to_score(self.pdb_folder), all_pdbqt_vina_files)
        
    ### NN2
    def test_NN2_class(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        NN2_obj = NN2.NN2(self.vars, self.smile_dict)

        self.assertEqual(NN2.run_nn_rescoring(self.vars, ['/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt']),['/home/jacob/Desktop/UnittestExample/vina_tests/Run_5/generation_0/PDBs/Gen_0_Cross_452996__1.pdbqt.nn2'])

        
        # 2nd item in list should be -0.490909090909 from demo
        self.assertEqual(float(NN2_obj.run_scoring(self.nn2_file)[-1]), float(-0.481))

        # return type should be list
        self.assertEqual(type(NN2_obj.run_scoring(self.nn2_file)), type([]))

    def test_NN2_class_fails(self):
        """
        Make a vina class object and retrieve a single score for a file

        The return object from run_scoring is supposed to be a list containing 
        short smile name and the docking score
        """
        NN2_obj = NN2.NN2(self.vars, self.smile_dict)
