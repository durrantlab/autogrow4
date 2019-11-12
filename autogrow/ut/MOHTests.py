from __future__ import absolute_import
import unittest
import os
import glob

import rdkit
from rdkit import Chem

import autogrow.operators.ConvertFiles.gypsum_dl.gypsum_dl.MolObjectHandling as MOH

def fun(x):
    y = x + 1
    return y


def funq(x):
    return x + 1

class MOHTests(unittest.TestCase):
        # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        # bad_SMILES_list are the wrong data type for the MOH functions. 
        # MOH functions mostly take rdkit.Chem.rdchem.Mol instead of strings
        # None of these strings can be sanitize in rdkit
        # c1ccccc1(C)(C) is invalid because it isnt viable to be an aromatic ring with a dimethlys on the same carbon.
        self.bad_SMILES_list = ["None","c1ccccc1(C)(C)", "CCC[N]=[N]=[N]", 123,"123",1.2, "1.2",["CCC"],("CCC"),{"Mol":["CCC"]}]

        # This is a list of strings which are capable of being sanitized
        self.valid_SMILES_list = ['CC', 'CCC', 'CCCCC', 'CCCCCC', 'CCCCCCC', 'CC', 'CCCCOC',"CCC[N]=[N-]=[N+]"]
        
        self.Ampicillin = Chem.MolFromSmiles("CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C")
        self.bad_azide = Chem.MolFromSmiles("CCC[N]=[N]=[N]", sanitize=False)
        self.idx_to_remove = [20,22,21,13, 14, 15, 16, 17, 18, 11]

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.bad_SMILES_list = None
        self.valid_SMILES_list = None
        self.Ampicillin = None
        self.idx_to_remove = None
        self.bad_azide = None

    def test_check_sanitization_raise_except_w_wrong_data_types(self):
        """ tests bad mols in the MOH.check_sanitization function
        input for MOH.check_sanitization should be an rdkit.mol.obj 
        but we are giving it the wrong data types...
        This will test the handleH function. These use the wrong data types and thus fails.
        Nones should return None while all other wrong data types will trigger a raise Exception
        """
        # This should raise an Error
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.check_sanitization(i)
        
        # None shouldn't raise an error, but rather return None
        self.assertEqual(MOH.check_sanitization(None), None)
        
        # # Good SMILES but in string format rather than rdkit Mol obj
        for i in self.valid_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.check_sanitization(i)

    def test_check_sanitization_w_valid_mols(self):
        """ tests valid mols in the MOH.check_sanitization function"""
        for i in self.valid_SMILES_list: 
            test_mol = Chem.MolFromSmiles(i, sanitize=False)
            self.assertEqual(type(MOH.check_sanitization(test_mol)), rdkit.Chem.rdchem.Mol)
        
      
    def test_handleH_raise_except_w_wrong_data_types(self):
        """
        This will test the handleH function. These use the wrong data types and thus fails.
        Nones should return None while all other wrong data types will trigger a raise Exception

        handleH requires an rdkit.mol.obj and a bol.
        """
        # This should raise an Error
        #   -with protanation on
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.handleHs(i, True)
        # This should raise an Error
        #   -with protanation off
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.handleHs(i, False)
        
        # None shouldn't raise an error, but rather return None
        self.assertEqual(MOH.handleHs(None, True), None)
        # None shouldn't raise an error, but rather return None
        self.assertEqual(MOH.handleHs(None, False), None)
        
        # # Good SMILES but in string format rather than rdkit Mol obj
        for i in self.valid_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.handleHs(i, False)


    def test_handleH_w_valid_mols(self):
        """
        handleH requires an rdkit.mol.obj and a bol.
        """
        
        for i in self.valid_SMILES_list: 
            mol = Chem.MolFromSmiles(i, sanitize=False)
            result = MOH.handleHs(mol, True)
            self.assertEqual(type(result), rdkit.Chem.rdchem.Mol)


    def test_try_deprotanation_raise_except_w_wrong_data_types(self):
        """Test deprotanation which will result in a raise exception """

        # This should raise an Error
        #   -with protanation on
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.try_deprotanation(i)
        # This should raise an Error
        #   -with protanation off
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.try_deprotanation(i)
        
        # None shouldn't raise an error, but rather return None
        self.assertEqual(MOH.try_deprotanation(None), None)
        # None shouldn't raise an error, but rather return None
        self.assertEqual(MOH.try_deprotanation(None), None)
        
        # # Good SMILES but in string format rather than rdkit Mol obj
        for i in self.valid_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.try_deprotanation(i)

    def test_remove_atom(self):
        mol = self.Ampicillin
        for atom in mol.GetAtoms():
            atom.SetIsotope(atom.GetIdx())
        
        # This tests the removal process and the sorting of the list
        self.assertEqual(MOH.remove_atoms(mol, self.idx_to_remove).GetNumAtoms(), 14)

    def test_remove_atom_two(self):
        mol = self.Ampicillin
        for atom in mol.GetAtoms():
            atom.SetIsotope(atom.GetIdx())
        mol_remove_list = []
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                mol_remove_list.append(atom.GetIsotope())

        self.assertEqual(Chem.MolToSmiles(MOH.remove_atoms(mol, mol_remove_list)),'C[1C]1([23CH3])[5S][4CH]2[6CH]([9NH][10C](=[11O])[12CH2][19NH2])[7C](=[8O])[3N]2[2CH]1[20C](=[21O])[22OH]')
        
    def test_N_fix(self):

        # This should raise an Error
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.adjust_Nitrogen_charges(i)
        

        for i in self.valid_SMILES_list: 
            mol = MOH.adjust_Nitrogen_charges(Chem.MolFromSmiles(i,sanitize=False))

            self.assertEqual(type(Chem.MolToSmiles(mol)), str)

        self.assertEqual(Chem.MolToSmiles(MOH.adjust_Nitrogen_charges(self.bad_azide)), "CCCN=[N+]=N")

        em1 = Chem.EditableMol(Chem.MolFromSmiles("CC"))
        
        # This should raise an Error
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.adjust_Nitrogen_charges(em1)
        
        # give it valid sanitized CCCC have it return a valid sanitize CCC
        #   -assertEqual... apply to most MOH functions

    def test_check_for_unassigned_atom(self):
        
        
        # This should raise an Error
        em1 = Chem.EditableMol(Chem.MolFromSmiles("CC"))
        with self.assertRaises(Exception) as cm:
            MOH.check_for_unassigned_atom(em1)

        mol_fail = Chem.MolFromSmiles("CCCC[*]")
        self.assertEqual(MOH.check_for_unassigned_atom(mol_fail), None)

        mol_pass = Chem.MolFromSmiles("CCCC[C]")
        self.assertNotEqual(MOH.check_for_unassigned_atom(mol_pass), None)

    def test_handle_frag_check(self):
        # This should raise an Error
        for i in self.bad_SMILES_list: 
            with self.assertRaises(Exception) as cm:
                MOH.handle_frag_check(em1)
        em1 = Chem.EditableMol(Chem.MolFromSmiles("CC"))
        with self.assertRaises(Exception) as cm:
            MOH.handle_frag_check(em1)
    
        mol_fail = Chem.MolFromSmiles("CCCC[*]")
        self.assertEqual(MOH.handle_frag_check(mol_fail), mol_fail)

        mol_pass = Chem.MolFromSmiles("CCCC[C]")
        self.assertEqual(MOH.handle_frag_check(mol_pass), mol_pass)

        mol_frag = Chem.MolFromSmiles("CCC.C[C]")
        self.assertEqual(Chem.MolToSmiles(MOH.handle_frag_check(mol_frag)), "CCC")

    