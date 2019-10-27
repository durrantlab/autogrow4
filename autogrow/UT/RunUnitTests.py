#run the unittests

from __future__ import absolute_import
import warnings
import unittest


import autogrow.UT.ScoringTests as ST
import autogrow.UT.MOHTests as MOHT
import autogrow.UT.RankingSelectorTests as RST

import autogrow.UT.CrossoverTests as CT


class UnitTests(object):
    """
    Unit testing object for autogrow.
    """
    def __init__(self):
        """
        Initalizes the unit tests.
        """
        print("init")
        self._suite = unittest.TestSuite()
        self._runner = unittest.TextTestRunner()
        print(self._runner)
    # Running Suite

    def run(self):
        """
        Runs the currently queued suite of tests.
        """
        print("run")
        print("")
        self._runner.run(self._suite)
        print("fin")

    def run_all(self):
        """
        Quickly runs all unit tests.
        """
        print("run all")
        self.add_all_tests()
        self.run()

    # Add specific module tests

    def add_all_tests(self):
        """
        Adds all available tests to the suite.
        """

        
        self.add_CrossoverTests()
        self.add_scoring_tests()
        self.add_MOH_tests()
        self.add_RankingSelectorTests()



    # Crossover Scripts
    def add_CrossoverTests(self):
        """
        Adds the information tests.
        """
        CrossoverTests = unittest.makeSuite(CT.CrossoverTests)
        
        print('CrossoverTests')
        print(CrossoverTests)
        self._suite.addTests(CrossoverTests)


    # These work # out to save time
    def add_RankingSelectorTests(self):
        """
        Adds the information tests.
        """
        RankingSelectorTests = unittest.makeSuite(RST.RankingSelectorTests)
        
        print('RankingSelectorTests')
        print(RankingSelectorTests)
        self._suite.addTests(RankingSelectorTests)

    def add_scoring_tests(self):
        """
        Adds the information tests.
        """
        scoring_tests = unittest.makeSuite(ST.ScoringTests)
        
        print('scoring_tests')
        print(scoring_tests)
        self._suite.addTests(scoring_tests)

    def add_MOH_tests(self):
        """
        Adds the information tests.
        """
        MOH_tests = unittest.makeSuite(MOHT.MOHTests)
        
        print('MOH_tests')
        print(MOH_tests)
        self._suite.addTests(MOH_tests)

if __name__ == '__main__':
    print("begin")
    unittest_obj = UnitTests()
    unittest_obj.run_all()

    print("end")
