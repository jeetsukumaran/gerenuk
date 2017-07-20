#! /usr/bin/env python

import os
import unittest
import time
from collections import Counter
from gerenuk import simulate
from gerenuk.test import TESTS_DATA_DIR

class Fsc2SiteFilepathTestCase(unittest.TestCase):

    def setUp(self):
        self.fsc = simulate.Fsc2Handler(
                name="test-data-1",
                fsc2_path="fsc25",
                working_directory=TESTS_DATA_DIR,
                )
    def test_parameter_filepath(self):
        self.assertEqual(self.fsc.parameter_filepath, "test-data-1.par")

    def test_results_dirpath(self):
        self.assertEqual(self.fsc.results_dirpath, "test-data-1")

    def test_deme0_derived_alllele_frequency_filepath(self):
        self.assertEqual(self.fsc.deme0_derived_alllele_frequency_filepath,
                os.path.join("test-data-1", "test-data-1_DAFpop0.obs"))

    def test_deme0_derived_alllele_frequency_filepath(self):
        self.assertEqual(self.fsc.deme0_derived_alllele_frequency_filepath,
                os.path.join("test-data-1", "test-data-1_DAFpop0.obs"))

    def test_deme1_derived_alllele_frequency_filepath(self):
        self.assertEqual(self.fsc.deme1_derived_alllele_frequency_filepath,
                os.path.join("test-data-1", "test-data-1_DAFpop1.obs"))

    def test_joint_derived_alllele_frequency_filepath(self):
        self.assertEqual(self.fsc.joint_derived_alllele_frequency_filepath,
                os.path.join("test-data-1", "test-data-1_jointDAFpop1_0.obs"))

if __name__ == "__main__":
    unittest.main()

