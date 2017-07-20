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
                name="test-one",
                fsc2_path="fsc25",
                working_directory=TESTS_DATA_DIR,
                )

    def test_parameter_filepath(self):
        self.assertEqual(self.fsc.parameter_filepath, "test-one.par")

    def test_results_dirpath(self):
        self.assertEqual(self.fsc.results_dirpath, "test-one")

    def test_deme0_derived_allele_frequency_filepath(self):
        self.assertEqual(self.fsc.deme0_derived_allele_frequency_filepath,
                os.path.join("test-one", "test-one_DAFpop0.obs"))

    def test_deme0_derived_allele_frequency_filepath(self):
        self.assertEqual(self.fsc.deme0_derived_allele_frequency_filepath,
                os.path.join("test-one", "test-one_DAFpop0.obs"))

    def test_deme1_derived_allele_frequency_filepath(self):
        self.assertEqual(self.fsc.deme1_derived_allele_frequency_filepath,
                os.path.join("test-one", "test-one_DAFpop1.obs"))

    def test_joint_derived_allele_frequency_filepath(self):
        self.assertEqual(self.fsc.joint_derived_allele_frequency_filepath,
                os.path.join("test-one", "test-one_jointDAFpop1_0.obs"))

class Fsc2DataExtractionTestCase(unittest.TestCase):

    def setUp(self):
        self.fsc = simulate.Fsc2Handler(
                name="test-one",
                fsc2_path="fsc25",
                working_directory=TESTS_DATA_DIR,
                )

    def test_deme_derived_allele_frequencies(self):
        fixtures = [
                ("test-one_DAFpop0.obs",
                    (364, 613, 264, 0, 243, 1467, 18, 0, 0, 0, 49, 81, 0, 0, 0, 1293, 0, 0, 0, 0, 0,)),
                ("test-one_DAFpop1.obs",
                    (423, 419, 1869, 155, 22, 6, 24, 51, 0, 0, 0, 49, 0, 0, 14, 0, 67, 0, 1293, 0, 0,)),
                ]
        for fidx, fixture in enumerate(fixtures):
            field_name_prefix = "testv{}".format(fidx+1)
            self.fsc._parse_deme_derived_allele_frequencies(
                    filepath=os.path.join(TESTS_DATA_DIR, "test-one", fixture[0]),
                    field_name_prefix="{}")

if __name__ == "__main__":
    unittest.main()

