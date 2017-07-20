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
                {"filename": "test-one_DAFpop0.obs",
                 "expected_values": (364, 613, 264, 0, 243, 1467, 18, 0, 0, 0, 49, 81, 0, 0, 0, 1293, 0, 0, 0, 0, 0,)},
                {"filename": "test-one_DAFpop1.obs",
                 "expected_values": (423, 419, 1869, 155, 22, 6, 24, 51, 0, 0, 0, 49, 0, 0, 14, 0, 67, 0, 1293, 0, 0,)},
                ]
        for fidx, fixture in enumerate(fixtures):
            field_name_prefix = "stats.testv{}".format(fidx+1)
            data = self.fsc._parse_deme_derived_allele_frequencies(
                    filepath=os.path.join(TESTS_DATA_DIR, "test-one", fixture["filename"]),
                    field_name_prefix=field_name_prefix)
            expected_values = fixture["expected_values"]
            self.assertEqual(len(expected_values), len(data))
            for v1, v2 in zip(expected_values, data.values()):
                self.assertEqual(v1, v2)

if __name__ == "__main__":
    unittest.main()

