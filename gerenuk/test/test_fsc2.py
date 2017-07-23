#! /usr/bin/env python

import os
import unittest
import time
import collections
from gerenuk import simulate
from gerenuk.test import TESTS_DATA_DIR

class Fsc2SiteFilepathTestCase(unittest.TestCase):

    def get_fsc_handler(self, is_folded_site_frequency_spectrum):
        fsc_handler = simulate.Fsc2Handler(
                name="test-one",
                fsc2_path="fsc25",
                working_directory=TESTS_DATA_DIR,
                is_folded_site_frequency_spectrum=is_folded_site_frequency_spectrum,
                )
        return fsc_handler

    def test_parameter_filepath(self):
        fsc_handler = self.get_fsc_handler(True)
        self.assertEqual(fsc_handler.parameter_filepath, "test-one.par")

    def test_results_dirpath(self):
        fsc_handler = self.get_fsc_handler(True)
        self.assertEqual(fsc_handler.results_dirpath,
                os.path.join(TESTS_DATA_DIR, "test-one"))

    def test_deme0_derived_allele_frequency_filepath(self):
        for folded, expected in (
                (True, "M"),
                (False, "D"),
                ):
            fsc_handler = self.get_fsc_handler(folded)
            self.assertEqual(fsc_handler.deme0_derived_allele_frequency_filepath,
                    os.path.join(TESTS_DATA_DIR, "test-one", "test-one_{}AFpop0.obs".format(expected)))

    def test_deme0_derived_allele_frequency_filepath(self):
        for folded, expected in (
                (True, "M"),
                (False, "D"),
                ):
            fsc_handler = self.get_fsc_handler(folded)
            self.assertEqual(fsc_handler.deme0_derived_allele_frequency_filepath,
                    os.path.join(TESTS_DATA_DIR, "test-one", "test-one_{}AFpop0.obs".format(expected)))

    def test_deme1_derived_allele_frequency_filepath(self):
        for folded, expected in (
                (True, "M"),
                (False, "D"),
                ):
            fsc_handler = self.get_fsc_handler(folded)
            self.assertEqual(fsc_handler.deme1_derived_allele_frequency_filepath,
                    os.path.join(TESTS_DATA_DIR, "test-one", "test-one_{}AFpop1.obs".format(expected)))

    def test_joint_derived_allele_frequency_filepath(self):
        for folded, expected in (
                (True, "M"),
                (False, "D"),
                ):
            fsc_handler = self.get_fsc_handler(folded)
            self.assertEqual(fsc_handler.joint_derived_allele_frequency_filepath,
                    os.path.join(TESTS_DATA_DIR, "test-one", "test-one_joint{}AFpop1_0.obs".format(expected)))

class Fsc2DataExtractionTestCase(unittest.TestCase):

    def setUp(self):
        self.fsc = simulate.Fsc2Handler(
                name="test-one",
                fsc2_path="fsc25",
                working_directory=TESTS_DATA_DIR,
                is_folded_site_frequency_spectrum=True,
                )

    def test_deme_derived_allele_frequencies(self):
        fixtures = [
                {"filename": "test-one_MAFpop0.obs",
                 "expected_values": (929, 110, 238, 101, 0, 216)},
                {"filename": "test-one_MAFpop1.obs",
                 "expected_values": (39, 1100, 98, 40, 0, 101, 214, 2, 0)},
                ]
        for fidx, fixture in enumerate(fixtures):
            field_name_prefix = "stats.testv{}".format(fidx+1)
            data = collections.OrderedDict()
            self.fsc._parse_deme_derived_allele_frequencies(
                    filepath=os.path.join(TESTS_DATA_DIR, "test-one", fixture["filename"]),
                    field_name_prefix=field_name_prefix,
                    results_d=data)
            expected_values = fixture["expected_values"]
            self.assertEqual(len(expected_values), len(data))
            for v1, v2 in zip(expected_values, data.values()):
                self.assertEqual(v1, v2)

    def test_joint_derived_allele_frequencies(self):
        fixtures = [
                {"filename": "test-one_jointMAFpop1_0.obs",
                 "expected_values": (
                        0   , 39 , 0   , 0   , 0 , 0   ,
                        918 , 11 , 171 , 0   , 0 , 0   ,
                        11  , 60 , 27  , 0   , 0 , 0   ,
                        0   , 0  , 40  , 0   , 0 , 0   ,
                        0   , 0  , 0   , 0   , 0 , 0   ,
                        0   , 0  , 0   , 101 , 0 , 0   ,
                        0   , 0  , 0   , 0   , 0 , 214 ,
                        0   , 0  , 0   , 0   , 0 , 2   ,
                        0   , 0  , 0   , 0   , 0 , 0   ,
                     )},
                ]
        for fidx, fixture in enumerate(fixtures):
            field_name_prefix = "stats.testv{}".format(fidx+1)
            data = collections.OrderedDict()
            self.fsc._parse_joint_derived_allele_frequencies(
                    filepath=os.path.join(TESTS_DATA_DIR, "test-one", fixture["filename"]),
                    field_name_prefix=field_name_prefix,
                    results_d=data)
            expected_values = fixture["expected_values"]
            self.assertEqual(len(expected_values), len(data))
            for v1, v2 in zip(expected_values, data.values()):
                self.assertEqual(v1, v2)

if __name__ == "__main__":
    unittest.main()

