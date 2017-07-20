#! /usr/bin/env python

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


if __name__ == "__main__":
    unittest.main()

