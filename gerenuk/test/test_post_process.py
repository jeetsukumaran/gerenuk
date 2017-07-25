#! /usr/bin/env python

import sys
import csv
import unittest
import os
from gerenuk import utility
from gerenuk.utility import StringIO
from gerenuk.test import TESTS_DATA_DIR

class FilterColumnsTestCase(unittest.TestCase):

    def test_filter_columns_from_template_file(self):
        master = (
        "c1|c2|c3|c4|c5|c6|c7|c8",
        "11|12|13|14|15|16|17|18",
        "21|22|23|24|25|26|27|28",
        "31|32|33|34|35|36|37|38",
        "41|42|43|44|45|46|47|48",
        "51|52|53|54|55|56|57|58",
        "61|62|63|64|65|66|67|68",
        "71|72|73|74|75|76|77|78",
        "81|82|83|84|85|86|87|78",
        )
        target = (
        "c1|x2|c3|x4|c5|x6|c7|x8",
        "11|12|13|14|15|16|17|18",
        "21|22|23|24|25|26|27|28",
        "31|32|33|34|35|36|37|38",
        "41|42|43|44|45|46|47|48",
        "51|52|53|54|55|56|57|58",
        "61|62|63|64|65|66|67|68",
        "71|72|73|74|75|76|77|78",
        "81|82|83|84|85|86|87|78",
        )
        expected = (
        "c1|c3|c5|c7",
        "11|13|15|17",
        "21|23|25|27",
        "31|33|35|37",
        "41|43|45|47",
        "51|53|55|57",
        "61|63|65|67",
        "71|73|75|77",
        "81|83|85|87",
        )
        dest = StringIO()
        master_file = StringIO("\n".join(master).replace("|", "\t"))
        target_file = StringIO("\n".join(target).replace("|", "\t"))
        utility.filter_columns_using_master_template_file(
                dest=dest,
                master_file=master_file,
                source_file=target_file)
        result = dest.getvalue().strip()
        expected_str = "\n".join(expected).replace("|", "\t")
        self.assertEqual(result, expected_str)

if __name__ == "__main__":
    unittest.main()


