#! /usr/bin/env python

import unittest
import os
from gerenuk import utility
from gerenuk.test import TESTS_DATA_DIR
CONFIG_DATA_DIR = os.path.join(TESTS_DATA_DIR, "configuration-files")

class Fsc2DataExtractionTestCase(unittest.TestCase):

    def test_parse_configuration(self):
        config_filepath = os.path.join(CONFIG_DATA_DIR, "sample1.cfg")
        config_d = utility.parse_legacy_configuration(config_filepath)
        self.assertIn("params", config_d)
        expected_params = {
            "concentrationShape": 1000.0,
            "concentrationScale": 0.00437,
            "thetaShape": 4.0,
            "thetaScale": 0.001,
            "ancestralThetaShape": 0,
            "ancestralThetaScale": 0,
            "thetaParameters": 000,
            "tauShape": 1.0,
            "tauScale": 0.02,
            "timeInSubsPerSite": 1,
            "bottleProportionShapeA": 0,
            "bottleProportionShapeB": 0,
            "bottleProportionShared": 0,
            "migrationShape": 0,
            "migrationScale": 0,
            "numTauClasses": 0,}
        self.assertEqual(len(expected_params), len(config_d["params"]))
        for key in expected_params:
            self.assertIn(key, config_d["params"])
            self.assertEqual(expected_params[key], config_d["params"][key])
        for key in config_d["params"]:
            self.assertIn(key, expected_params)
        self.assertIn("locus_info", config_d)
        expected_locus_info = [
            {'taxon_label': 'species-1', 'locus_label': 'locus-1', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 10, 'num_genes_deme1': 8, 'ti_tv_rate_ratio': 32.42, 'num_sites': 389, 'freq_a': 0.27, 'freq_c': 0.24, 'freq_g': 0.26, 'alignment_filepath': 'species-1-locus-1.fasta', },
            {'taxon_label': 'species-1', 'locus_label': 'locus-2', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 8, 'num_genes_deme1': 6, 'ti_tv_rate_ratio': 5.51, 'num_sites': 500, 'freq_a': 0.25, 'freq_c': 0.22, 'freq_g': 0.24, 'alignment_filepath': 'species-1-locus-2.fasta', },
            {'taxon_label': 'species-1', 'locus_label': 'locus-3', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 6, 'num_genes_deme1': 8, 'ti_tv_rate_ratio': 8.38, 'num_sites': 524, 'freq_a': 0.26, 'freq_c': 0.23, 'freq_g': 0.26, 'alignment_filepath': 'species-1-locus-3.fasta', },
            {'taxon_label': 'species-1', 'locus_label': 'locus-4', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 8, 'num_genes_deme1': 10, 'ti_tv_rate_ratio': 5.20, 'num_sites': 345, 'freq_a': 0.25, 'freq_c': 0.23, 'freq_g': 0.24, 'alignment_filepath': 'species-1-locus-4.fasta', },
            {'taxon_label': 'species-1', 'locus_label': 'locus-5', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 8, 'num_genes_deme1': 8, 'ti_tv_rate_ratio': 29.59, 'num_sites': 417, 'freq_a': 0.27, 'freq_c': 0.23, 'freq_g': 0.21, 'alignment_filepath': 'species-1-locus-5.fasta', },
            {'taxon_label': 'species-1', 'locus_label': 'mito-1', 'ploidy_factor': 0.25, 'mutation_rate_factor': 4.0, 'num_genes_deme0': 5, 'num_genes_deme1': 5, 'ti_tv_rate_ratio': 8.15, 'num_sites': 600, 'freq_a': 0.22, 'freq_c': 0.24, 'freq_g': 0.27, 'alignment_filepath': 'species-1-mito-1.fasta', },
            {'taxon_label': 'species-2', 'locus_label': 'locus-1', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 6, 'num_genes_deme1': 10, 'ti_tv_rate_ratio': 7.53, 'num_sites': 400, 'freq_a': 0.25, 'freq_c': 0.24, 'freq_g': 0.26, 'alignment_filepath': 'species-2-locus-1.fasta', },
            {'taxon_label': 'species-2', 'locus_label': 'locus-3', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 10, 'num_genes_deme1': 8, 'ti_tv_rate_ratio': 11.14, 'num_sites': 550, 'freq_a': 0.27, 'freq_c': 0.22, 'freq_g': 0.24, 'alignment_filepath': 'species-2-locus-3.fasta', },
            {'taxon_label': 'species-2', 'locus_label': 'locus-4', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 8, 'num_genes_deme1': 8, 'ti_tv_rate_ratio': 9.39, 'num_sites': 350, 'freq_a': 0.24, 'freq_c': 0.24, 'freq_g': 0.23, 'alignment_filepath': 'species-2-locus-4.fasta', },
            {'taxon_label': 'species-2', 'locus_label': 'locus-5', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 10, 'num_genes_deme1': 10, 'ti_tv_rate_ratio': 13.32, 'num_sites': 450, 'freq_a': 0.26, 'freq_c': 0.24, 'freq_g': 0.22, 'alignment_filepath': 'species-2-locus-5.fasta', },
            {'taxon_label': 'species-2', 'locus_label': 'mito-1', 'ploidy_factor': 0.25, 'mutation_rate_factor': 4.0, 'num_genes_deme0': 4, 'num_genes_deme1': 5, 'ti_tv_rate_ratio': 7.59, 'num_sites': 549, 'freq_a': 0.23, 'freq_c': 0.26, 'freq_g': 0.23, 'alignment_filepath': 'species-2-mito-1.fasta', },
            {'taxon_label': 'species-3', 'locus_label': 'locus-1', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 10, 'num_genes_deme1': 6, 'ti_tv_rate_ratio': 17.03, 'num_sites': 367, 'freq_a': 0.25, 'freq_c': 0.23, 'freq_g': 0.27, 'alignment_filepath': 'species-3-locus-1.fasta', },
            {'taxon_label': 'species-3', 'locus_label': 'locus-3', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 8, 'num_genes_deme1': 10, 'ti_tv_rate_ratio': 59.17, 'num_sites': 541, 'freq_a': 0.26, 'freq_c': 0.22, 'freq_g': 0.25, 'alignment_filepath': 'species-3-locus-3.fasta', },
            {'taxon_label': 'species-3', 'locus_label': 'locus-4', 'ploidy_factor': 1.0, 'mutation_rate_factor': 1.0, 'num_genes_deme0': 6, 'num_genes_deme1': 8, 'ti_tv_rate_ratio': 6.90, 'num_sites': 333, 'freq_a': 0.28, 'freq_c': 0.23, 'freq_g': 0.21, 'alignment_filepath': 'species-3-locus-4.fasta', },
            {'taxon_label': 'species-3', 'locus_label': 'mito-1', 'ploidy_factor': 0.25, 'mutation_rate_factor': 4.0, 'num_genes_deme0': 5, 'num_genes_deme1': 4, 'ti_tv_rate_ratio': 11.42, 'num_sites': 587, 'freq_a': 0.22, 'freq_c': 0.22, 'freq_g': 0.25, 'alignment_filepath': 'species-3-mito-1.fasta', },
        ]
        self.assertEqual(len(expected_locus_info), len(config_d["locus_info"]))
        for exp_locus_definition, obs_locus_definition in zip(expected_locus_info, config_d["locus_info"]):
            self.assertEqual(len(exp_locus_definition), len(obs_locus_definition))
            for key in exp_locus_definition:
                self.assertIn(key, obs_locus_definition)
                self.assertEqual(exp_locus_definition[key], obs_locus_definition[key])
            for key in obs_locus_definition:
                self.assertIn(key, exp_locus_definition)

if __name__ == "__main__":
    unittest.main()
