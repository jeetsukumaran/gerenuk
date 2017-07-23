
import unittest
import time
from collections import Counter
from gerenuk import simulate

class TestWorker(simulate.SimulationWorker):

    def simulate(self):
        time.sleep(0.2)
        return {
                "name": self.name,
                "task_count": self.num_tasks_received,
                "rand_int": self.model.rng.randint(1, 1E6),
                }

class MpArchitectureTests(unittest.TestCase):

    def test_workers_used(self):
        config_d = {
                "name": "test",
                "standard_error_logging_level": "warning",
                "log_to_file": False,
                "params": {
                    "concentrationShape": 10,
                    "concentrationScale": 0.3766,
                    "thetaShape": 1,
                    "thetaScale": 0.03,
                    # "popsizeShape": 2,
                    # "popsizeScale": 1E7,
                    # "mutRateShape": 2,
                    # "mutRateScale": 2E-8,
                    "ancestralThetaShape": 0,
                    "ancestralThetaScale": 0,
                    "thetaParameters": "012",
                    "tauShape": 1.0,
                    "tauScale": 0.007,
                    "timeInSubsPerSite": 1,
                    "bottleProportionShapeA": 0,
                    "bottleProportionShapeB": 0,
                    "bottleProportionShared": 0,
                    "migrationShape": 0,
                    "migrationScale": 0,
                    "numTauClasses": 0
                    },
                "locus_info": [
                    {'taxon_label': 'S1',
                        'locus_label': 'LocusS1M1',
                        'ploidy_factor': 1,
                        'mutation_rate_factor': 1,
                        'num_genes_deme0': 38,
                        'num_genes_deme1': 38,
                        'ti_tv_rate_ratio': 3,
                        'num_sites': 80,
                        'freq_a': 0.263,
                        'freq_c': 0.258,
                        'freq_g': 0.255,
                        'alignment_filepath': 'S1M1.fasta',
                        },

                    {'taxon_label': 'S1',
                        'locus_label': 'LocusS1M2',
                        'ploidy_factor': 1,
                        'mutation_rate_factor': 1,
                        'num_genes_deme0': 38,
                        'num_genes_deme1': 36,
                        'ti_tv_rate_ratio': 3,
                        'num_sites': 80,
                        'freq_a': 0.149,
                        'freq_c': 0.146,
                        'freq_g': 0.226,
                        'alignment_filepath': 'S1M2.fasta',
                        },
                    ]
                }
        nreps = 10
        num_processes = 2
        gs = simulate.GerenukSimulator(
                config_d=config_d,
                num_processes=num_processes,
                is_verbose_setup=False)
        gs.worker_class = TestWorker
        results = gs.execute(nreps)
        self.assertEqual(len(results), nreps)
        counter = Counter()
        for result_idx, result in enumerate(results):
            counter[result["name"]] += 1
        self.assertEqual(len(counter), num_processes)

if __name__ == "__main__":
    unittest.main()
