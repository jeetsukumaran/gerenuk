
import unittest
import time
from collections import Counter
from gerenuk import simulate

class TestWorker(simulate.SimulationWorker):

    def simulate(self):
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
                }
        nreps = 10000
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
        # print(counter.keys()) # heisenbug wave-collapse!!
        self.assertEqual(len(counter), num_processes)

if __name__ == "__main__":
    unittest.main()
