
import unittest
import time
from collections import Counter
from gerenuk import simulate

class TestWorker(simulate.SimulationWorker):

    def simulate(self):
        time.sleep(0.0001)
        return {
                "name": self.name,
                "task_count": self.num_tasks_received,
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
        for result in results:
            counter[result["name"]] += 1
        self.assertEqual(len(counter), num_processes)

if __name__ == "__main__":
    unittest.main()
