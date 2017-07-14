
import unittest
from gerenuk import simulate

class TestWorker(simulate.SimulationWorker):

    def simulate(self):
        return [self.name, self.num_tasks_received]

class MpArchitectureTests(unittest.TestCase):

    def test_run(self):
        config_d = {
                "name": "test",
                "standard_error_logging_level": "warning",
                "log_to_file": False,
                }
        nreps = 1000
        num_processes = 10
        gs = simulate.GerenukSimulator(
                config_d=config_d,
                num_processes=num_processes,
                is_verbose_setup=False)
        gs.worker_class = TestWorker
        results = gs.execute(nreps)
        self.assertEqual(len(results), nreps)

if __name__ == "__main__":
    unittest.main()
