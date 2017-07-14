
import unittest
from gerenuk import simulate

class MpArchitectureTests(unittest.TestCase):

    def test_run(self):
        config_d = {
                "name": "test",
                "standard_error_logging_level": "warning",
                "log_to_file": False,
                }
        gs = simulate.GerenukSimulator(
                config_d=config_d,
                num_processes=3,
                is_verbose_setup=False)
        gs.execute(10000)

if __name__ == "__main__":
    unittest.main()
