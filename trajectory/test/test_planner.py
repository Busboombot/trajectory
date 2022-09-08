import unittest
from trajectory.iplanner import *

class TestPlanner(unittest.TestCase):



    def setUp(self) -> None:
        self.a_max = 50_000
        self.v_max = 5_000

        self.j = Joint(self.v_max, self.a_max)

    def test_something(self):
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
