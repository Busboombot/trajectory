import unittest
import sys
sys.path.append('/Users/eric/Documents/proj/trajectory/src/cmake-build-debug/python/')
import pyplanlib

class MyTestCase(unittest.TestCase):

    def test_something(self):
        print(dir(pyplanlib))


if __name__ == '__main__':
    unittest.main()
