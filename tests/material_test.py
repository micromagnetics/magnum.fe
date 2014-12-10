import unittest
from dolfin import *
from magnumfe import *

set_log_active(False)

class MaterialTest(unittest.TestCase):

  def test_constructor(self):
    mat = Material(ms = 8e5)
    self.assertAlmostEqual(8e5, mat.ms)

if __name__ == '__main__':
    unittest.main()
