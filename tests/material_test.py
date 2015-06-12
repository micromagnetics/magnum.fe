import unittest
from magnumfe import *

set_log_active(False)

class MaterialTest(unittest.TestCase):

  def test_init(self):
    mat = Material(ms = 8e5)
    self.assertAlmostEqual(8e5, mat.ms)

  def test_template_init(self):
    mat1 = Material(ms = 8e5)
    mat2 = Material(mat1, alpha = 1.0)
    self.assertAlmostEqual(8e5, mat2.ms)
    self.assertAlmostEqual(1.0, mat2.alpha)

if __name__ == '__main__':
    unittest.main()
