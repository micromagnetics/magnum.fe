import unittest
from dolfin import *
from magnumfe import *
import numpy


class MesherTest(unittest.TestCase):
    def test_create_shell(self):
      mesher = Mesher()
      # create and mesh sample
      mesher.create_cuboid((3.0, 2.0, 1.0), (15, 10, 5))

      # create and mesh shell
      mesher.create_shell(1);
      mesh_with_shell = mesher.mesh()
      self.assertEqual(mesh_with_shell.num_cells(), 5640)
      self.assertEqual(mesh_with_shell.num_vertices(), 1188)

    def test_get_sample_size(self):
      mesher = Mesher()
      mesher.create_cuboid((3.0, 2.0, 1.0), (15, 10, 5))
      size = mesher.get_sample_size()

      self.assertEqual(size[0], 3.0)
      self.assertEqual(size[1], 2.0)

if __name__ == '__main__':
    unittest.main()
