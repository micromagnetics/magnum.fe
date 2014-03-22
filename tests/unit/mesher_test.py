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
      mesher.create_shell(1)
      mesh_with_shell = mesher.mesh()
      self.assertEqual(mesh_with_shell.num_cells(), 5640)
      self.assertEqual(mesh_with_shell.num_vertices(), 1188)

    def test_get_sample_size(self):
      mesher = Mesher()
      mesher.create_cuboid((3.0, 2.0, 1.0), (15, 10, 5))
      size = mesher.get_sample_size()

      self.assertEqual(size[0], 3.0)
      self.assertEqual(size[1], 2.0)

    def test_get_scaled_sample_size(self):
      mesher = Mesher()
      mesher.create_cuboid((3.0, 2.0, 1.0), (15, 10, 5))
      size = mesher.get_sample_size(scale = 10)

      self.assertEqual(size[0], 30.0)
      self.assertEqual(size[1], 20.0)

    def test_scale(self):
      mesher = Mesher()
      mesher.create_cuboid((3.0, 2.0, 1.0), (15, 10, 5))

      mesh1 = mesher.mesh()
      mesh2 = mesher.mesh(10.0)

      self.assertAlmostEqual(mesh1.hmax()*10.0, mesh2.hmax())
      self.assertAlmostEqual(mesh1.hmin()*10.0, mesh2.hmin())

    def test_read_mesh_with_domains(self):
      # TODO test domain stuff
      #mesher = Mesher()
      #mesher.read_file("mesh/sphere.msh")
      #mesher.create_cuboid((1,1,1),(1,1,1))
      #mesher.create_shell(1, n=(1,1,1))
      #mesh = mesher.mesh()
      #f = File("mesh.xml")
      #f << mesh

      #mesher = Mesher()
      #mesher.read_file("mesh/multidomain.msh")
      #mesher.create_shell(1, n=(10,10,10), margin=0.1)
      #mesh = mesher.mesh()
      #f = File("mesh2.xml")
      #f << mesh
      #domains = MeshFunction("size_t", mesh, 3, mesh.domains())
      #f = File("domains.pvd")
      #f << domains
      #mesher.write_file("bla.msh")

if __name__ == '__main__':
    unittest.main()
