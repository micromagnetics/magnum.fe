import unittest
from dolfin import *
from magnumfe import *
import numpy

mesher = Mesher()
mesher.create_cuboid((1,1,1), (7,7,7))
mesher.create_shell(2)
complete_mesh = mesher.mesh()

mesh = WrappedMesh.create(complete_mesh, 0)

class WrappedMeshTest(unittest.TestCase):

    def test_cut(self):
      V    = FunctionSpace(mesh.with_shell, "Lagrange", 2)
      expr = Expression("sin(x[0])")
      f    = interpolate(expr, V)


      f_cut = mesh.cut(f)
      
      self.assertOnSuperMesh(f)
      self.assertOnSubMesh(f_cut)

      self.assertEqualAtPoint(f, f_cut, (0.1, 0.2, 0.3))
      self.assertEqualAtPoint(f, f_cut, (0.3, 0.5, 0.8))
      self.assertEqualAtPoint(f, f_cut, (-0.1, -0.2, 0.3))

    def test_expand(self):
      V    = FunctionSpace(mesh, "Lagrange", 2)
      expr = Expression("sin(x[0])")
      f    = interpolate(expr, V)

      f_expanded = mesh.expand(f)
      
      self.assertOnSubMesh(f)
      self.assertOnSuperMesh(f_expanded)

      self.assertEqualAtPoint(f, f_expanded, (0.1, 0.2, 0.3))
      self.assertEqualAtPoint(f, f_expanded, (0.3, 0.5, 0.8))
      self.assertEqualAtPoint(f, f_expanded, (-0.1, -0.2, 0.3))
      self.assertZeroAtPoint(f_expanded, (1.6, 1.6, 1.6))
      self.assertZeroAtPoint(f_expanded, (-1.6, -1.6, -1.6))

    def assertEqualAtPoint(self, f1, f2, point):
      v1 = numpy.zeros((1,), dtype="d")
      v2 = numpy.zeros((1,), dtype="d")
      f1.eval(v1, numpy.array(point))
      f2.eval(v2, numpy.array(point))
      self.assertAlmostEqual(v1[0], v2[0])

    def assertZeroAtPoint(self, f, point):
      v = numpy.zeros((1,), dtype="d")
      f.eval(v, numpy.array(point))
      self.assertEqual(0, v[0])

    def assertOnSuperMesh(self, f):
      self.assertEqual(mesh.with_shell.size(0), f.function_space().mesh().size(0))

    def assertOnSubMesh(self, f):
      self.assertEqual(mesh.size(0), f.function_space().mesh().size(0))

if __name__ == '__main__':
    unittest.main()
