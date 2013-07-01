import unittest
from dolfin import *
from magnumfe import *
import numpy

class DemagFieldCRTest(unittest.TestCase):

    def test_create_mesh(self):
      mesh = DemagFieldCR.create_mesh((0.5, 0.5, 0.5), (11, 11, 11))
      self.assertEqual(mesh.size(0), 1331)
      self.assertEqual(mesh.with_shell.size(0), 1933)
      self.assertEqual(mesh.data['sample_size'], [0.5, 0.5, 0.5])

    def test_energy_unit_cube(self):
      mesh = DemagFieldCR.create_mesh((0.5, 0.5, 0.5), (10, 10, 10), d = 3)
      VS = FunctionSpace(mesh, "Lagrange", 2)
      VV = VectorFunctionSpace(mesh, "Lagrange", 1)

      m = interpolate(Constant((0.0, 0.0, 1.0)), VV)
      demag = DemagFieldCR(mesh, order = 1)
      u = demag.calculate(m)
      M = 0.5 * inner(m, grad(u)) * dx
      energy1 = assemble(M)

      demag = DemagField(mesh, order = 2)
      u = demag.calculate(m)
      M = 0.5 * inner(m, grad(u)) * dx
      energy2 = assemble(M)

      print "CR: %f ; CG: %f" % (energy1, energy2)
      self.assertTrue(abs(energy1 - 1.0/6.0) < 0.01)

if __name__ == '__main__':
    unittest.main()
