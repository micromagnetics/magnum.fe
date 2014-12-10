import unittest
from dolfin import *
from magnumfe import *
import numpy
import os

#set_log_active(False)

class DemagFieldTest(unittest.TestCase):

  def test_create_mesh(self):
    mesh, sample_size = DemagField.create_mesh((0.5, 0.5, 0.5), (10, 10, 10))
    self.assertEqual(mesh.with_shell.size(0), 1933)
    self.assertEqual(sample_size, [0.5, 0.5, 0.5])

  def test_create_mesh_with_domains(self):
    class TestDomain(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (0.1, 0.4))

    mesh, sample_size = DemagField.create_mesh((0.5, 0.5, 0.5), (10, 10, 10), domains = {2: TestDomain()})
    state = State(mesh, {'test': 2})
    self.assertAlmostEqual(0.3, assemble(Constant(1.0)*state.dx('test')))

  def test_energy_unit_cube(self):
    energy = self.energy_cube()
    self.assertTrue(abs(energy - 1.0/6.0) < 0.01)

  def test_energy_scaled_big_cube(self):
    energy = self.energy_cube(2)
    self.assertTrue(abs(energy - 8.0/6.0) < 0.08)

  def test_energy_scaled_small_cube(self):
    energy = self.energy_cube(0.5)
    self.assertTrue(abs(energy - 1.0/48.0) < 0.001)

  def energy_cube(self, scale=1.0):
    mesh, sample_size = DemagField.create_mesh((0.5, 0.5, 0.5), (10, 10, 10), d = 3, scale=scale)
    demag_field = DemagField("ST", sample_size, 2)

    state = State(mesh, m = Constant((0.0, 0.0, 1.0)))
    u = demag_field.calculate_potential(state)

    M = Constant(0.5) * inner(state.m, grad(u)) * dx(mesh)
    return assemble(M)
  
  def ttest_energy_sphere(self):
    sphere_mesh = os.path.dirname(os.path.realpath(__file__)) + "/mesh/sphere.msh"
    mesh, sample_size = DemagField.create_mesh(sphere_mesh, d=5, n=(10,10,10), margin=0.2)
    demag_field = DemagField("ST", sample_size, 2)

    state = State(mesh, m = Constant((0.0, 0.0, 1.0)))
    u = demag_field.calculate_potential(state)

    M = Constant(0.5) * inner(state.m, grad(u)) * state.dx('magnetic')
    energy = assemble(M)

    self.assertTrue(abs(energy - pi/36.0) < 0.01)

if __name__ == '__main__':
    unittest.main()
