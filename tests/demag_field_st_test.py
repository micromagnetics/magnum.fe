import unittest
from magnumfe import *
import numpy
import os

set_log_active(False)

class DemagFieldTest(unittest.TestCase):

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
    mesh = DemagField.wrap_mesh(os.path.dirname(os.path.realpath(__file__)) + "/mesh/demag_st_cube_%d.xml.gz" % round(scale*10))
    demag_field = DemagField("ST", 2)

    state = State(mesh, m = Constant((0.0, 0.0, 1.0)))
    u = demag_field.u(state)

    M = Constant(0.5) * inner(state.m, grad(u)) * dx(mesh)
    return assemble(M)
  
  def test_energy_sphere(self):
    mesh = DemagField.wrap_mesh(Mesh(os.path.dirname(os.path.realpath(__file__)) + "/mesh/demag_st_sphere.xml.gz"))
    demag_field = DemagField("ST", 2)

    state = State(mesh, m = Constant((0.0, 0.0, 1.0)))
    u = demag_field.u(state)

    M = Constant(0.5) * inner(state.m, grad(u)) * state.dx('magnetic')
    energy = assemble(M)

    self.assertTrue(abs(energy - pi/36.0) < 0.01)

if __name__ == '__main__':
    unittest.main()
