import unittest
from dolfin import *
from magnumfe import *

set_log_active(False)

class DemagFieldTestFK(unittest.TestCase):

    def test_energy_unit_cube(self):
      mesh = UnitCubeMesh(10, 10, 10)
      demag_field = DemagFieldFK()
      state = State(mesh, m = Constant((0.0, 0.0, 1.0)))
      u = demag_field.calculate_potential(state)
      energy = assemble(Constant(0.5) * inner(state.m, grad(u)) * dx(mesh))
      self.assertTrue(abs(energy - 1.0/6.0) < 0.01)

if __name__ == '__main__':
    unittest.main()
