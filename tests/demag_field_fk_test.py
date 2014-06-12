import unittest
from dolfin import *
from magnumfe import *

set_log_active(False)

class DemagFieldTestFK(unittest.TestCase):

  def test_energy_unit_cube(self):
    mesh = UnitCubeMesh(10, 10, 10)
    demag_field = DemagField("FK")
    state = State(mesh, m = Constant((0.0, 0.0, 1.0)))
    u = demag_field.calculate_potential(state)
    energy = assemble(Constant(0.5) * inner(state.m, grad(u)) * dx(mesh))
    self.assertTrue(abs(energy - 1.0/6.0) < 0.01)

  def test_multi_domain(self):
    class UnitCubeDomain(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (-0.50, 0.50)) and \
               between(x[1], (-0.50, 0.50)) and \
               between(x[2], (-0.50, 0.50))

    cube_domain = UnitCubeDomain()

    mesher = Mesher()
    mesher.create_cuboid((1.0, 1.0, 1.0), (20, 20, 20))
    mesher.create_celldomain(cube_domain, 2)
    mesh = mesher.mesh()

    state = State(mesh, {'magnetic': 2}, m = Constant((0.0, 0.0, 1.0)))
    demag_field = DemagField("FK")
    u = demag_field.calculate_potential(state)
    energy = assemble(Constant(0.5) * inner(state.m, grad(u)) * state.dx('magnetic'))
    self.assertTrue(abs(energy - 1.0/6.0) < 0.01)

if __name__ == '__main__':
    unittest.main()
