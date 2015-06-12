import unittest
from magnumfe import *

set_log_active(False)

class DemagFieldTestFK(unittest.TestCase):

  def test_energy_unit_cube(self):
    mesh = UnitCubeMesh(10, 10, 10)
    demag_field = DemagField("FK", loglevel = 'low')
    state = State(mesh, m = Constant((0.0, 0.0, 1.0)))
    u = demag_field.u(state)
    energy = assemble(Constant(0.5) * inner(state.m, grad(u)) * dx(mesh))
    self.assertTrue(abs(energy - 1.0/6.0) < 0.01)

  def test_state_switch(self):
    demag_field = DemagField("FK", loglevel = 'low')
    mesh1 = UnitCubeMesh(5,5,5)
    mesh2 = UnitCubeMesh(6,6,6)

    state1 = State(mesh1, m = Constant((0.0, 0.0, 1.0)))
    state2 = State(mesh2, m = Constant((0.0, 1.0, 0.0)))

    u1 = demag_field.u(state1)
    u2 = demag_field.u(state2)

    energy1 = assemble(Constant(0.5) * inner(state1.m, grad(u1)) * dx(mesh1))
    energy2 = assemble(Constant(0.5) * inner(state2.m, grad(u2)) * dx(mesh2))

    self.assertAlmostEqual(energy1, energy2, 2)

  def test_multi_domain(self):
    class UnitCubeDomain(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (-0.50, 0.50)) and \
               between(x[1], (-0.50, 0.50)) and \
               between(x[2], (-0.50, 0.50))

    cube_domain = UnitCubeDomain()

    mesh = Mesh("mesh/demag_fk_multi.xml.gz")

    state = State(mesh, {'magnetic': 2}, m = Constant((0.0, 0.0, 1.0)))
    demag_field = DemagField("FK", loglevel = 'low')
    u = demag_field.u(state)
    energy = assemble(Constant(0.5) * inner(state.m, grad(u)) * state.dx('magnetic'))
    self.assertTrue(abs(energy - 1.0/6.0) < 0.01)

  def test_field_with_lumping(self):
    class UnitCubeDomain(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (-0.50, 0.50)) and \
               between(x[1], (-0.50, 0.50)) and \
               between(x[2], (-0.50, 0.50))

    cube_domain = UnitCubeDomain()

    mesh = Mesh("mesh/demag_fk_lumping.xml.gz")

    state = State(mesh, {'magnetic': 2},
              m = Constant((0.0, 0.0, 1.0)),
              material = Material(ms = 1.0)
            )
    demag_field = DemagField("FK", loglevel = 'low')

    h = demag_field.field(state, lump = True)
    energy = assemble(Constant(0.5) * inner(state.m, - h) * state.dx('magnetic'))
    self.assertTrue(abs(energy - 1.0/6.0) < 0.01)

if __name__ == '__main__':
    unittest.main()
