import unittest
from dolfin import *
from magnumfe import *

set_log_active(False)

class StateTest(unittest.TestCase):
  def test_domains_ids(self):
    state = State(UnitCubeMesh(1,1,1), regions = {'magnetic': 1, 'conducting': (1, 2), 'air': (3, 4)})

    self.assertEqual((1,),      state.domain_ids('magnetic'))
    self.assertEqual((1,2),     state.domain_ids('conducting'))
    self.assertEqual((1,2,3,4), state.domain_ids('all'))
    self.assertEqual((1,2),     state.domain_ids('!air'))
    self.assertEqual((2,3,4),   state.domain_ids('!magnetic'))

  def test_named_regions(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, regions = {'magnetic': 1, 'air': (2, 3)})

    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('magnetic')), 6.4)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('air')), 1.6)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('!magnetic')), 1.6)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx(2)), 0.8)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('all')), 8.0)

  def test_attribute_init(self):
    mesh = UnitCubeMesh(1,1,1)
    m = Constant((1.0, 0.0, 0.0))
    state = State(mesh, m = m)

    self.assertTrue(isinstance(state.m, Function))
    self.assertAlmostEqual(1.0, assemble(inner(state.m, m)*state.dx()))

  def test_dx_if_no_domains_defined(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh)

    self.assertAlmostEqual(1.0, assemble(Constant(1.0) * state.dx()))

  def test_material_assignment(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, regions = {'magnetic': 1, 'air': (2, 3)})
    state.material['magnetic'] = Material(alpha = 1.0, k_axis = (0.0, 1.0, 0.0))
    #state.material['magnetic'] = Material() # TODO
    state.material['air'] = Material(alpha = 2.0, k_axis = (1.0, 0.0, 0.0))

    self.assertAlmostEqual(9.6, assemble(state.material.alpha * state.dx('all')))
    self.assertAlmostEqual(1.6, assemble(inner(state.material.k_axis, Constant((1.0, 0.0, 0.0))) * state.dx('all')))
    self.assertAlmostEqual(6.4, assemble(inner(state.material.k_axis, Constant((0.0, 1.0, 0.0))) * state.dx('all')))

  def test_set_global_material(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh)
    state.material = Material(alpha = 1)
    self.assertTrue(isinstance(state.material.alpha, Constant))

  def test_mesh_has_domains(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh)
    self.assertFalse(state.mesh_has_domains())

    mesh = self.mesh_with_subdomains()
    state = State(mesh)
    self.assertTrue(state.mesh_has_domains())

  def mesh_with_subdomains(self):
    class TestDomain1(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (0.0, 0.25))

    class TestDomain2(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (-0.25, 0.0))

    test_domain1 = TestDomain1()
    test_domain2 = TestDomain2()

    mesher = Mesher()
    mesher.create_cuboid((1,1,1), (10,10,10))
    mesher.create_shell(1)
    mesher.create_celldomain(test_domain1, 2)
    mesher.create_celldomain(test_domain2, 3)
    return mesher.mesh()

if __name__ == '__main__':
    unittest.main()
