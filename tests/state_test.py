import unittest
from dolfin import *
from magnumfe import *

set_log_active(False)

class StateTest(unittest.TestCase):

  def test_domains_ids(self):
    state = State(UnitCubeMesh(1,1,1), celldomains = {'magnetic': 1, 'conducting': (1, 2), 'air': (3, 4)})

    self.assertEqual((1,),      state.domain_ids('magnetic'))
    self.assertEqual((1,2),     state.domain_ids('conducting'))
    self.assertEqual((1,2,3,4), state.domain_ids('all'))
    self.assertEqual((1,2),     state.domain_ids('!air'))
    self.assertEqual((2,3,4),   state.domain_ids('!magnetic'))

  def test_named_celldomains(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)})

    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('magnetic')), 6.4)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('air')), 1.6)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('!magnetic')), 1.6)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx(2)), 0.8)
    self.assertAlmostEqual(assemble(Constant(1.0)*state.dx('all')), 8.0)

  def test_named_facetdomains(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, facetdomains = {'outermagnet': 5})
    self.assertAlmostEqual(assemble(Constant(1.0)*state.ds('outermagnet')), 16.0)

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
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)})
    state.material['magnetic'] = Material(alpha = 1.0, k_axis = (0.0, 1.0, 0.0))
    #state.material['magnetic'] = Material() # TODO
    state.material['air'] = Material(alpha = 2.0, k_axis = (1.0, 0.0, 0.0))

    self.assertAlmostEqual(9.6, assemble(state.material.alpha * state.dx('all')))
    self.assertAlmostEqual(1.6, assemble(inner(state.material.k_axis, Constant((1.0, 0.0, 0.0))) * state.dx('all')))
    self.assertAlmostEqual(6.4, assemble(inner(state.material.k_axis, Constant((0.0, 1.0, 0.0))) * state.dx('all')))

  def test_volume(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)})
    self.assertAlmostEqual(6.4, state.volume('magnetic'))

  def test_average(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)}, m = Expression(("x[0]*x[0]", "0.0", "0.0")))

    avg = state.m.average('magnetic')
    self.assertAlmostEqual(0.42, avg[0])
    self.assertAlmostEqual(0.00, avg[1])
    self.assertAlmostEqual(0.00, avg[2])

  def test_crop(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)}, m = Expression(("x[0]*x[0]", "0.0", "0.0")))

    cropped = state.m.crop('magnetic')

    self.assertTrue(mesh.size(3) > cropped.function_space().mesh().size(3))
    self.assertAlmostEqual(0.42, assemble(cropped[0]*dx(cropped.function_space().mesh())) / state.volume('magnetic'))

  def test_normalize(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, {'everything': (1,2,3)})
    state.f = state.interpolate({
      1: Constant((0.0, 2.0, 0.0))
    })

    self.assertAlmostEqual(0.0, assemble(inner(Constant((1.0, 0.0, 0.0)), state.f) * state.dx()))
    self.assertAlmostEqual(14.4, assemble(inner(Constant((0.0, 1.0, 0.0)), state.f) * state.dx()))

    state.f.normalize()

    self.assertAlmostEqual(0.8, assemble(inner(Constant((1.0, 0.0, 0.0)), state.f) * state.dx()))
    self.assertAlmostEqual(7.2, assemble(inner(Constant((0.0, 1.0, 0.0)), state.f) * state.dx()))


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

  def test_interpolate(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh)
    f = state.interpolate({
      2: Constant(1.0),
      3: Constant(2.0)
    })
    self.assertAlmostEqual(assemble(f*state.dx(2)), 1.6)
    self.assertAlmostEqual(assemble(f*state.dx(3)), 2.0)

  def test_scale_parameter(self):
    mesh = UnitCubeMesh(1,1,1)

    state = State(mesh)
    self.assertAlmostEqual(state.scale, 1.0)

    state = State(mesh, scale=1e-9)
    self.assertAlmostEqual(state.scale, 1e-9)

  def test_step(self):
    class Integrator(object):
      def __init__(self):
        self.called = False
      def step(self, state, dt):
        self.called = True
    integrator = Integrator()

    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh)

    self.assertAlmostEqual(0.0, state.t)
    self.assertFalse(integrator.called)
    state.step(integrator, 1e-9)
    self.assertTrue(integrator.called)
    self.assertAlmostEqual(1e-9, state.t)

  def test_M_inv_diag(self):
    mesh = UnitCubeMesh(2,2,2)
    state = State(mesh)
    v = TestFunction(state.VectorFunctionSpace())
    r = assemble(inner(Constant((1.0, 0.0, 0.0)), v) * state.dx())
    result = state.M_inv_diag() * r
    self.assertAlmostEqual(27.0, result.array().sum())

  def test_M_inv_diag_on_subdomain(self):
    mesh = UnitCubeMesh(2,2,2)
    state = State(mesh)
    v = TestFunction(state.VectorFunctionSpace())
    r = assemble(inner(Constant((1.0, 0.0, 0.0)), v) * state.dx())
    result = state.M_inv_diag() * r
    self.assertAlmostEqual(27.0, result.array().sum())

  def mesh_with_subdomains(self):
    class TestDomain1(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (0.0, 0.25))

    class TestDomain2(SubDomain):
      def inside(self, x, on_boundary):
        return between(x[0], (-0.25, 0.0))

    class FacetDomain(SubDomain):
      def inside(self, x, on_boundary):
        return near(x[0], 2.0)

    mesher = Mesher()
    mesher.create_cuboid((1,1,1), (10,10,10))
    mesher.create_shell(1)
    mesher.create_celldomain(TestDomain1(), 2)
    mesher.create_celldomain(TestDomain2(), 3)
    mesher.create_facetdomain(FacetDomain(), 5)
    return mesher.mesh()

if __name__ == '__main__':
    unittest.main()
