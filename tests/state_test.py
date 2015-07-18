import unittest
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

  def test_ds_intersection(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, facetdomains = {'A': 5, 'B': 6, 'AB': (5,6)})
    ab           = assemble(Constant(1.0)*state.ds('AB'))
    a_intersect  = assemble(Constant(1.0)*state.ds('AB', intersect='A'))
    a            = assemble(Constant(1.0)*state.ds('A'))
    zero         = assemble(Constant(1.0)*state.ds('A', intersect='B'))

    self.assertAlmostEqual(2*a, ab)
    self.assertAlmostEqual(a, a_intersect)
    self.assertAlmostEqual(0.0, zero)

  def test_attribute_init(self):
    mesh = UnitCubeMesh(1,1,1)
    m = Constant((1.0, 0.0, 0.0))
    state = State(mesh, m = m)

    self.assertAlmostEqual(1.0, assemble(inner(state.m, m)*state.dx()))

  def test_attribute_setter_should_reuse_function(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh, m = Constant((1.0, 0.0, 0.0)))
    self.assertAlmostEqual(1.0, state.m.average()[0])
    m = state.m
    state.m = Constant((2.0, 0.0, 0.0))
    self.assertEqual(m, state.m)
    self.assertAlmostEqual(2.0, state.m.average()[0])

  def test_attribute_uuid(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh, m = Constant((1.0, 0.0, 0.0)))
    uuid = state.m.uuid
    self.assertEqual(uuid, state.m.uuid)
    state.m = Constant((2.0, 0.0, 0.0))
    self.assertNotEqual(uuid, state.m.uuid)

  def test_multi_uuid(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh, a1 = 1.0, a2 = 2.0)
    self.assertEqual([1.0, 2.0], state.uuid("a1", "a2"))
    state.a2 = 3.0
    self.assertEqual([1.0, 3.0], state.uuid("a1", "a2"))

  def test_lambda_attribute_cycle(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh,
        a = lambda state: (state.b, "b"),
        b = lambda state: (state.a, "a")
    )
    self.assertRaises(AttributeCycleDependency, lambda:state.a)

  def test_lambda_attribute(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh, t2 = lambda state: 2.0 * state.t)
    state.t = 1.0
    self.assertAlmostEqual(2.0, state.t2)
    state.t = 2.0
    self.assertAlmostEqual(4.0, state.t2)

  def test_lambda_attribute_decoration(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh,
        celldomains = {'a': 2, 'b': 3},
        f           = lambda state: Expression("x[0]")
    )
    f = state.f.crop('a')
    f_int = assemble(f*dx(f.function_space().mesh()))
    self.assertAlmostEqual(0.1, f_int/state.volume('a'))
    self.assertAlmostEqual(0.1, state.f.average('a'))

  def test_lambda_attribute_reset(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh)
    state.a = lambda state: state.interpolate(Constant(state.t))
    self.assertTrue(isinstance(state.a, Function))
    state.t = 1.0
    self.assertTrue(isinstance(state.a, Function))

  def test_lambda_attribute_caching(self):
    global count
    count = 0
    def lamb(state):
      global count
      count += 1
      return (2.0 * state.t, "t", "m")

    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh, lamb = lamb, m = Constant((1.0, 0.0, 0.0)))
    state.lamb
    state.lamb
    self.assertEqual(1, count)
    state.t = 1.0
    state.lamb
    state.lamb
    self.assertEqual(2, count)
    state.m = Constant((0.0, 1.0, 0.0))
    state.lamb
    state.lamb
    self.assertEqual(3, count)

  def test_lambda_scalar_attribute_caching(self):
    global count
    count = 0
    def l1(state):
      global count
      count += 1
      return (2.0 * state.t, "l2")

    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh, l1=l1, l2=lambda state: int(state.t))

    state.t = 1.0
    state.l1
    self.assertAlmostEqual(2.0, state.l1)
    self.assertEqual(1, count)
    state.t = 1.5
    self.assertAlmostEqual(2.0, state.l1)
    self.assertEqual(1, count)
    state.t = 2.5
    self.assertAlmostEqual(5.0, state.l1)
    self.assertEqual(2, count)

  def test_lambda_attribute_uuid(self):
    global count
    count = 0
    def lamb(state):
      global count
      count += 1
      return (Constant(2.0 * state.t), "t", "m")

    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh, lamb = lamb, m = Constant((1.0, 0.0, 0.0)))
    uuid = state.uuid('lamb')
    self.assertEqual(uuid, state.uuid('lamb'))
    state.m = Constant((0.0, 1.0, 0.0))
    self.assertNotEqual(uuid, state.uuid('lamb'))

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

  def test_cascaded_material_assignment(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': (1, 2, 3), 'm1': 2})
    state.material['magnetic'] = Material(ms=8e5,alpha=1.0)
    state.material['m1'] = Material(ms=4e5)

    self.assertAlmostEqual(state.volume(), assemble(state.material.alpha*state.dx()))

  def test_dP_functional(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)}, a = Constant(1.0))

    self.assertAlmostEqual(state.a.crop('magnetic').function_space().mesh().size(0), assemble(state.a*state.dP('magnetic')))
    self.assertAlmostEqual(state.a.crop('air').function_space().mesh().size(0), assemble(state.a*state.dP('air')))

  def test_dP_vector(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)}, a = Constant((1.0, 2.0, 3.0)))

    v = TestFunction(state.VectorFunctionSpace())
    result = assemble(inner(v, state.a)*state.dP('magnetic')).array()

    self.assertAlmostEqual(state.a.crop('magnetic').function_space().mesh().size(0)*6, result.sum())

  def test_material_assignment_with_base(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains = {'magnetic': 1, 'air': (2, 3)}, material = Material(alpha = 0.1))
    state.material['magnetic'] = Material(alpha = 1.0)

    self.assertAlmostEqual(6.56, assemble(state.material.alpha*state.dx()))

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
    state.material = Material(alpha = 0.1)
    self.assertTrue(isinstance(state.material.alpha, Constant))
    self.assertAlmostEqual(0.1, assemble(state.material.alpha*state.dx()))
    state.material.alpha = 0.2
    self.assertAlmostEqual(0.2, assemble(state.material.alpha*state.dx()))

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
    self.assertAlmostEqual(assemble(f*state.dx(2)), 1.2)
    self.assertAlmostEqual(assemble(f*state.dx(3)), 1.6)

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
    return Mesh("mesh/state.xml.gz")

if __name__ == '__main__':
    unittest.main()
