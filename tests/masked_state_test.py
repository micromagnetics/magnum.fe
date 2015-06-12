import unittest
from magnumfe import *

set_log_active(False)

class MaskedStateTest(unittest.TestCase):

  def test_mesh_size(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh)
    masked = MaskedState(state, 2)
    self.assertEqual(9*9*3, masked.mesh.size(0))
    self.assertEqual(9*9*9, state.mesh.size(0))

  def test_function_getter(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh)
    state.a = Constant(1.0)

    masked = MaskedState(state, 2)
    self.assertEqual(9*9*3, masked.a.function_space().mesh().size(0))

  def test_function_getter_cache(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh)
    state.a = Constant(1.0)

    masked = MaskedState(state, 2)
    a1 = masked.a
    a2 = masked.a

    self.assertEqual(a1, a2)

  def test_function_getter_uuid(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh)
    state.a = Constant(1.0)

    masked = MaskedState(state, 2)
    uuid = masked.uuid("a")

    state.a = Constant(2.0)
    self.assertFalse(uuid == masked.uuid("a"))

  def test_cut(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'domain': (1,2,3)})
    state.a = Constant(1.0)

    masked = MaskedState(state, 2)
    a = state.interpolate(Constant(2.0))
    masked.cut(a)

  def test_function_setter(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'domain': (1,2,3)})
    state.a = Constant(1.0)
    self.assertAlmostEqual(8.0, assemble(state.a*state.dx()))

    masked = MaskedState(state, 2)
    masked.a = masked.interpolate(Constant(2.0))
    self.assertAlmostEqual(4.0, assemble(masked.a*masked.dx()))
    self.assertAlmostEqual(11.0, assemble(state.a*state.dx()))
    self.assertEqual("a", masked.a.name())
    self.assertEqual("a", state.a.name())

  def test_time(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'domain': (1,2,3)})
    state.t = 1.0
    masked = MaskedState(state, 2)
    self.assertAlmostEqual(1.0, state.t)
    self.assertAlmostEqual(1.0, masked.t)
    masked.t = 2.0
    self.assertAlmostEqual(2.0, state.t)
    self.assertAlmostEqual(2.0, masked.t)

  def test_material(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'1': 1, '2': 2, '3': 3})
    state.material['1'] = Material(a = 1.0)
    state.material['2'] = Material(a = 2.0)
    state.material['3'] = Material(a = 3.0)

    masked = MaskedState(state, 2)
    self.assertAlmostEqual(4.0,  assemble(masked.material.a*masked.dx()))
    self.assertAlmostEqual(14.0, assemble(state.material.a * state.dx()))

  def test_forbird_material_setter(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'1': 1, '2': 2, '3': 3})
    masked = MaskedState(state, 2)
    with self.assertRaises(AttributeError):
      masked.material = Material(a = 1.0)
    with self.assertRaises(AttributeError):
      masked.material['1'] = Material(a = 1.0)

  def test_normalize(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'1': 1, '2': 2, '3': 3})
    v = TestFunction(state.FunctionSpace())
    state.a = Constant((2.0, 0.0, 0.0))

    norm = assemble(inner(state.a,state.a)*v*dP)
    self.assertAlmostEqual(4.0, norm.min())
    self.assertAlmostEqual(4.0, norm.max())

    masked = MaskedState(state, 2)
    masked.a.normalize()

    norm = assemble(inner(state.a,state.a)*v*dP)
    self.assertAlmostEqual(1.0, norm.min())
    self.assertAlmostEqual(4.0, norm.max())
    self.assertAlmostEqual(2187.0, norm.sum())

  def test_full(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'1': 1, '2': 2, '3': 3})
    masked = MaskedState(state, 2)
    self.assertEqual(masked.full, state)
    self.assertEqual(state.full, state)

  def test_dP(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'1': 1, '2': 2, '3': 3})
    masked = MaskedState(state, 2)
    a = assemble(Constant(1.0)*state.dP(2))
    b = assemble(Constant(1.0)*masked.dP())
    self.assertAlmostEqual(a, b)

  def test_M_inv_diag(self):
    mesh = self.mesh_with_subdomains()
    state = State(mesh, celldomains={'1': 1, '2': 2, '3': 3})
    masked = MaskedState(state, 2)

    A = masked.M_inv_diag()
    self.assertEqual(A.size(0), masked.mesh.size(0) * 3)

  def mesh_with_subdomains(self):
    return Mesh("mesh/masked_state.xml.gz")

if __name__ == '__main__':
    unittest.main()
