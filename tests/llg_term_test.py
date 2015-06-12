import unittest
from magnumfe import *

set_log_active(False)

class LLGTermTest(unittest.TestCase):

  def test_external_field(self):
    mesh = UnitCubeMesh(2, 2, 2)
    state = State(mesh)
    extfield = ExternalField(Constant((1.0, 0.0, 0.0)))
    f = extfield.field(state)

    self.assertAlmostEqual(1.0, assemble(inner(Constant((1.0, 0.0, 0.0)), f) * dx))
    self.assertAlmostEqual(0.0, assemble(inner(Constant((0.0, 1.0, 0.0)), f) * dx))

  def test_time_varying_external_field(self):
    mesh = UnitCubeMesh(2, 2, 2)
    state = State(mesh)

    extfield = ExternalField(lambda state: Constant((state.t, 0.0, 0.0)))
    state.t = 0.0
    self.assertAlmostEqual(0.0, assemble(inner(Constant((1.0, 0.0, 0.0)), extfield.field(state)) * state.dx()))
    state.t = 1.0
    self.assertAlmostEqual(1.0, assemble(inner(Constant((1.0, 0.0, 0.0)), extfield.field(state)) * state.dx()))

  def test_exchange_field(self):
    mesh = UnitCubeMesh(5, 5, 5)
    state = State(mesh, material = Material(ms = 1.0, Aex = 1.0), m = self.m_test(), scale = 1e-9)
    field = ExchangeField()
    f = field.field(state)
    self.assertAlmostEqual(-3.769911732099959e+25, assemble(inner(state.m, f) * dx))

  def test_exchange_field(self):
    mesh = UnitCubeMesh(5, 5, 5)
    state = State(mesh, material = Material(ms = 1.0, K_uni = 1.0, K_uni_axis = (1, 0, 0)), m = self.m_test(), scale = 1e-9)
    field = UniaxialAnisotropyField()
    f = field.field(state)
    self.assertAlmostEqual(479218.7944399755, assemble(inner(state.m, f) * dx))

  def m_test(sefl):
    return Expression((
      "cos(2*pi*x[0])",
      "sin(2*pi*x[0])",
      "0.0"
    ))

  def error(self, s):
    return assemble(inner(ref-s, ref-s) / inner(ref, ref) * state.dx('magnetic')) / volume

if __name__ == '__main__':
    unittest.main()
