import unittest
from dolfin import *
from magnumfe import *

set_log_active(False)

class LLGTermTest(unittest.TestCase):

  def test_external_field(self):
    mesh = UnitCubeMesh(2, 2, 2)
    state = State(mesh, m = Constant((1.0, 0.0, 0.0)))
    extfield = ExternalField(Constant((1.0, 0.0, 0.0)))
    f = extfield.field(state)

    self.assertAlmostEqual(1.0, assemble(inner(state.m, f) * dx))

  # TODO test something more sophisticated than a constant
  def test_exchange_field(self):
    mesh = UnitCubeMesh(2, 2, 2)
    state = State(mesh, material = Material(ms = 1.0, Aex = 1.0), m = Constant((1.0, 0.0, 0.0)))
    exfield = ExchangeField()
    f = exfield.field(state)

    self.assertAlmostEqual(0.0, assemble(inner(state.m, f) * dx))


if __name__ == '__main__':
    unittest.main()
