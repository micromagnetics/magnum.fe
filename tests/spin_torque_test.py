import unittest
from magnumfe import *
import numpy

set_log_active(False)

mesh = Mesh("mesh/spin_torque.xml.gz")

permalloy = Material(
  alpha      = 0.1,
  ms         = 8e5,
  Aex        = 1.3e-11,
  D0         = 1e-3,
  C0         = -7.5e27 * Constants.e * 1e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 2.236068e-09,
  c          = 3.157346e-03
)

# vortex parametrization
norm    = "sqrt(x[0]*x[0] + x[1]*x[1] + 10*10)"
m_start = Expression(("-x[1]/r", "x[0]/r", "10/r"), r = Expression(norm))

state = State(mesh,
  facetdomains = {'outermagnet': (1,2,3), 'left': 2, 'right': 3},
  material = permalloy,
  scale = 1e-9,
  m = m_start
)

volume = state.volume('all')
ref    = Function(state.VectorFunctionSpace(), "ref/s.xml")

class SpinTorqueTest(unittest.TestCase):

  def test_spin_torque_integration(self):
    state.m   = m_start
    state.j   = Constant((1e12, 0.0, 0.0))
    state.s   = Constant((0.0, 0.0, 0.0))
    spin_diff = SpinDiffusion()

    for i in range(10): spin_diff.step(state, 1e-12)
    self.assertTrue(self.error(state.s) < 1e-14)

  def error(self, s):
    return assemble(inner(ref-s, ref-s) / inner(ref, ref) * state.dx('magnetic')) / volume

if __name__ == '__main__':
    unittest.main()
