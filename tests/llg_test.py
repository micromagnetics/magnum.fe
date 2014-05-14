import unittest
from dolfin import *
from magnumfe import *
import os

set_log_active(False)

mesh, sample_size = DemagField.create_mesh((50.0/2.0, 50.0/2.0, 3.0/2.0), (50, 50, 1), d=3)
volume   = 50.0 * 50.0 * 3.0
arg      = "sqrt((3.141592*(x[0]/1e1))*(3.141592*(x[0]/1e1)))"
m_expr   = Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0"))
ref_file = os.path.dirname(os.path.realpath(__file__)) + "/ref/llg_v.xml"

class LlgTest(unittest.TestCase):

  def test_llg(self):
    state = State(mesh, material = Material.py(), m = m_expr)

    llg = LLG([DemagField(sample_size, 2)], scale = 1e-9)
    v   = llg.calculate_v(state, 1e-12)

    ref = Function(state.VectorFunctionSpace(), ref_file)

    error = assemble(inner(ref - v, ref - v) / inner(ref, ref) * state.dx('magnetic')) / volume
    self.assertTrue(error < 0.3)

  # TODO currently broken (segfault probably caused by cbc.block)
  def ttest_llg4(self):
    llg = LLG4(mesh, material, scale=1e-9, demag_order=1)
    m   = llg.interpolate(m_expr)
    dm  = llg.calculate_dm(m, 1e-12)

    ref = Function(VV, ref_file)

    f = File("ref/dm.pvd")
    f << dm

    error = assemble(inner(ref - dm, ref - dm) / inner(ref, ref) * dx(mesh)) / volume
    print error
    self.assertTrue(error < 0.3)

if __name__ == '__main__':
    unittest.main()
