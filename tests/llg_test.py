import unittest
from dolfin import *
from magnumfe import *
import os

set_log_active(False)

mesh, sample_size = DemagField.create_mesh((50.0/2.0, 50.0/2.0, 3.0/2.0), (49, 49, 1), d=3)
volume   = 50.0 * 50.0 * 3.0
arg      = "sqrt((3.141592*(x[0]/1e1))*(3.141592*(x[0]/1e1)))"
m_expr   = Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0"))
ref_file = os.path.dirname(os.path.realpath(__file__)) + "/ref/llg_v.xml"

class LlgTest(unittest.TestCase):

  def test_llg_alouges_lagrange(self):
    parameters["linear_algebra_backend"] = "Epetra"
    state = State(mesh, material = Material.py(), m = m_expr)

    llg = LLGAlougesLagrange([], scale = 1e-9)
    v   = llg.calculate_v(state, 1e-12)

    ref = Function(state.VectorFunctionSpace(), ref_file)

    error = assemble(inner(ref - v, ref - v) / inner(ref, ref) * state.dx('magnetic')) / volume
    self.assertTrue(error < 0.1)

  def test_llg_alouges_project(self):
    parameters["linear_algebra_backend"] = "PETSc"
    state = State(mesh, material = Material.py(), m = m_expr)

    llg = LLGAlougesProject([], scale = 1e-9)
    v   = llg.calculate_v(state, 1e-12)

    ref = Function(state.VectorFunctionSpace(), ref_file)

    error = assemble(inner(ref - v, ref - v) / inner(ref, ref) * state.dx('magnetic')) / volume
    self.assertTrue(error < 0.1)

if __name__ == '__main__':
    unittest.main()
