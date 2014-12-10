import unittest
from dolfin import *
from magnumfe import *
import os

set_log_active(False)

mesh, sample_size = DemagField.create_mesh((50.0/2.0, 50.0/2.0, 3.0/2.0), (49, 49, 1), d=3)
V        = VectorFunctionSpace(mesh, "CG", 1)
ref      = Function(V, os.path.dirname(os.path.realpath(__file__)) + "/ref/llg_v.xml")
volume   = 50.0 * 50.0 * 3.0
arg      = "sqrt((3.141592*(x[0]/1e1))*(3.141592*(x[0]/1e1)))"
m_expr   = Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0"))
demag    = DemagField("ST", sample_size, 1)

class LlgTest(unittest.TestCase):

<<<<<<< HEAD
  # TODO test different backends explicitly
=======
>>>>>>> d9a478c... remove Epetra code from LLG integrators
  def test_llg_alouges_lagrange(self):
    state = State(mesh, material = Material.py(), scale = 1e-9, m = m_expr)

    llg = LLGAlougesLagrange([demag, ExchangeField()])
    v   = llg.calculate_v(state, 1e-12)

    error = assemble(inner(ref - v, ref - v) / inner(ref, ref) * state.dx('magnetic')) / volume * 1e-12
    self.assertTrue(error < 1e-13)

  def test_llg_alouges_project(self):
    state = State(mesh, material = Material.py(), scale = 1e-9, m = m_expr)

    llg = LLGAlougesProject([demag, ExchangeField()])
    v   = llg.calculate_v(state, 1e-12)

    error = assemble(inner(ref - v, ref - v) / inner(ref, ref) * state.dx('magnetic')) / volume * 1e-12
    self.assertTrue(error < 1e-13)

  def test_llg_alouges_full_implicit(self):
    state = State(mesh.with_shell, material = Material.py(), scale = 1e-9, m = m_expr)

    llg = LLGAlougesFullImplicit([], sample_size, demag_order = 1)
    v   = llg.calculate_v(state, 1e-12)

    error = assemble(inner(ref - v, ref - v) / inner(ref, ref) * dx(mesh)) / volume * 1e-12
    self.assertTrue(error < 1e-8)

if __name__ == '__main__':
    unittest.main()
