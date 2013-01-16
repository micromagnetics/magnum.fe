import unittest
from dolfin import *
from magnumfe import *
import numpy

mesh           = DemagField.create_mesh((50.0/2.0, 50.0/2.0, 3.0/2.0), (10, 10, 1), d=3)
VV             = VectorFunctionSpace(mesh, "CG", 1, 3)

material       = Material.py()
material.alpha = 1.0

arg            = "sqrt((3.141592*(x[0]/1e1))*(3.141592*(x[0]/1e1)))"
m_expr         = Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0"))
dt             = 1e-12

class LlgTest(unittest.TestCase):

  def test_llg2(self):
    llg = LLG2(mesh, material, scale=1e-9, demag_order=1)
    m   = llg.interpolate(m_expr)
    dm  = llg.calculate_dm(m, dt)

    f  = File("data3/m.pvd")
    f << m

    f  = File("data3/dm_a.pvd")
    f << dm

  def test_llg3(self):
    llg = LLG3(mesh, material, scale=1e-9, demag_order=1)
    m   = llg.interpolate(m_expr)
    dm  = llg.calculate_dm(m, dt)

    f  = File("data3/dm_b.pvd")
    f << dm

if __name__ == '__main__':
    unittest.main()
