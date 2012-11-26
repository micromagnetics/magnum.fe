import unittest
from dolfin import *
from magnumfe import *
import numpy

class LlgTest(unittest.TestCase):

  def setUp(self):
    pass

#  def test_llg(self):
#    mesh  = DemagField.create_mesh((500, 125, 3), (100, 25, 3), d=3)
#    #mesh  = DemagField.create_mesh((500, 125, 3), (40, 10, 1), d=2)
#    VV = VectorFunctionSpace(mesh, "CG", 1, 3)
#    VS = FunctionSpace(mesh, "CG", 1)
#
#    m0 = interpolate(Constant((1.0, 0.01, 0.0)), VV)
#
#    demag   = DemagField(mesh, order=2)
#    u_demag = demag.calculate(m0)
#
#    f = File("u.pvd")
#    f << u_demag
#
#    llg     = LLG(mesh, Material.py())
#    dm      = llg.calculate(m0, 1e-9, u_demag)
#
#    f = File("dm.pvd")
#    f << dm
  def test_integrate(self):
    #mesh  = DemagField.create_mesh((500, 125, 3), (100, 25, 3), d=3)
    mesh  = DemagField.create_mesh((500, 125, 3), (40, 10, 1), d=2)
    VV    = VectorFunctionSpace(mesh, "CG", 1, 3)
    m0    = interpolate(Constant((1.0, 0.01, 0.0)), VV)

    llg = LLG(mesh, Material.py())
    m = llg.integrate(m0, 3e-9)

    f = File("m.pvd")
    f << m


if __name__ == '__main__':
    unittest.main()
