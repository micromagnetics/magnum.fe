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
    mesh  = DemagField.create_mesh((500.0/2.0, 125.0/2.0, 3.0/2.0), (80, 20, 1), d=2)
    VV    = VectorFunctionSpace(mesh, "CG", 1, 3)

    arg = "sqrt((3.141592*(x[0]/1000))*(3.141592*(x[0]/1000)))"
    m   = interpolate(Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0")), VV)

    #f = File("data/start.pvd")
    #f << m

    #m = interpolate(Constant((0.95, 0.3122, 0)), VV)
    #m = Constant((0.9, 0.1, 0))

    llg = LLG2(mesh, Material.py())

    for i in range(100):
      m = llg.step(m, 1e-11)
      f = File("data/m_%d.pvd" % i)
      f << m

if __name__ == '__main__':
    unittest.main()
