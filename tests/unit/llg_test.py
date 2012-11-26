import unittest
from dolfin import *
from magnumfe import *
import numpy

class LlgTest(unittest.TestCase):

  def test_integrate(self):
    mesh  = DemagField.create_mesh((500.0/2.0, 125.0/2.0, 3.0/2.0), (80, 20, 1), d=2, scale=1e-9)
    VV    = VectorFunctionSpace(mesh, "CG", 1, 3)

    arg = "sqrt((3.141592*(x[0]/1e-6))*(3.141592*(x[0]/1e-6)))"
    m   = interpolate(Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0")), VV)

    llg = LLG2(mesh, Material.py())

    for i in range(1000):
      m = llg.step(m, 1e-10)
      if (i % 10 == 0):
        f = File("data/m_%d.pvd" % i)
        f << m

if __name__ == '__main__':
    unittest.main()
