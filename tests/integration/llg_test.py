import unittest
from dolfin import *
from magnumfe import *
import numpy

mesh  = DemagField.create_mesh((500.0/2.0, 125.0/2.0, 3.0/2.0), (100, 25, 1), d=4)
VV    = VectorFunctionSpace(mesh, "CG", 1, 3)

class LlgTest(unittest.TestCase):

  def prepare_s_state(self):
    #mesh  = DemagField.create_mesh((500.0/2.0, 125.0/2.0, 3.0/2.0), (100, 25, 1), d=2)
    #VV    = VectorFunctionSpace(mesh, "CG", 1, 3)

    filename = "data/s-state.xml"
    try:
      with open(filename) as f: pass
      # File exists, read and return
      return Function(VV, filename)

    except IOError as e:
      material = Material.py()
      material.alpha = 1.0
      llg = LLG(mesh, material, scale=1e-9, demag_order=2)

      # File does not exists, calculate s-state and return
      arg = "sqrt((3.141592*(x[0]/1e3))*(3.141592*(x[0]/1e3)))"
      m   = llg.interpolate(Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0")))


      volume = 187500
      for i in range(400):
        m = llg.step(m, 1e-11)
        m_x = assemble(m[0] / volume * dx)
        print "Mx: %f" % m_x


      f = File(filename)
      f << m

      return m

  def test_sp4(self):
    #mesh  = DemagField.create_mesh((500.0/2.0, 125.0/2.0, 3.0/2.0), (100, 25, 1), d=2)
    #VV    = VectorFunctionSpace(mesh, "CG", 1, 3)

    m = self.prepare_s_state()

    llg = LLG(mesh, Material.py(), scale=1e-9, demag_order=2)
    field = Constant((-24.6e-3/Constants.mu0, +4.3e-3/Constants.mu0, 0.0))

    scalar_file = open("data/sp4.dat","w",0)
    dt = 1e-13
    T  = 1e-9

    t  = 0.0
    print "Total Steps: %d" % int(T / dt)
    for i in range(int(T / dt)):
      t = i * dt
      
      # write magnetization configuration (only each 20th step)
      nth_step = 10 
      if (i % nth_step == 0):
        f = File("data/m_%d.pvd" % int(i/nth_step))
        f << m

      # write scalar information
      volume = 187500
      m_x = assemble(m[0] / volume * dx)
      m_y = assemble(m[1] / volume * dx)
      m_z = assemble(m[2] / volume * dx)
      scalar_file.write("%.10f %f %f %f\n" % (t*1e9, m_x, m_y, m_z))

      # calculate next step
      m = llg.step(m, dt, h_ext = field)

    scalar_file.close()

if __name__ == '__main__':
    unittest.main()
