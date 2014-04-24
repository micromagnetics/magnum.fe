import unittest
from dolfin import *
from magnumfe import *
import numpy

mesh, sample_size  = DemagField.create_mesh((500.0/2.0, 125.0/2.0, 3.0/2.0), (100, 25, 1), d=4)
VV    = VectorFunctionSpace(mesh, "CG", 1, 3)

class LlgTest(unittest.TestCase):

  def prepare_s_state(self):
    filename = "data/s-state.xml"
    try:
      with open(filename) as f: pass
      # File exists, read and return
      return Function(VV, filename)

    except IOError as e:
      material = Material.py()
      material.alpha = 1.0

      # File does not exists, calculate s-state and return
      arg = "sqrt((3.141592*(x[0]/1e3))*(3.141592*(x[0]/1e3)))"
      m_start = Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0"))

      state = State(mesh, material = material, m = m_start)
      llg = LLG([DemagField(sample_size, 2)], scale = 1e-9)

      volume = 187500
      for i in range(400):
        llg.step(state, 1e-11)
        m_x = assemble(state.m[0] / volume * dx)
        print "Mx: %f" % m_x


      f = File(filename)
      f << state.m

      return state.m

  def test_sp4(self):
    m_start = self.prepare_s_state()
    state   = State(mesh, material = Material.py(), m = m_start)

    llg = LLG([
        ExternalField((-24.6e-3/Constants.mu0, +4.3e-3/Constants.mu0, 0.0)),
        DemagField(sample_size, 2)
    ], scale = 1e-9)

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
        f << state.m

      # write scalar information
      volume = 187500
      m_x = assemble(state.m[0] / volume * dx)
      m_y = assemble(state.m[1] / volume * dx)
      m_z = assemble(state.m[2] / volume * dx)
      scalar_file.write("%.10f %f %f %f\n" % (t*1e9, m_x, m_y, m_z))

      # calculate next step
      llg.step(state, dt)

    scalar_file.close()

if __name__ == '__main__':
    unittest.main()
