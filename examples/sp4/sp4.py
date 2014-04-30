from dolfin import *
from magnumfe import *

#######################################
#### GENERATE MESH WITH SHELL
#######################################

mesh, sample_size = DemagField.create_mesh((500.0/2.0, 125.0/2.0, 3.0/2.0), (100, 25, 1), d=6)

#######################################
#### RELAX SYSTEM TO S-STATE
#######################################

# define start magnetization
arg     = "sqrt((3.141592*(x[0]/1e3))*(3.141592*(x[0]/1e3)))"
m_start = Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0"))

state   = State(mesh, material = Material.py(), m = m_start)
llg     = LLG([DemagField(sample_size, 2)], scale = 1e-9)

state.material.alpha = 1.0
for i in range(200): llg.step(state, 1e-11)

#######################################
#### SIMULATE SWITCHING
#######################################

state.material.alpha = 0.02

llg = LLG([
    ExternalField((-24.6e-3/Constants.mu0, +4.3e-3/Constants.mu0, 0.0)),
    DemagField(sample_size, 2)
], scale = 1e-9)

logfile = open("sp4.dat", "w", 0)
dt, T = 2e-13, 1e-9

t  = 0.0
for i in range(int(T / dt)):
  t = i * dt
  
  #if (i % nth_step == 10):
  #  f = File("data/m_%d.pvd" % int(i/nth_step))
  #  f << state.m

  # write scalar information
  volume = 187500
  m_x = assemble(state.m[0] / volume * dx)
  m_y = assemble(state.m[1] / volume * dx)
  m_z = assemble(state.m[2] / volume * dx)
  scalar_file.write("%.10f %f %f %f\n" % (t*1e9, m_x, m_y, m_z))

  # calculate next step
  llg.step(state, dt)

logfile.close()
