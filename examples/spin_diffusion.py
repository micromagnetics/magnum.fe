from dolfin import *
from magnumfe import *

#######################################
#### GENERATE MESH WITH SHELL
#######################################

class FM1(SubDomain):
  def inside(self, x, on_boundary):
    return between(x[2], (-90, 10))

class NM(SubDomain):
  def inside(self, x, on_boundary):
    return between(x[2], (10, 30))

class FM2(SubDomain):
  def inside(self, x, on_boundary):
    return between(x[2], (30, 90))

mesh, sample_size = DemagField.create_mesh(
    (128.0/2.0, 64.0/2.0, 180/2.0), (30, 15, 36),
    d = 4,
    domains = { 1: FM1(), 2: NM(), 3: FM2()}
  )


#######################################
#### RELAX SYSTEM TO S-STATE
#######################################

# define start magnetization
state   = State(mesh, {'magnetic': (1, 3), 'conducting': (1, 2, 3)}, m = (0, 0, 1), s = (0, 0, 0), j = (0, 0, 1e11))

state.material['magnetic']  = Material({
  alpha      = 0.1,
  ms         = 8e5,
  Aex        = 13e-12,
  K          = 500,
  K_axis     = (0, 0, 1),
  D0         = 1e-3,
})
state.material['!magnetic']  = Material({
  D0     = 5e-3
})
state.material['conducting'] = Material({
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 4e-9,
  c          = 3.125e-3 # J?
})

llg      = LLG([DemagField(sample_size, 2)], scale = 1e-9)
spindiff = SpinDiffusion()

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
  
  #if (i % 10 == 0):
  #  f = File("data/m_%d.pvd" % int(i/10))
  #  f << state.m

  # write scalar information
  volume = 187500
  m_x = assemble(state.m[0] / volume * dx)
  m_y = assemble(state.m[1] / volume * dx)
  m_z = assemble(state.m[2] / volume * dx)
  logfile.write("%.10f %f %f %f\n" % (t*1e9, m_x, m_y, m_z))

  # calculate next step
  llg.step(state, dt)

logfile.close()
