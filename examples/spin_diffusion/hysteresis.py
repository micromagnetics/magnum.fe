from dolfin import *
from magnumfe import *

#######################################
#### GENERATE MESH WITH SUBDOMAINS
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

#mesh, sample_size = DemagField.create_mesh(
#    (128.0/2.0, 64.0/2.0, 180/2.0), (30, 15, 36),
#    d = 4,
#    domains = { 1: FM1(), 2: NM(), 3: FM2()}
#  )

mesher = Mesher()
mesher.create_cuboid((128.0/2.0, 64.0/2.0, 180/2.0), (30, 15, 36))
mesher.create_celldomain(FM1(), 1)
mesher.create_celldomain(NM(),  2)
mesher.create_celldomain(FM2(), 3)
mesh = mesher.mesh()

#######################################
#### DEFINE STATE AND MATERIAL
#######################################

state = State(mesh, {'magnetic': (1, 3), 'conducting': (1, 2, 3)},
    m = Constant((1, 0, 0)),
    s = Constant((0, 0, 0)),
    j = Constant((0, 0, 1e11))
  )

state.material['magnetic'] = Material(
  alpha      = 0.1,
  ms         = 8e5,
  Aex        = 13e-12,
  K_uni      = 5e2,
  K_uni_axis = (0, 0, 1),
  D0         = 1e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 4e-9,
  c          = 3.125e-3
)
state.material['!magnetic']  = Material(
  D0         = 5e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 4e-9,
  c          = 3.125e-3
)

external_field = ExternalField((0.0, 0.0, 0.0))
llg      = LLGAlougesProject([
  DemagFieldFK(),
  UniaxialAnisotropyField(),
  SpinCurrent(),
  external_field
], scale = 1e-9)

spindiff = SpinDiffusion()

v_magnetic = Constant(assemble(Constant(1.0)*state.dx('magnetic')))
logfile = open("hyst.dat", "w", 0)

steps = 100
for i in range(steps+1):
  external_field.set((0.0, 0.0, (1.0 - 2.0*i/steps) * 60e-3/Constants.mu0))

  # relax
  for j in range(2000):
    llg.step(state, 1e-12)
    spindiff.step(state, 1e-12)

    if j % 10 == 0:
      f = File("data/m_%d.pvd" % (j / 10))
      f << state.m

      f = File("data/s_%d.pvd" % (j / 10))
      f << state.s

  # XXX
  exit()

  m_x = assemble(state.m[0] / v_magnetic * state.dx('magnetic'))
  m_y = assemble(state.m[1] / v_magnetic * state.dx('magnetic'))
  m_z = assemble(state.m[2] / v_magnetic * state.dx('magnetic'))
  logfile.write("%f %f %f\n" % (m_x, m_y, m_z))

  f = File("data/m_%d.pvd" % int(i))
  f << state.m

  f = File("data/s_%d.pvd" % int(i))
  f << state.s

logfile.close()
