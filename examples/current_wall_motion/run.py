from dolfin import *
from magnumfe import *
from math import sqrt

#######################################
#### GENERATE MESH WITH SUBDOMAINS
#######################################

class FMB(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 200.0) or near(x[0], 200.0)

mesher = Mesher()
mesher.create_cuboid((600.0/2.0, 100.0/2.0, 10.0/2.0), (100, 25, 1))
#mesher.create_celldomain(FM1(), 1)
mesher.create_facetdomain(FMB(), 2)
mesh = mesher.mesh()

#######################################
#### DEFINE STATE AND MATERIAL
#######################################

state = State(mesh,
    celldomains  = {'magnetic': 1, 'conducting': 1},
    facetdomains = {'outermagnet': 2},
    m = Expression(('1.0 - 2*(x[0] < 0.0)', '0.0', '0.0')),
    s = Constant((0, 0, 0)),
    j = Constant((0, 0, 0))
    #j = Constant((1e11, 0, 0))
  )

state.material['all'] = Material(
  alpha      = 1.0,
  ms         = 8e5,
  Aex        = 1.3e-11,
  D0         = 1e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 4e-9,
  c          = 3.125e-3
)

llg = LLGAlougesProject([
  DemagFieldFK(),
  SpinCurrent()
], scale = 1e-9)

spindiff = SpinDiffusion(scale = 1e-9)

v_py = Constant(assemble(Constant(1.0)*state.dx('all')))

# relax
for j in range(200):
  llg.step(state, 1e-12)
  spindiff.step(state, 1e-12)

# apply current
state.j = Constant((5e12, 0, 0))

for j in range(1000):
  llg.step(state, 1e-12)
  spindiff.step(state, 1e-12)

  if j % 10 == 0:
    f = File("data/m_%d.pvd" % (j / 10))
    f << state.m

    f = File("data/s_%d.pvd" % (j / 10))
    f << state.s
