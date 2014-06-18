"""
Current induced switching of a soft magnetic layer in a multilayer structure.
"""

# Copyright (C) 2011-2014 Claas Abert
#
# This file is part of magnum.fe. 
#
# magnum.fe is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# magnum.fe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with magnum.fe. If not, see <http://www.gnu.org/licenses/>.
# 
# Last modified by Claas Abert, 2014-06-18

from dolfin import *
from magnumfe import *
from math import sqrt

#######################################
#### GENERATE MESH WITH SUBDOMAINS
#######################################

class FM1(SubDomain):
  def inside(self, x, on_boundary):
    return between(x[2], (-12.5, -2.5))

class NM(SubDomain):
  def inside(self, x, on_boundary):
    return between(x[2], (-2.5, 2.5))

class FM2(SubDomain):
  def inside(self, x, on_boundary):
    return between(x[2], (2.5, 12.5))

class FMB(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[2], 12.5) or near(x[2], -12.5)

mesher = Mesher()
mesher.create_cuboid((100.0/2.0, 100.0/2.0, 25.0/2.0), (25, 25, 10))
mesher.create_celldomain(FM1(), 1)
mesher.create_celldomain(NM(),  2)
mesher.create_celldomain(FM2(), 3)
mesher.create_facetdomain(FMB(), 2)
mesh = mesher.mesh()

#######################################
#### DEFINE STATE AND MATERIAL
#######################################

state = State(mesh,
    celldomains  = {'pinned': 1, 'py': 3, 'magnetic': (1, 3), 'conducting': (1, 2, 3)},
    facetdomains = {'outermagnet': 2},
    m = Expression(('1.0 - 2*(x[2] < 0.0)', '0.0', '0.0')),
    s = Constant((0, 0, 0)),
    j = Constant((0, 0, -5e12))
  )

state.material['pinned'] = Material(
  alpha      = 1.0,
  ms         = 1.43/Constants.mu0,
  Aex        = 2.158e-11,
  #K_uni      = 6.6e6,     # gyro freq: 3.25*10^11
  K_uni      = 5e5,
  K_uni_axis = (1, 0, 0),
  D0         = 1e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 4e-9,
  c          = 3.125e-3
)

state.material['py'] = Material(
  alpha      = 1.0,
  ms         = 8e5,
  Aex        = 1.3e-11,
  K_uni      = 0.0,
  K_uni_axis = (1, 0, 0),
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

llg = LLGAlougesProject([
  DemagField("FK"),
  UniaxialAnisotropyField(),
  SpinCurrent()
], scale = 1e-9)

spindiff = SpinDiffusion(scale = 1e-9)

v_py = Constant(assemble(Constant(1.0)*state.dx('py')))

# relax
for j in range(5000):
  llg.step(state, 1e-12)
  spindiff.step(state, 1e-12)

  if j % 10 == 0:
    f = File("data/m_%d.pvd" % (j / 10))
    f << state.m.crop('magnetic')

    f = File("data/s_%d.pvd" % (j / 10))
    f << state.s
