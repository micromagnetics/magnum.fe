"""
Current driven domain-wall motion with constant current and spin accumulation.
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
# Last modified by Claas Abert, 2014-10-02

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
  DemagField("FK"),
  SpinCurrent()
], scale = 1e-9)

spindiff = SpinDiffusion(scale = 1e-9)

# relax
for j in range(200):
  llg.step(state, 1e-12)
  spindiff.step(state, 1e-12)

# apply current
state.j = Constant((3e12, 0, 0))

for j in range(1000):
  llg.step(state, 1e-12)
  spindiff.step(state, 1e-12)

  if j % 10 == 0:
    f = File("data/m_%d.pvd" % (j / 10))
    f << state.m

    f = File("data/s_%d.pvd" % (j / 10))
    f << state.s
