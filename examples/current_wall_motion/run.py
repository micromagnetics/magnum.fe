"""
Current driven domain-wall motion with constant current and spin accumulation.
"""

# Copyright (C) 2011-2015 Claas Abert
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
# Last modified by Claas Abert, 2015-01-05

from dolfin import *
from magnumfe import *

#######################################
#### DEFINE MESH, STATE AND MATERIAL
#######################################
mesh = BoxMesh(-600.0/2, -100.0/2, -10.0/2, 600.0/2, 100.0/2, 10.0/2, 120, 20, 1)

state = State(mesh, scale = 1e-9,
    material = Material(
      alpha      = 0.1,
      ms         = 8e5,
      Aex        = 1.3e-11,
      D0         = 1e-3,
      beta       = 0.9,
      beta_prime = 0.8,
      lambda_sf  = 10e-9,
      lambda_j   = 4e-9,
      c          = 3.125e-3
    ),
    m = Expression(('1.0 - 2*(x[0] < 0.0)', 'x[0] > -10.0 && x[0] < 10.0', '0.0')),
    s = Constant((0.0, 0.0, 0.0)),
    j = Constant((0.0, 0.0, 0.0))
)

# normalize since initial configuration is not normalized
state.m.normalize()

# setup integrators
llg = LLGAlougesProject([
  ExchangeField(),
  DemagField("FK"),
  SpinCurrent()
])

spindiff = SpinDiffusion()

# relax
for j in range(200): state.step(llg, 1e-12)

# apply constant current
state.j = Constant((3e12, 0, 0))

for j in range(1000):

  # save fields every 10th step
  if j % 10 == 0:
    f = File("data/m_%d.pvd" % (j / 10))
    f << state.m

    f = File("data/s_%d.pvd" % (j / 10))
    f << state.s

  # calculate next step
  state.step([llg, spindiff], 1e-12)
