"""
MuMag Standard Problem #5 computed with Hybrid FEM/BEM  method by 
Fredkin and Koehler for the demagnetization-field computation.
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
# Last modified by Claas Abert, 2014-09-30

from dolfin import *
from magnumfe import *

#######################################
#### CREATE MESH AND MATERIAL
#######################################

mesh = BoxMesh(-100.0/2, -100.0/2, -10.0/2, 100.0/2, 100.0/2, 10.0/2, 40, 40, 1)

permalloy = Material(
  alpha      = 0.1,
  ms         = 8e5,
  Aex        = 1.3e-11,
  D0         = 1e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 2.236068e-09,
  c          = 3.157346e-03
)

# define start magnetization
norm    = "sqrt(x[0]*x[0] + x[1]*x[1] + 10*10)"
m_start = Expression(("-x[1]/%s" % norm, "x[0]/%s" % norm, "10/%s" % norm))
state   = State(mesh, material = permalloy, scale = 1e-9,
    m = m_start,
    j = Constant((1e12, 0.0, 0.0)),
    s = Constant((0.0, 0.0, 0.0))
)

#######################################
#### RELAX SYSTEM TO VORTEX STATE
#######################################

llg = LLGAlougesProject([ExchangeField(), DemagField("FK")])

state.material.alpha = 1.0

for i in range(400): llg.step(state, 2e-11)

#######################################
#### APPLY CURRENT
#######################################

state.material.alpha = 0.1

llg      = LLGAlougesProject([ExchangeField(), DemagField("FK"), SpinCurrent()])
spindiff = SpinDiffusion()

logfile = open("sp5.dat", "w", 0)
dt, T = 1e-13, 14e-9

for i in range(int(T / dt)):

  # write field every 100th step
  if (i % 100 == 0):
    f = File("data/m_%d.pvd" % int(i/100))
    f << state.m

  # write scalar information every 10th step
  if (i % 10 == 0):
    logfile.write("%.10f %f %f %f\n" % ((state.t*1e9,) + state.m.average()))

  # calculate next step
  state.step([llg, spindiff], dt)

logfile.close()
