"""
Switching of a permalloy multilayer structure with an electric current.
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
# Last modified by Claas Abert, 2014-12-11

from dolfin import *
from magnumfe import *

mesher = Mesher()
mesher.read_file("multilayer.msh")
mesh = mesher.mesh()

#######################################
#### DEFINE STATE AND MATERIAL
#######################################

state = State(mesh,
    celldomains  = {'magnetic': (1, 3), 'conducting': (1, 2, 3)},
    facetdomains = {'outermagnet': 1, 'interface': 2},
    scale = 1e-9,
    m = Constant((1, 0, 0)),
    s = Constant((0, 0, 0)),
    j = Constant((0, 0, -1e12))
  )

state.material['magnetic'] = Material(
  alpha      = 1.0,
  ms         = 8e5,
  Aex        = 1.3e-11,
  D0         = 1e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 2.236e-09,
  c          = 3.155e-3
)

state.material['!magnetic']  = Material(
  D0         = 5e-3,
  beta       = 0.9,
  beta_prime = 0.8,
  lambda_sf  = 10e-9,
  lambda_j   = 2.236e-09,
  c          = 3.155e-3
)

state.m.assign(state.interpolate({1: Constant((1.0, 0.0, 0.0)), 3: Constant((-1.0, 0.0, 0.0))}))
state.m.normalize()

llg      = LLGAlougesProject([ExchangeField(), DemagField("FK"), SpinCurrent()])
spindiff = SpinDiffusion()

for i in range(300):

  # save magnetization
  f = File("data/m_%d.pvd" % i)
  f << state.m.crop('magnetic')

  state.step([llg, spindiff], 1e-11)
