"""
MuMag Standard Problem #4 computed with Hybrid FEM/BEM  method by 
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
# Last modified by Claas Abert, 2014-06-18

from dolfin import *
from magnumfe import *

#######################################
#### GENERATE MESH WITH SHELL
#######################################

mesh = BoxMesh(-500.0/2, -125.0/2, -3.0/2, 500.0/2, 125.0/2, 3.0/2, 100, 25, 1)
volume = assemble(Constant(1.0)*dx(mesh))

#######################################
#### RELAX SYSTEM TO S-STATE
#######################################

# define start magnetization
arg     = "sqrt((3.141592*(x[0]/1e3))*(3.141592*(x[0]/1e3)))"
m_start = Expression(("cos(%s)" % arg, "sin(%s)" % arg, "0.0"))

state   = State(mesh, material = Material.py(), m = m_start)
llg     = LLGAlougesProject([DemagField("FK")], scale = 1e-9)

state.material.alpha = 1.0
for i in range(200): llg.step(state, 2e-11)

#######################################
#### SIMULATE SWITCHING
#######################################

state.material.alpha = 0.02

llg = LLGAlougesProject([
    ExternalField((-24.6e-3/Constants.mu0, +4.3e-3/Constants.mu0, 0.0)),
    DemagField("FK")
], scale = 1e-9)

logfile = open("sp4_fk.dat", "w", 0)
dt, T = 2e-13, 1e-9

t  = 0.0
for i in range(int(T / dt)):
  t = i * dt
  
  #if (i % 10 == 0):
  #  f = File("data/m_%d.pvd" % int(i/10))
  #  f << state.m

  # write scalar information
  logfile.write("%.10f %f %f %f\n" % ((t*1e9,) + state.m.average()))

  # calculate next step
  llg.step(state, dt)

logfile.close()
