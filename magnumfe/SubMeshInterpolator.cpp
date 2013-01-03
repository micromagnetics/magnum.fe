// Copyright (C) 2011-2012 Claas Abert
//
// This file is part of magnum.fe.
//
// magnum.fe is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// magnum.fe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with magnum.fe. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2012-11-30
// Last changed: 2012-11-30

#include "SubMeshInterpolator.h"
#include <iostream>

using namespace magnumfe;

SubMeshInterpolator::SubMeshInterpolator(dolfin::FunctionSpace& Vsuper, dolfin::FunctionSpace& Vsub): mapping(Vsub.dim())
{
  // TODO assert that Vsuper and Vsub use the same elements
  // TODO assert that mesh of Vsub is submesh of Vsuper mesh

  // do it right
  dolfin::Function fsuper(Vsuper);
  dolfin::Function fsub(Vsub);

  std::pair<size_t, size_t> local_range = Vsuper.dofmap()->ownership_range();
  for (size_t i = local_range.first; i < local_range.second; ++i) {
    const double value = i;
    const dolfin::la_index idx = i;
    fsuper.vector()->set(&value, 1, &idx);
  }
  fsuper.vector()->apply("insert");

  Vsub.interpolate(*fsub.vector(), fsuper);

  // TODO MPIifize the following
  for (size_t i = 0; i<mapping.size(); ++i) {
    const dolfin::la_index idx = i;
    mapping[i] = floor((*fsub.vector())[idx] + 0.5);
  }
}
//-----------------------------------------------------------------------------
void SubMeshInterpolator::cut(const dolfin::GenericVector& src, dolfin::GenericVector& target) {
  for (size_t i=0; i<mapping.size(); ++i) {
    const double value = src[mapping[i]];
    const dolfin::la_index idx = i;
    target.set(&value, 1, &idx);
  }
  target.apply("insert");
}
//-----------------------------------------------------------------------------
void SubMeshInterpolator::expand(const dolfin::GenericVector& src, dolfin::GenericVector& target) {
  target.zero();
  for (uint i=0; i<mapping.size(); ++i) {
    const double value = src[i];
    const dolfin::la_index j = mapping[i];
    target.set(&value, 1, &j);
  }
  target.apply("insert");
}
