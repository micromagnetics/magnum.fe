// Copyright (C) 2011-2015 Claas Abert
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
// Last modified by Claas Abert, 2015-01-05

#ifndef _SUB_MESH_INTERPOLATOR_H_
#define _SUB_MESH_INTERPOLATOR_H_

#include <dolfin.h>
#include <vector>

namespace magnumfe {

  /// The class offers methods for the fast mapping from functions
  /// defined on a mesh to a function defined on a submesh. The
  /// dof mapping is computed once and then cashed.
  /// 
  /// This works only for the same function spaces on a mesh/submesh
  /// pair.

  class SubMeshInterpolator
  {
  public:
    /// Create a _SubMeshInterpolator_.
    ///
    /// *Arguments*
    ///     superspace (dolfin::FunctionSpace)
    ///         The function space on the super mesh
    ///     subspace (dolfin::FunctionSpace)
    ///         The function space on the sub mesh
    SubMeshInterpolator(dolfin::FunctionSpace& superspace,
                        dolfin::FunctionSpace& subspace);

    /// Takes a vector of a super space function and returns a vector
    /// of the sub space with values cut.
    ///
    /// *Arguments*
    ///     src (dolfin::GenericVector)
    ///         The source vector from super space.
    ///     target (dolfin::GenericVector)
    ///         The taraget vector from sub space.
    void cut(const dolfin::GenericVector& src,
             dolfin::GenericVector& target);

    /// Takes a vector of a sub space function and returns a vector
    /// of the super space with unknown values set to zero.
    ///
    /// *Arguments*
    ///     src (dolfin::GenericVector)
    ///         The source vector from sub space.
    ///     target (dolfin::GenericVector)
    ///         The taraget vector from super space.
    void expand(const dolfin::GenericVector& src,
                dolfin::GenericVector& target);

  private:

    // index mapping beween sub and super space
    std::vector<int> mapping;
  };

}

#endif
