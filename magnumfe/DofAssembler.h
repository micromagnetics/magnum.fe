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

#ifndef _DOF_ASSEMBLER_H_
#define _DOF_ASSEMBLER_H_

#include <dolfin.h>
#include "DofForm.h"

namespace magnumfe {

  /// This class contains static methods for the assembly of
  /// instances of _DofForm_. Dof forms are form with arbitrary cell
  /// sparsity patterns and are usually used to represent forms
  /// which are defined on nodes (and not as intergrals over cells).

  class DofAssembler {
  public:

    /// Assemble a _DofForm_
    ///
    /// *Arguments*
    ///     A (dolfin::GenericTensor&)
    ///         The resulting tensor.
    ///     a (_DofForm_)
    ///         The _DofForm_ to be assembled.
    ///     reset_sparsity (bool)
    ///         If set to true, the sparsity pattern of
    ///         A is reset.
    static void assemble(dolfin::GenericTensor& A,
        const DofForm& a,
        bool reset_sparsity=true);

    /// Initialize the sparsity pattern of a tensor for a _DofForm_.
    ///
    /// *Arguments*
    ///     tensor_layout (dolfin::TensorLayout&)
    ///         The layout whose sparsity pattern is reset.
    ///     a (_DofForm_)
    ///         The _DofForm_.
    static void init_tensor_layout(dolfin::TensorLayout& tensor_layout,
        const DofForm& a);
  };
}
#endif
