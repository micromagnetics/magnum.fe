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

#ifndef _CG_DOF_FORM_H_
#define _CG_DOF_FORM_H_

#include "DofForm.h"
#include <dolfin.h>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace magnumfe {

  /// This class represents a user-defined DofForm suited for CG
  /// FunctionSpaces. The sparsity pattern is setup for node-wise
  /// assembly.
  ///
  /// An CGDofForm is defined by overloading the nodeEval() method.
  /// 
  /// The rank of the form and the number of coefficients have to
  /// be supplied to the constructor.

  class CGDofForm: public DofForm
  {
  public:

    /// Create CGDofForm
    ///
    /// *Arguments*
    ///     rank (uint)
    ///         Rank of the form.
    ///     num_coefficients (uint)
    ///         Number of coefficients of the form.
    CGDofForm(uint rank, uint num_coefficients): DofForm(rank, num_coefficients) {};

    /// Evaluate form at a single CG node.
    ///
    /// *Arguments*
    ///     A (double*)
    ///         Linearized result of the evaluation.
    ///     w (const double* double*)
    ///         Linearized values of all coefficients at CG node.
    virtual void nodeEval(double *A, const double * const *w) const = 0;

    /// Evaluate form for a given cell.
    ///
    /// *Arguments*
    ///     A (double*)
    ///         Linearized result of the evaluation.
    ///     w (const double* double*)
    ///         Linearized values of all coefficients in the cell.
    void eval(double* A, const double * const * w) const;

    /// Return cell coordinates of non-zero entries of the form.
    ///
    /// *Arguments*
    ///     entries (boost::multi_array<uint, 2>&)
    ///         Coordinates of non-zero entries of the form.
    void cell_sparsity(boost::multi_array<uint, 2>& entries) const;

    /// Number of non-zero tensor entries per cell.
    ///
    /// *Returns*
    ///     uint
    ///         Number of non-zero tensor entries.
    uint non_zero_entries() const;

    /// Number of non-zero tensor entrues per node.
    ///
    /// *Returns*
    ///     uint
    ///         Number of non-zero tensor entries per node.
    uint entries_per_node() const;

    /// Number of nodes per cell.
    ///
    /// *Returns*
    ///     uint
    ///         Number of nodes per cell.
    uint nodes_per_cell() const;
  };
}

#endif
