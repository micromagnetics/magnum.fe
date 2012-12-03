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

#ifndef _SCALAR_PRODUCT_MATRIX2_H_
#define _SCALAR_PRODUCT_MATRIX2_H_

#include "DofForm.h"
#include <dolfin.h>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace magnumfe {

  /// The class represents a matrix that is used to compute the
  /// nodewise scalar product with a given vector.

  class ScalarProductMatrix2: public DofForm
  {
  public:

    /// Create a scalar product matrix.
    ///
    /// *Arguments*
    ///     V1 (boost::shared_ptr<const dolfin::FunctionSpace>)
    ///         The function space for the nodewise scalar product.
    ///         (scalar field)
    ///     V2 (boost::shared_ptr<const dolfin::FunctionSpace>)
    ///         The vector function space of the functions to be
    ///         multiplied.
    ///     a (boost::shared_ptr<const dolfin::GenercFunction>)
    ///         The function for which the nodewise scalarproduct is
    ///         presented in matrix form.
    ScalarProductMatrix2(
        boost::shared_ptr<const dolfin::FunctionSpace> V1,
        boost::shared_ptr<const dolfin::FunctionSpace> V2,
        boost::shared_ptr<const dolfin::GenericFunction> a) : DofForm(2, 1)
    {
      std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > function_spaces;
      _function_spaces[0] = V1;
      _function_spaces[1] = V2;
      _function_spaces[2] = V2; // first coefficient
      set_coefficient(0, a);
    }

    // Override functions from DofForm
    virtual uint coefficient_number(const std::string & name) const
    {
      if (name == "a") 
        return 0;

      // TODO throw error
      return 0;
    }

    virtual std::string coefficient_name(uint i) const
    {
      switch (i)
      {
        case 0:
          return "a";
      }

      // TODO throw error
      return "unnamed";
    }

    virtual void eval(double* A, const double * const * w) const
    {
      for (size_t i=0; i<non_zero_entries(); ++i) {
        A[i] = w[0][i];
      }
    }

    virtual void cell_sparsity(boost::multi_array<uint, 2>& entries) const
    {
      const uint dim = _function_spaces[0]->dofmap()->max_cell_dimension();
      if (entries.shape()[0] != non_zero_entries() ||
          entries.shape()[1] != rank())
      {
        boost::multi_array<uint, 2>::extent_gen extents;
        entries.resize(extents[non_zero_entries()][rank()]);
      }

      for (uint i=0; i<dim; ++i) {
        for (uint j=0; j<3; ++j) { // TODO get dimension from function space
          entries[i + j*dim][0] = i;
          entries[i + j*dim][1] = i + j*dim;
        }
      }
    }

    virtual uint non_zero_entries() const
    {
      const uint dim = _function_spaces[1]->dofmap()->max_cell_dimension();
      return dim; // TODO get dimension from function space
    }

  };
}

#endif
