// Copyright (C) 2011-2014 Claas Abert
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

#ifndef _SCALAR_PRODUCT_MATRIX_H_
#define _SCALAR_PRODUCT_MATRIX_H_

#include "CGDofForm.h"
#include <dolfin.h>
#include <vector>

namespace magnumfe {

  /// The class represents a matrix that is used to compute the
  /// nodewise scalar product with a given vector.

  class ScalarProductMatrix: public CGDofForm
  {
  public:

    /// Create a scalar product matrix.
    ///
    /// *Arguments*
    ///     V1 (std::shared_ptr<const dolfin::FunctionSpace>)
    ///         The function space for the nodewise scalar product.
    ///         (scalar field)
    ///     V2 (std::shared_ptr<const dolfin::FunctionSpace>)
    ///         The vector function space of the functions to be
    ///         multiplied.
    ///     a (std::shared_ptr<const dolfin::GenercFunction>)
    ///         The function for which the nodewise scalarproduct is
    ///         presented in matrix form.
    ScalarProductMatrix(
        std::shared_ptr<const dolfin::FunctionSpace> V1,
        std::shared_ptr<const dolfin::FunctionSpace> V2,
        std::shared_ptr<const dolfin::GenericFunction> a) : CGDofForm(2, 1)
    {
      std::vector<std::shared_ptr<const dolfin::FunctionSpace> > function_spaces;
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

    virtual void nodeEval(double* A, const double * const * w) const
    {
      const uint dim = _function_spaces[1]->element()->value_dimension(0);
      for (size_t i=0; i<dim; ++i) {
        A[i] = w[0][i];
      }
      //std::cout << "Vec: " << A[0] << "," << A[1] << "," << A[2] << std::endl;
    }
  };
}

#endif
