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

#ifndef _NORMALIZED_VECTOR_H_
#define _NORMALIZED_VECTOR_H_

#include "CGDofForm.h"
#include <dolfin.h>
#include <vector>

namespace magnumfe {

  /// The class defines a _DofForm_ that represents a nodewise
  /// normalized vector field.

  class NormalizedVector: public CGDofForm
  {
  public:

    /// Create a normalizing form.
    ///
    /// *Arguments*
    ///     V1 (std::shared_ptr<const dolfin::FunctionSpace>)
    ///         The function space for which the form is defined.
    ///     a (std::shared_ptr<const dolfin::GenercFunction>)
    ///         The function to be normalized.
    ///     length (double)
    ///         The length, the vector fuction is normalized to.
    NormalizedVector(
        std::shared_ptr<const dolfin::FunctionSpace> V1,
        std::shared_ptr<const dolfin::GenericFunction> a,
        double length = 1.0) : CGDofForm(1, 1), _length(length)
    {
      _function_spaces[0] = V1;
      _function_spaces[1] = V1;
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

      double norm = 0.0;
      for (size_t i=0; i<dim; ++i) {
        norm += w[0][i] * w[0][i];
      }
      norm = sqrt(norm);

      for (size_t i=0; i<dim; ++i) {
        A[i] = w[0][i] / norm * _length;
      }
    }

  protected:

    // length for normalization
    double _length;
  };
}

#endif
