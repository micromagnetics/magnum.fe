#ifndef _NORMALIZED_VECTOR_H_
#define _NORMALIZED_VECTOR_H_

#include "CGDofForm.h"
#include <dolfin.h>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace magnumfe {
  class NormalizedVector: public CGDofForm
  {
    public:

    NormalizedVector(
        boost::shared_ptr<const dolfin::FunctionSpace> V1,
        boost::shared_ptr<const dolfin::GenericFunction> a,
        double length = 1.0) : CGDofForm(1, 1), _length(length)
    {
      _function_spaces[0] = V1;
      _function_spaces[1] = V1;
      set_coefficient(0, a);
    }

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
    double _length;
  };
}

#endif
