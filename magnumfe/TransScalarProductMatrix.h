#ifndef _TRANS_SCALAR_PRODUCT_MATRIX_H_
#define _TRANS_SCALAR_PRODUCT_MATRIX_H_

#include "CGDofForm.h"
#include <dolfin.h>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace magnumfe {
  class TransScalarProductMatrix: public CGDofForm
  {
    public:

    TransScalarProductMatrix(
        boost::shared_ptr<const dolfin::FunctionSpace> V1,
        boost::shared_ptr<const dolfin::FunctionSpace> V2,
        boost::shared_ptr<const dolfin::GenericFunction> a) : CGDofForm(2, 1)
    {
      std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > function_spaces;
      _function_spaces[0] = V2;
      _function_spaces[1] = V1;
      _function_spaces[2] = V2; // first coefficient
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
      const uint dim = _function_spaces[0]->element()->value_dimension(0);
      for (size_t i=0; i<dim; ++i) {
        A[i] = w[0][i];
      }
      //std::cout << "Vec: " << A[0] << "," << A[1] << "," << A[2] << std::endl;
    }
  };
}

#endif
