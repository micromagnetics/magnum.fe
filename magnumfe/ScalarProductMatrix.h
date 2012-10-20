#ifndef _SCALAR_PRODUCT_H_
#define _SCALAR_PRODUCT_H_

#include "DofForm.h"
#include <dolfin.h>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace magnumfe {
  class ScalarProductMatrix: public DofForm
  {
    public:

    ScalarProductMatrix(
        boost::shared_ptr<const dolfin::FunctionSpace> V1,
        boost::shared_ptr<const dolfin::FunctionSpace> V2,
        boost::shared_ptr<const dolfin::GenericFunction> a) : DofForm(2, 1)
    {
      std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > function_spaces;
      _function_spaces[0] = V1;
      _function_spaces[1] = V2;
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

    virtual void eval(double* A, const double * const * w) const
    {
      // do something
    }

    virtual void cell_sparsity(std::vector<std::vector<uint> >& entries) const
    {
      // do something
    }

  };
}

#endif
