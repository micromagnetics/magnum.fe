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
