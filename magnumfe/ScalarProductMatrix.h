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
      A[0]  = w[0][0];
      A[1]  = w[0][1];
      A[2]  = w[0][2];

      A[3]  = w[0][3];
      A[4]  = w[0][4];
      A[5]  = w[0][5];

      A[6]  = w[0][6];
      A[7]  = w[0][7];
      A[8]  = w[0][8];

      A[9]  = w[0][9];
      A[10] = w[0][10];
      A[11] = w[0][11];
    }

    //virtual void cell_sparsity(uint** entries) const
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
