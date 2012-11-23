#include "CGDofForm.h"

namespace magnumfe {
  void CGDofForm::eval(double* A, const double * const * w) const
  {
    // TODO implement
  }

  void CGDofForm::cell_sparsity(boost::multi_array<uint, 2>& entries) const
  {
    // TODO implement
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

  uint CGDofForm::non_zero_entries() const
  {
    // TODO implement
    const uint dim = _function_spaces[1]->dofmap()->max_cell_dimension();
    return dim; // TODO get dimension from function space
  }

}
