#ifndef _CG_DOF_FORM_H_
#define _CG_DOF_FORM_H_

#include "DofForm.h"
#include <dolfin.h>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace magnumfe {
  class CGDofForm: public DofForm
  {
    public:
    CGDofForm(uint rank, uint num_coefficients): DofForm(rank, num_coefficients) {};
    virtual void nodeEval(double *A, const double * const *w) const = 0;

    // TODO really need to redeclare without virtual?
    void eval(double* A, const double * const * w) const;
    void cell_sparsity(boost::multi_array<uint, 2>& entries) const;
    uint non_zero_entries() const;
    uint entries_per_node() const;
    uint nodes_per_cell() const;
  };
}

#endif
