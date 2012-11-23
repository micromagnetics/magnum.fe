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
    virtual void pointEval(double *A, const double * const *w) const = 0;

    // TODO really need to redeclare without virtual?
    void eval(double* A, const double * const * w) const;
    void cell_sparsity(boost::multi_array<uint, 2>& entries) const;
    uint non_zero_entries() const;
  };
}

#endif
