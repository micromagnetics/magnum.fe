#ifndef _DOF_ASSEMBLER_H_
#define _DOF_ASSEMBLER_H_

#include <dolfin.h>
#include "DofForm.h"

namespace magnumfe {
  class DofAssembler {
    public:
      static void assemble(dolfin::GenericTensor& A,
          const DofForm& a,
          bool reset_sparsity=true);

      static void init_tensor_layout(dolfin::TensorLayout& tensor_layout,
          const DofForm& a);
  };
}
#endif
