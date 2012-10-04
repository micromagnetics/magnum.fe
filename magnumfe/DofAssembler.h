#ifndef _DOF_ASSEMBLER_H_
#define _DOF_ASSEMBLER_H_

#include <dolfin.h>

namespace magnumfe {
  class DofAssembler {
    public:
      //static void assemble(dolfin::GenericTensor& A,
      //    const DofForm& a,
      //    bool reset_sparsity=true);

      static void getMapping(boost::unordered_map<uint, uint>& dofmap,
          const dolfin::FunctionSpace V1,
          const dolfin::FunctionSpace V2);
  };

  /*
  class DofForm {

    public:
      dolfin::FunctionSpace V1, V2;

      DofForm(dolfin::FunctionSpace V1, dolfin::FunctionSpace V2);

      virtual void eval(Array<double>& values,
          const Array<double>& a,
          const Array<double>& b);

    private:
  }
  */
}
#endif
