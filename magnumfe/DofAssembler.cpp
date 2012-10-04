#include "DofAssembler.h"
#include <dolfin.h>

namespace magnumfe {
  /*
  void DofAssembler::assemble(dolfin::GenericTensor& A,
      const DofForm& a,
      bool reset_sparsity=true)
  {
    // reset and build sparsity pattern
    // call eval on DofForm to set matrix entries
  }
  */

  void DofAssembler::getMapping(boost::unordered_map<uint, uint>& dofmap,
      const dolfin::FunctionSpace V1,
      const dolfin::FunctionSpace V2)
  {
    //const std::pair<uint, uint> local_range = V1.dofmap().ownership_range();
    dolfin::Function f1(V1);

    // prepare vector
    for (uint i=0; i<f1.vector()->size(); ++i) {
      const double value = i;
      f1.vector()->set(&value, 1, &i);
    }
    f1.vector()->apply("insert");

    // interpolate
    dolfin::Vector v2;
    V2.interpolate(v2, f1);

    // create map
    std::vector<uint> map(v2.size());
    for (uint i=0; i<v2.size(); ++i) map[i] = floor(v2[i] + 0.5);
  }
}
