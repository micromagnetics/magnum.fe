#ifndef _SUB_MESH_INTERPOLATOR_H_
#define _SUB_MESH_INTERPOLATOR_H_

#include <dolfin.h>
#include <vector>

namespace magnumfe {

  class SubMeshInterpolator {
    private:
      std::vector<int> mapping;
    public:
      SubMeshInterpolator(dolfin::FunctionSpace& superspace, dolfin::FunctionSpace& subspace);
      void cut(const dolfin::GenericVector& src, dolfin::GenericVector& target);
      void expand(const dolfin::GenericVector& src, dolfin::GenericVector& target);
  };

}

#endif
