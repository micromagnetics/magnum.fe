#include "SubMeshInterpolator.h"
#include <iostream>

namespace magnumfe {

  SubMeshInterpolator::SubMeshInterpolator(dolfin::FunctionSpace& Vsuper, dolfin::FunctionSpace& Vsub): mapping(Vsub.dim())
  {
    // TODO assert that Vsuper and Vsub use the same elements
    // TODO assert that mesh of Vsub is submesh of Vsuper mesh

    //std::cout << "Setup Interpolator ... " << std::flush;

    // do it right
    dolfin::Function fsuper(Vsuper);
    dolfin::Function fsub(Vsub);

    std::pair<size_t, size_t> local_range = Vsuper.dofmap()->ownership_range();
    for (size_t i = local_range.first; i < local_range.second; ++i) {
      const double value = i;
      const dolfin::DolfinIndex idx = i;
      fsuper.vector()->set(&value, 1, &idx);
    }
    fsuper.vector()->apply("insert");

    Vsub.interpolate(*fsub.vector(), fsuper);

    // TODO MPIifize the following
    for (size_t i = 0; i<mapping.size(); ++i) {
      const dolfin::DolfinIndex idx = i;
      mapping[i] = floor((*fsub.vector())[idx] + 0.5);
    }

    //std::cout << "done" << std::endl;

    /*
    for (size_t i=0; i<mapping.size(); ++i) {
      std::cout << mapping[i] << ", ";
    }
    std::cout << std::endl;
    */
  }

  void SubMeshInterpolator::cut(const dolfin::GenericVector& src, dolfin::GenericVector& target) {
    for (size_t i=0; i<mapping.size(); ++i) {
      const double value = src[mapping[i]];
      const dolfin::DolfinIndex idx = i;
      target.set(&value, 1, &idx);
    }
    target.apply("insert");
  }

  void SubMeshInterpolator::expand(const dolfin::GenericVector& src, dolfin::GenericVector& target) {
    target.zero();
    for (uint i=0; i<mapping.size(); ++i) {
      const double value = src[i];
      const dolfin::DolfinIndex j = mapping[i];
      target.set(&value, 1, &j);
    }
    target.apply("insert");
  }
}
