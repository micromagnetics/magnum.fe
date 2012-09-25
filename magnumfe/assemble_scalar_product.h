#ifndef __ASSEMBLE_SCALAR_PRODUCT_H
#define __ASSEMBLE_SCALAR_PRODUCT_H

#include <dolfin.h>

namespace magnumfe {

  void init_tensor_layout(dolfin::TensorLayout& tensor_layout,
    const dolfin::FunctionSpace& V1,
    const dolfin::FunctionSpace& V2);

  void init_mapping(const dolfin::FunctionSpace& VV,
      const dolfin::FunctionSpace& VS);

  void assemble_scalar_product(const dolfin::GenericVector& x,
      const dolfin::FunctionSpace& VV,
      const dolfin::FunctionSpace& VS,
      dolfin::GenericMatrix& A);

  void assemble_transp_scalar_product(const dolfin::GenericVector& x,
      const dolfin::FunctionSpace& VV,
      const dolfin::FunctionSpace& VS,
      dolfin::GenericMatrix& A);

}

#endif
