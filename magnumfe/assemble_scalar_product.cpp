#include <dolfin.h>
#include "assemble_scalar_product.h"
#include <iterator>
#include <algorithm>

namespace magnumfe {

  void init_tensor_layout(dolfin::TensorLayout& tensor_layout,
    const dolfin::FunctionSpace& V1,
    const dolfin::FunctionSpace& V2)
  {
    dolfin_assert(tensor_layout);
    const uint rank = 2;

    // create list of dofmaps
    std::vector<const dolfin::GenericDofMap*> dofmaps;
    dofmaps.push_back(V1.dofmap().get());
    dofmaps.push_back(V2.dofmap().get());


    std::vector<uint> global_dimensions(rank);
    std::vector<std::pair<uint, uint> > local_range(rank);
    for (uint i = 0; i < rank; i++)
    {
      dolfin_assert(dofmaps[i]);
      global_dimensions[i] = dofmaps[i]->global_dimension();
      local_range[i]       = dofmaps[i]->ownership_range();
    }

    tensor_layout.init(global_dimensions, local_range);                                                     
    if (tensor_layout.sparsity_pattern())
    {
      dolfin::GenericSparsityPattern& pattern = *tensor_layout.sparsity_pattern();
      const std::vector<std::pair<std::pair<uint, uint>, std::pair<uint, uint> > > periodic_master_slave_dofs;
      dolfin::SparsityPatternBuilder::build(pattern, *(V1.mesh()), dofmaps, periodic_master_slave_dofs, false, false, false, false);
    }
  }

  void init_mapping(const dolfin::FunctionSpace& VV,
      const dolfin::FunctionSpace& VS)
  {
    //dolfin::FunctionSpace Vx();
    //boost::unordered_map<uint, uint> Mx, My, Mz;

    //dolfin::FunctionSpace Vx = VV.sub(0).collapse(Mx);
    //Vy = VV.sub(1).collapse(My);
    //Vz = VV.sub(2).collapse(Mz);
  }

  void assemble_scalar_product(const dolfin::GenericVector& x, const dolfin::FunctionSpace& VV, const dolfin::FunctionSpace& VS, dolfin::GenericMatrix& A) {
    //TODO assert VV.mesh == VS.mesh
    boost::shared_ptr<dolfin::TensorLayout> tensor_layout = A.factory().create_layout(2);
    init_tensor_layout(*tensor_layout, VS, VV);
    dolfin::GenericSparsityPattern& pattern = *tensor_layout->sparsity_pattern();

    // get local range of matrix
    const std::pair<uint, uint> row_range = pattern.local_range(0);

    // create sparsity pattern
    for (uint i = row_range.first; i < row_range.second; ++i) {
      for (uint j = 0; j < 3; ++j) {
        std::vector<uint> k(1, i);
        std::vector<uint> l(1, i+j*VS.dim());

        std::vector<const std::vector<uint>*> entries;
        entries.push_back(&k);
        entries.push_back(&l);

        pattern.insert(entries);
      }
    }
    pattern.apply();

    // initialize matrix
    A.init(*tensor_layout);
    A.zero();

    // gather values of vector
    std::vector<uint>   indices;
    std::vector<double> local_x;
    for (uint i = 0; i<VS.dim()*3; ++i)
      indices.push_back(i);

   // for (uint i = row_range.first; i < row_range.second; ++i)
   //   for (uint j = 0; j < 3; ++j)
   //     indices.push_back(i+j*VS.dim());
    x.gather(local_x, indices);
    std::copy(indices.begin(), indices.end(), std::ostream_iterator<uint>(std::cout, ", "));
    std::cout << std::endl;
    std::copy(local_x.begin(), local_x.end(), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
    std::vector<double> bla;
    x.get_local(bla);
    std::copy(bla.begin(), bla.end(), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;

    // set matrix values
    uint local_idx = 0;
    for (uint i = row_range.first; i < row_range.second; ++i) {
      for (uint j = 0; j < 3; ++j) {
        uint k   = i+j*VS.dim();
        double v = local_x[local_idx];
        local_idx++;

        A.add(&v, 1, &i, 1, &k);
      }
    }
    A.apply("add");
  }

  void assemble_transp_scalar_product(const dolfin::GenericVector& x, const dolfin::FunctionSpace& VV, const dolfin::FunctionSpace& VS, dolfin::GenericMatrix& A) {
    //TODO assert VV.mesh == VS.mesh

    boost::shared_ptr<dolfin::TensorLayout> tensor_layout = A.factory().create_layout(2);
    init_tensor_layout(*tensor_layout, VV, VS);
    dolfin::GenericSparsityPattern& pattern = *tensor_layout->sparsity_pattern();

    for (uint i=0; i<VV.dim(); ++i) {
      std::vector<uint> j(1, i);
      std::vector<uint> k(1, i%VS.dim());

      std::vector<const std::vector<uint>*> entries;
      entries.push_back(&j);
      entries.push_back(&k);

      pattern.insert(entries);
    }
    pattern.apply();

    // initialize matrix
    A.init(*tensor_layout);
    A.zero();

    // set matrix values
    for (uint i=0; i<VV.dim(); ++i) {
      uint j   = i;
      uint k   = i%VS.dim();
      double v = x[i];

      A.add(&v, 1, &j, 1, &k);
    }
    A.apply("add");
  }
}
