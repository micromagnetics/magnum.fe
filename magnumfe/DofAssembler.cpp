#include "DofAssembler.h"
#include "DofForm.h"
#include <dolfin.h>
#include <memory>
#include <boost/multi_array.hpp>

namespace magnumfe {
  void DofAssembler::assemble(dolfin::GenericTensor& A,
      const DofForm& a,
      bool reset_sparsity)
  {
    const uint form_rank = a.rank();

    // INITIALIZE SPARSITY PATTERN
    boost::shared_ptr<dolfin::TensorLayout> tensor_layout = A.factory().create_layout(form_rank);
    init_tensor_layout(*tensor_layout, a);
    dolfin::GenericSparsityPattern& pattern = *tensor_layout->sparsity_pattern();

    // fill sparsity pattern
    std::vector<const std::vector<uint>* > dofs(form_rank);
    boost::multi_array<uint, 2> cell_sparsity(boost::extents[1][1]);
    a.cell_sparsity(cell_sparsity);

    // prepare for filling sparsity pattern
    std::vector<std::vector<uint> > entries(form_rank);
    for (uint i=0; i<form_rank; ++i) entries[i].resize(1);
    std::vector<const std::vector<uint>* > ptr_entries(form_rank);

    // iterate over cells
    for (dolfin::CellIterator cell(*a.mesh()); !cell.end(); ++cell) {
      // get cell dofs
      for (uint i=0; i<form_rank; ++i)
        dofs[i] = &(a.function_space(i)->dofmap()->cell_dofs(cell->index()));

      // write global sparsity pattern for cell
      for (uint i=0; i<a.non_zero_entries(); ++i) {
        for (uint j=0; j<form_rank; ++j) {
          const uint value = (*dofs[j])[cell_sparsity[i][j]];
          entries[j][0]  = value;
          ptr_entries[j] = &entries[j];
        }
        pattern.insert(ptr_entries);
      }
    }
    pattern.apply();
    A.init(*tensor_layout);
    A.zero();

    // WRITE VALUES
    for (dolfin::CellIterator cell(*a.mesh()); !cell.end(); ++cell) {
      for (uint i=0; i<form_rank; ++i)
        dofs[i] = &(a.function_space(i)->dofmap()->cell_dofs(cell->index()));

      //for (uint i=0; i<a.non_zero_entries(); ++i) {
      //  for (uint j=0; j<form_rank; ++j) {
      //    const uint value = (*dofs[j])[cell_sparsity[i][j]];
      //  }
      //}
    }
  }

  void DofAssembler::init_tensor_layout(dolfin::TensorLayout& tensor_layout,
      const DofForm& a)
  {
    dolfin_assert(tensor_layout);
    const uint rank = a.rank();

    // create list of dofmaps
    std::vector<const dolfin::GenericDofMap*> dofmaps;
    for (uint i=0; i<rank; ++i) {
      dofmaps.push_back(a.function_space(i)->dofmap().get());
    }

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
      dolfin::SparsityPatternBuilder::build(pattern, *(a.mesh()), dofmaps, periodic_master_slave_dofs, false, false, false, false);
    }
  }

}
