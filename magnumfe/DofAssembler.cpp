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
    const uint num_coefficients = a.num_coefficients();

    // INITIALIZE SPARSITY PATTERN
    boost::shared_ptr<dolfin::TensorLayout> tensor_layout = A.factory().create_layout(form_rank);
    init_tensor_layout(*tensor_layout, a);
    dolfin::GenericSparsityPattern& pattern = *tensor_layout->sparsity_pattern();

    // fill sparsity pattern
    std::vector<const std::vector<size_t>* > dofs(form_rank);
    boost::multi_array<uint, 2> cell_sparsity(boost::extents[1][1]);
    a.cell_sparsity(cell_sparsity);

    // prepare for filling sparsity pattern
    std::vector<std::vector<size_t> > entries(form_rank);
    for (size_t i=0; i<form_rank; ++i) entries[i].resize(1);
    std::vector<const std::vector<size_t>* > ptr_entries(form_rank);

    // iterate over cells
    for (dolfin::CellIterator cell(*a.mesh()); !cell.end(); ++cell) {
      // get cell dofs
      for (size_t i=0; i<form_rank; ++i)
        dofs[i] = &(a.function_space(i)->dofmap()->cell_dofs(cell->index()));

      // write global sparsity pattern for cell
      for (size_t i=0; i<a.non_zero_entries(); ++i) {
        for (size_t j=0; j<form_rank; ++j) {
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
    std::vector<size_t> num_rows(a.non_zero_entries(), 1);
    std::vector<const size_t*> rows(form_rank);

    std::vector<double> values(a.non_zero_entries());
    std::vector<double *> w(num_coefficients);

    std::vector<std::vector<double> > vector_w(num_coefficients);
    for (size_t i=0; i<num_coefficients; ++i)
      vector_w[i].resize(a.function_space(i + form_rank)->dofmap()->max_cell_dimension());

    for (dolfin::CellIterator cell(*a.mesh()); !cell.end(); ++cell) {
      // get cell dofs for each rank
      for (size_t i=0; i<form_rank; ++i)
        dofs[i] = &(a.function_space(i)->dofmap()->cell_dofs(cell->index()));

      // create w
      for (size_t i=0; i<num_coefficients; ++i) {
        dolfin::Function coeff(a.function_space(i + form_rank));
        coeff.interpolate(*a.coefficient(i));
        const size_t* coeff_dofs = &(coeff.function_space()->dofmap()->cell_dofs(cell->index())[0]);
        coeff.vector()->get_local(&vector_w[i][0], coeff.function_space()->dofmap()->cell_dimension(cell->index()), coeff_dofs);
        w[i] = &vector_w[i][0];
      }
      a.eval(&values[0], &w[0]);

      //const double value = a.eval() TODO should take dofs somehow
      for (size_t i=0; i<a.non_zero_entries(); ++i) {
        const double value = values[i];
        for (size_t j=0; j<form_rank; ++j) {
          const size_t* row = &(*dofs[j])[cell_sparsity[i][j]];
          rows[j] = row;
        }
        //(const double* block, const uint* num_rows, const uint * const * rows)
        A.set(&value, &num_rows[0], &rows[0]);
      }
    }
    A.apply("insert");
  }

  void DofAssembler::init_tensor_layout(dolfin::TensorLayout& tensor_layout,
      const DofForm& a)
  {
    dolfin_assert(tensor_layout);
    const uint rank = a.rank();

    // create list of dofmaps
    std::vector<const dolfin::GenericDofMap*> dofmaps;
    for (size_t i=0; i<rank; ++i) {
      dofmaps.push_back(a.function_space(i)->dofmap().get());
    }

    std::vector<size_t> global_dimensions(rank);
    std::vector<std::pair<size_t, size_t> > local_range(rank);
    for (size_t i = 0; i < rank; i++)
    {
      dolfin_assert(dofmaps[i]);
      global_dimensions[i] = dofmaps[i]->global_dimension();
      local_range[i]       = dofmaps[i]->ownership_range();
    }

    tensor_layout.init(global_dimensions, local_range);                                                     
    if (tensor_layout.sparsity_pattern())
    {
      dolfin::GenericSparsityPattern& pattern = *tensor_layout.sparsity_pattern();
      const std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t> > > periodic_master_slave_dofs;
      dolfin::SparsityPatternBuilder::build(pattern, *(a.mesh()), dofmaps, periodic_master_slave_dofs, false, false, false, false);
    }
  }

}
