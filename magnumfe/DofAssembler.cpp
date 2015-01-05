// Copyright (C) 2011-2015 Claas Abert
//
// This file is part of magnum.fe.
//
// magnum.fe is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// magnum.fe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with magnum.fe. If not, see <http://www.gnu.org/licenses/>.
//
// Last modified by Claas Abert, 2015-01-05

#include "DofAssembler.h"
#include "DofForm.h"
#include <dolfin.h>
#include <memory>
#include <boost/multi_array.hpp>

using namespace magnumfe;

void DofAssembler::assemble(dolfin::GenericTensor& A,
    const DofForm& a,
    bool reset_sparsity)
{
  const uint form_rank = a.rank();
  const uint num_coefficients = a.num_coefficients();

  // INITIALIZE SPARSITY PATTERN
  std::shared_ptr<dolfin::TensorLayout> tensor_layout = A.factory().create_layout(form_rank);
  init_tensor_layout(*tensor_layout, a);

  // fill sparsity pattern
  std::vector<const std::vector<dolfin::la_index>* > dofs(form_rank);
  boost::multi_array<uint, 2> cell_sparsity(boost::extents[1][1]);
  a.cell_sparsity(cell_sparsity);


  //std::cout << "Build sparsity ... " << std::flush;
  if (tensor_layout->sparsity_pattern()) {
    dolfin::GenericSparsityPattern& pattern = *tensor_layout->sparsity_pattern();

    // prepare for filling sparsity pattern
    std::vector<std::vector<dolfin::la_index> > entries(form_rank);
    for (size_t i=0; i<form_rank; ++i) entries[i].resize(1);
    std::vector<const std::vector<dolfin::la_index>* > ptr_entries(form_rank);

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
  }
  A.init(*tensor_layout);
  //std::cout << "done" << std::endl;
  A.zero();

  // WRITE VALUES
  std::vector<dolfin::la_index> num_rows(a.non_zero_entries(), 1);
  std::vector<const dolfin::la_index*> rows(form_rank);

  std::vector<double> values(a.non_zero_entries());
  std::vector<double *> w(num_coefficients);

  std::vector<std::vector<double> > vector_w(num_coefficients);
  std::vector<std::shared_ptr<dolfin::Function> > coefficients(num_coefficients);
  for (size_t i=0; i<num_coefficients; ++i) {
    // init value containers
    vector_w[i].resize(a.function_space(i + form_rank)->dofmap()->max_cell_dimension());

    // interpolate coefficients
    std::shared_ptr<dolfin::Function> coeff(new dolfin::Function(a.function_space(i + form_rank)));
    coeff->interpolate(*a.coefficient(i));
    coefficients[i] = coeff;
  }

  //std::cout << "Calculate values ... " << std::flush;
  for (dolfin::CellIterator cell(*a.mesh()); !cell.end(); ++cell) {
    // get cell dofs for each rank
    for (size_t i=0; i<form_rank; ++i)
      dofs[i] = &(a.function_space(i)->dofmap()->cell_dofs(cell->index()));

    // create w
    for (size_t i=0; i<num_coefficients; ++i) {
      const dolfin::la_index* coeff_dofs = &(coefficients[i]->function_space()->dofmap()->cell_dofs(cell->index())[0]);
      coefficients[i]->vector()->get_local(&vector_w[i][0], coefficients[i]->function_space()->dofmap()->cell_dimension(cell->index()), coeff_dofs);
      w[i] = &vector_w[i][0];
    }
    a.eval(&values[0], &w[0]);

    //const double value = a.eval() TODO should take dofs somehow
    for (size_t i=0; i<a.non_zero_entries(); ++i) {
      const double value = values[i];
      for (size_t j=0; j<form_rank; ++j) {
        const dolfin::la_index* row = &(*dofs[j])[cell_sparsity[i][j]];
        rows[j] = row;
      }
      //(const double* block, const uint* num_rows, const uint * const * rows)
      A.set(&value, &num_rows[0], &rows[0]);
    }
  }
  A.apply("insert");
  //std::cout << " done" << std::endl;
}
//-----------------------------------------------------------------------------
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
  std::vector<std::size_t> block_sizes;
  for (size_t i = 0; i < rank; i++)
  {
    dolfin_assert(dofmaps[i]);
    global_dimensions[i] = dofmaps[i]->global_dimension();
    local_range[i]       = dofmaps[i]->ownership_range();
    block_sizes.push_back(dofmaps[i]->block_size);
  }

  // Set block size for sparsity graphs
  std::size_t block_size = 1;
  if (a.rank() == 2)
  {
    const std::vector<std::size_t> _bs(a.rank(), dofmaps[0]->block_size);
    block_size = (block_sizes == _bs) ? dofmaps[0]->block_size : 1;
  }

  tensor_layout.init(a.mesh()->mpi_comm(), global_dimensions, block_size, local_range);
  if (tensor_layout.sparsity_pattern())
  {
    dolfin::GenericSparsityPattern& pattern = *tensor_layout.sparsity_pattern();
    //const std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t> > > periodic_master_slave_dofs;
    dolfin::SparsityPatternBuilder::build(pattern, *(a.mesh()), dofmaps, false, false, false, false);
  }
}
//-----------------------------------------------------------------------------
