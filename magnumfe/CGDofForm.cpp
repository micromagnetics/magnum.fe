// Copyright (C) 2011-2014 Claas Abert
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
// First added:  2012-11-30
// Last changed: 2012-11-30

#include "CGDofForm.h"
#include <vector>

using namespace magnumfe;

//-----------------------------------------------------------------------------
void CGDofForm::eval(double* A, const double * const * w) const
{
  // setup node specific A and w
  std::vector<double *> w_node(num_coefficients());
  std::vector<std::vector<double> > vector_w_node(num_coefficients());
  std::vector<double> A_node(entries_per_node());

  for (uint i=0; i<num_coefficients(); ++i) {
    vector_w_node[i].resize(_function_spaces[i+rank()]->element()->value_dimension(0));
    w_node[i] = &vector_w_node[i][0];
  }

  // call nodeEval for each node
  for (uint node=0; node<nodes_per_cell(); ++node) {
    // setup w_node
    for (uint i=0; i<num_coefficients(); ++i) {
      for (uint j=0; j<vector_w_node[i].size(); ++j) {
        w_node[i][j] = w[i][node + j*nodes_per_cell()];
      }
    }
    // call nodeEval
    nodeEval(&A_node[0], &w_node[0]);

    // copy A_node into A
    // TODO use some copy method?
    for (uint i=0; i<entries_per_node(); ++i) {
      A[i + entries_per_node() * node] = A_node[i];
    }
  }
}
//-----------------------------------------------------------------------------
void CGDofForm::cell_sparsity(boost::multi_array<uint, 2>& entries) const
{
  // check and fix shape of entries multiarray
  if (entries.shape()[0] != non_zero_entries() ||
      entries.shape()[1] != rank())
  {
    boost::multi_array<uint, 2>::extent_gen extents;
    entries.resize(extents[non_zero_entries()][rank()]);
  }

  // build stride array
  std::vector<uint> strides(rank() + 1);
  uint value = 1;
  for (uint i=0; i<rank(); ++i) {
    strides[i] = value;
    value *= _function_spaces[i]->element()->value_dimension(0);
  }
  strides[rank()] = value;

  // build array of non-zero entries
  for (uint node=0; node<nodes_per_cell(); ++node) {
    for (uint i=0; i<entries_per_node(); ++i) {
      const uint idx = i + entries_per_node() * node;
      for (uint r=0; r<rank(); ++r) {
        entries[idx][r] = (i % strides[r+1]) / strides[r] * nodes_per_cell() + node;
        //std::cout << idx << "," << r << "," << entries[idx][r] << std::endl;
      }
    }
  }
}
//-----------------------------------------------------------------------------
uint CGDofForm::non_zero_entries() const
{
  return entries_per_node() * nodes_per_cell();
}
//-----------------------------------------------------------------------------
uint CGDofForm::entries_per_node() const
{
  uint result = 1;
  for (uint i=0; i<rank(); ++i) {
    result *= _function_spaces[i]->element()->value_dimension(0);
  }
  return result;
}
//-----------------------------------------------------------------------------
uint CGDofForm::nodes_per_cell() const
{
  return _function_spaces[0]->element()->space_dimension() /
         _function_spaces[0]->element()->value_dimension(0);

}
//-----------------------------------------------------------------------------
