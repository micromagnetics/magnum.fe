// Copyright (C) 2011-2012 Claas Abert
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

#ifndef _DOF_FORM_H_
#define _DOF_FORM_H_

#include <dolfin.h>
#include <map>
#include <vector>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

namespace magnumfe {

  /// This class represents a user-defined DofForm. A DofForm is a
  /// userdefined form that may have an arbitrary cell sparsity
  /// pattern as well as an arbitrary eval method.
  ///
  /// This class is meant to be used for forms not arising from
  /// integration e.g. forms declared on the nodes of CG functions.
  ///
  /// A DofForm is defined by overloading the eval() method as well
  /// as the 
  /// An CGDofForm is defined by overloading a number of methods
  /// including the eval() and sparsity_pattern() method.
  /// 
  /// The rank of the form and the number of coefficients have to
  /// be supplied to the constructor.

  class DofForm
  {
  public:

    /// Create DofForm
    ///
    /// *Arguments*
    ///     rank (uint)
    ///         Rank of the form.
    ///     num_coefficients (uint)
    ///         Number of coefficients of the form.
    DofForm(uint rank, uint num_coefficients);

    /// Destructor
    virtual ~DofForm();

    /// Return function space for given argument
    ///
    /// *Arguments*
    ///     i (uint)
    ///         Number of function space.
    /// *Returns*
    ///     boost::shared_ptr<const dolfin::FunctionSpace>
    ///         Function space.
    boost::shared_ptr<const dolfin::FunctionSpace> function_space(uint i) const;

    /// Return all function spaces
    ///
    /// *Returns*
    ///     std::vector<boost::shared_ptr<const dolfin::FunctionSpace> >
    ///         Function spaces.
    std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > function_spaces() const; 

    /// Set coefficient with given number.
    ///
    /// *Arguments*
    ///     i (uint)
    ///         Number of coefficient.
    ///     coefficient (boost::shared_ptr<const dolfin::GenericFunction> coefficient)
    ///         The coefficient.
    void set_coefficient(uint i,
                         boost::shared_ptr<const dolfin::GenericFunction> coefficient);

    /// Set coefficient with given name.
    ///
    /// *Arguments*
    ///     name (std::string)
    ///         Name of coefficient.
    ///     coefficient (boost::shared_ptr<const dolfin::GenericFunction> coefficient)
    ///         The coefficient.
    void set_coefficient(std::string name,
                         boost::shared_ptr<const dolfin::GenericFunction> coefficient);

    /// Set all coefficients at once.
    ///
    /// *Arguments*
    ///     coefficients (std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction> >)
    ///         The coefficients.
    void set_coefficients(std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction> > coefficients);

    /// Get coefficient with given number.
    ///
    /// *Arguments*
    ///     i (uint)
    ///         Number of coefficient.
    /// *Returns*
    ///     boost::shared_ptr<const dolfin::GenericFunction>
    ///         The coefficient.
    boost::shared_ptr<const dolfin::GenericFunction> coefficient(uint i) const;

    /// Get coefficient with given name.
    ///
    /// *Arguments*
    ///     name (std::string)
    ///         Name of coefficient.
    /// *Returns*
    ///     boost::shared_ptr<const dolfin::GenericFunction>
    ///         The coefficient.
    boost::shared_ptr<const dolfin::GenericFunction> coefficient(std::string name) const;

    /// Get all coefficients.
    ///
    /// *Returns*
    ///     std::vector<boost::shared_ptr<const dolfin::GenericFunction> >
    ///         The coefficients.
    std::vector<boost::shared_ptr<const dolfin::GenericFunction> > coefficients() const;

    /// Returns the rank of the form.
    ///
    /// *Returns*
    ///     unsigned int
    ///         The rank.
    unsigned int rank() const;

    /// Returns the number of coefficients.
    ///
    /// *Returns*
    ///     uint
    ///         The number of coefficients.
    uint num_coefficients() const;

    /// Returns the number of a coefficient with given name.
    ///
    /// *Arguments*
    ///     name (std::string)
    ///         Name of coefficient.
    /// *Returns*
    ///     uint
    ///         Number of the coefficients.
    virtual uint coefficient_number(const std::string & name) const;

    /// Returns the name of a coefficient with given number.
    ///
    /// *Arguments*
    ///     i (uint)
    ///         Number of coefficient.
    /// *Returns*
    ///     std::string
    ///         Name of the coefficient.
    virtual std::string coefficient_name(uint i) const;

    /// Returns the mesh.
    ///
    /// *Returns*
    ///     boost::shared_ptr<const dolfin::Mesh>
    ///         The mesh.
    virtual boost::shared_ptr<const dolfin::Mesh> mesh() const;

    /// Evaluate form for a given cell.
    ///
    /// *Arguments*
    ///     A (double*)
    ///         Linearized result of the evaluation.
    ///     w (const double* double*)
    ///         Linearized values of all coefficients in the cell.
    // TODO add cell_index?
    virtual void eval(double* A, const double * const * w) const = 0;

    /// Return cell coordinates of non-zero entries of the form.
    ///
    /// *Arguments*
    ///     entries (boost::multi_array<uint, 2>&)
    ///         Coordinates of non-zero entries of the form.
    // TODO add cell_index?
    virtual void cell_sparsity(boost::multi_array<uint, 2>& entries) const = 0;

    /// Number of non-zero tensor entries per cell.
    ///
    /// *Returns*
    ///     uint
    ///         Number of non-zero tensor entries.
    // TODO add cell_index?
    virtual uint non_zero_entries() const = 0;

  protected:
    // Function spaces (one for each argument)
    std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > _function_spaces;

    // Coefficients
    std::vector<boost::shared_ptr<const dolfin::GenericFunction> > _coefficients;

    // The rank of the form
    const uint _rank;

  };
}
#endif
