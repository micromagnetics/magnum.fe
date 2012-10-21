#ifndef _DOF_FORM_H_
#define _DOF_FORM_H_

#include <dolfin.h>
#include <map>
#include <vector>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

namespace magnumfe {
  class DofForm
  {
  public:
    // Constructors and destructors
    DofForm(uint rank, uint num_coefficients);
    DofForm(std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > function_spaces,
            std::vector<boost::shared_ptr<const dolfin::GenericFunction> > coefficients);
    virtual ~DofForm();

    // Implemented functions
    boost::shared_ptr<const dolfin::FunctionSpace> function_space(uint i) const;
    std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > function_spaces() const; 
    void set_coefficient(uint i,
                         boost::shared_ptr<const dolfin::GenericFunction> coefficient);
    void set_coefficient(std::string name,
                         boost::shared_ptr<const dolfin::GenericFunction> coefficient);
    void set_coefficients(std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction> > coefficients);
    boost::shared_ptr<const dolfin::GenericFunction> coefficient(uint i) const;
    boost::shared_ptr<const dolfin::GenericFunction> coefficient(std::string name) const;
    std::vector<boost::shared_ptr<const dolfin::GenericFunction> > coefficients() const;
    uint rank() const;
    uint num_coefficients() const;

    // Virtual functions
    virtual uint coefficient_number(const std::string & name) const;
    virtual std::string coefficient_name(uint i) const;
    virtual boost::shared_ptr<const dolfin::Mesh> mesh() const;

    virtual void eval(double* A, const double * const * w) const = 0;
    //virtual void cell_sparsity(uint** entries) const = 0;                      // TODO add cell_index?
    virtual void cell_sparsity(boost::multi_array<uint, 2>& entries) const = 0; // TODO add cell_index?
    virtual uint non_zero_entries() const = 0;                                 // TODO add cell_index?

  protected:
    // Function spaces (one for each argument)
    std::vector<boost::shared_ptr<const dolfin::FunctionSpace> > _function_spaces;

    // Coefficients
    std::vector<boost::shared_ptr<const dolfin::GenericFunction> > _coefficients;

  private:

    const uint _rank;

  };
}
#endif
