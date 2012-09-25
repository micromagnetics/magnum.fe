#include "SubMeshInterpolator.h"
#include <iostream>

namespace magnumfe {

  SubMeshInterpolator::SubMeshInterpolator(dolfin::FunctionSpace& Vsuper, dolfin::FunctionSpace& Vsub): mapping(Vsub.dim())
  {
    // TODO assert that Vsuper and Vsub use the same elements
    // TODO assert that mesh of Vsub is submesh of Vsuper mesh

    const double one = 1.0;
    dolfin::Function fsub(Vsub);

    std::vector<double> cell_coefficients(Vsuper.dofmap()->max_cell_dimension());
    dolfin::UFCCell ufc_cell(*Vsuper.mesh(), false);

    // iterate over cells of submesh
    for (dolfin::CellIterator subcell(*Vsub.mesh()); !subcell.end(); ++subcell) {

      const int subindex            = subcell->index();
      const int superindex          = Vsuper.mesh()->intersected_cell(subcell->midpoint());
      const dolfin::Cell& supercell = dolfin::Cell(*Vsuper.mesh(), superindex);

      // iterate over dofs of submesh cell
      for (size_t i=0; i<Vsub.dofmap()->max_cell_dimension(); ++i) {
        const uint subdof   = Vsub.dofmap()->cell_dofs(subindex)[i];

        fsub.vector()->zero();
        fsub.vector()->set(&one, 1, &subdof);
        fsub.vector()->apply("insert");

        ufc_cell.update(supercell);
        fsub.restrict(&cell_coefficients[0], *Vsuper.element(), supercell, ufc_cell);

        // find matching dof in supermesh cell
        for (size_t j=0; i<Vsuper.dofmap()->cell_dimension(subindex); ++j) {
          if (cell_coefficients[j] > 0.5) {
            const uint superdof = Vsuper.dofmap()->cell_dofs(superindex)[j];

            mapping[subdof] = superdof;
            break;
          }
        }
      }
    }

    /*
    for (size_t i=0; i<mapping.size(); ++i) {
      std::cout << mapping[i] << ", ";
    }
    std::cout << std::endl;
    */
  }

  void SubMeshInterpolator::cut(const dolfin::GenericVector& src, dolfin::GenericVector& target) {
    for (uint i=0; i<mapping.size(); ++i) {
      const double value = src[mapping[i]];
      target.set(&value, 1, &i);
    }
    target.apply("insert");
  }

  void SubMeshInterpolator::expand(const dolfin::GenericVector& src, dolfin::GenericVector& target) {
    target.zero();
    for (uint i=0; i<mapping.size(); ++i) {
      const double value = src[i];
      const uint j = mapping[i];
      target.set(&value, 1, &j);
    }
    target.apply("insert");
  }
}
