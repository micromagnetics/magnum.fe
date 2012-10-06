#include "DofAssembler.h"
#include <dolfin.h>
#include <memory>

namespace magnumfe {
  /*
  void DofAssembler::assemble(dolfin::GenericTensor& A,
      const DofForm& a,
      bool reset_sparsity=true)
  {
    // reset and build sparsity pattern
    // call eval on DofForm to set matrix entries
  }
  */

  void DofAssembler::getMapping(boost::unordered_map<uint, uint>& dofmap,
      const dolfin::FunctionSpace& V1,
      const dolfin::FunctionSpace& V2)
  {
    // collapse function spaces and store mapping
    boost::unordered_map<uint, uint> cmap1;
    const dolfin::FunctionSpace V1_collapsed = V1.dofmap()->is_view() ? *V1.collapse(cmap1) : V1;
    const std::pair<uint, uint> local_range1 = V1_collapsed.dofmap()->ownership_range();

    boost::unordered_map<uint, uint> cmap2;
    const dolfin::FunctionSpace V2_collapsed = V2.dofmap()->is_view() ? *V2.collapse(cmap2) : V2;
    const std::pair<uint, uint> local_range2 = V2_collapsed.dofmap()->ownership_range();

    // setup function in V2
    dolfin::Function f2(V2_collapsed);
    for (uint i=local_range2.first; i<local_range2.second; ++i) {
      const double value = i;
      f2.vector()->set(&value, 1, &i);
    }
    f2.vector()->apply("insert");
    f2.update();

    // interpolate to V1
    dolfin::Vector v1(f2.vector()->size());
    V1_collapsed.interpolate(v1, f2);

    // create a vector with all indices for gathering
    std::vector<uint> gather_all(v1.size());
    for (uint i=0; i<v1.size(); ++i) gather_all[i] = i;

    // create map
    std::vector<uint> map(v1.size());

    // gather value in v1
    dolfin::Vector v1_local(v1.size());
    v1.gather(v1_local, gather_all);

    for (uint i=0; i<v1.size(); ++i) {
      // TODO throw error if value is not a natural number
      map[i] = floor(v1_local[i] + 0.5);
    }

    // modify map for collapsed spaces
    dofmap.clear();
    if (V1.dofmap()->is_view()) {
      // distribute cmap via a dolfin::Vector
      boost::shared_ptr<dolfin::GenericVector> cmap = dolfin::Function(V1_collapsed).vector();
      for (uint i=local_range1.first; i<local_range1.second; ++i) {
        const double value = cmap1[i];
        cmap->set(&value, 1, &i);
      }
      cmap->apply("insert");

      dolfin::Vector cmap_local(cmap->size());
      cmap->gather(cmap_local, gather_all);

      // Apply mapping
      for (uint i=0; i<cmap_local.size(); ++i) {
        dofmap[cmap_local[i]] = map[i];
      }
    }
    else {
      for (uint i=0; i<map.size(); ++i) {
        dofmap[i] = map[i];
      }
    }

    if (V2.dofmap()->is_view()) {
      // distribute cmap via a dolfin::Vector
      boost::shared_ptr<dolfin::GenericVector> cmap = dolfin::Function(V2_collapsed).vector();
      for (uint i=local_range2.first; i<local_range2.second; ++i) {
        const double value = cmap2[i];
        cmap->set(&value, 1, &i);
      }
      cmap->apply("insert");

      dolfin::Vector cmap_local(cmap->size());
      cmap->gather(cmap_local, gather_all);

      // Apply mapping
      for (boost::unordered_map<uint, uint>::iterator it = dofmap.begin(); it != dofmap.end(); ++it) {
        //it->second = map2[it->second];
        it->second = cmap_local[it->second];
      }
    }
  }
}
