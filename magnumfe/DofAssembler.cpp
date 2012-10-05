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
    boost::unordered_map<uint, uint> map1;

    //boost::shared_ptr<const dolfin::GenericDofMap> hallo = V1.dofmap();
    //dolfin::DofMap bla(map1, *static_cast<const dolfin::DofMap*>(hallo.get()), *V1.mesh(), false);

    //const dolfin::FunctionSpace V1_collapsed = V1.dofmap()->is_view() ? dolfin::DofMap(map1, static_cast<dolfin::DofMap>(*V1.dofmap()), *V1.mesh(), false) : V1;
    const dolfin::FunctionSpace V1_collapsed = V1.dofmap()->is_view() ? *V1.collapse(map1) : V1;

    std::cout << "MAP 1 size: " << map1.size() << std::endl;
    for (boost::unordered_map<uint, uint>::iterator it = map1.begin(); it != map1.end(); ++it) {
      std::cout << "map1: " << it->first << " -> " << it->second << std::endl;
    }

    boost::unordered_map<uint, uint> map2;
    const dolfin::FunctionSpace V2_collapsed = V2.dofmap()->is_view() ? *V2.collapse(map2) : V2;
    std::cout << "MAP 2 size: " << map1.size() << std::endl;

    // setup function in V2
    dolfin::Function f2(V2_collapsed);
    const std::pair<uint, uint> local_range = V2_collapsed.dofmap()->ownership_range();

    //for (uint i=0; i<f2.vector()->size(); ++i) {
    for (uint i=local_range.first; i<local_range.second; ++i) {
      const double value = i;
      f2.vector()->set(&value, 1, &i);
    }
    f2.vector()->apply("insert");
    f2.update();

    // interpolate to V1
    dolfin::Vector v1(f2.vector()->size());
    V1_collapsed.interpolate(v1, f2);

    // create map
    // TODO throw error if value is not a natural number
    std::vector<uint> map(v1.size());

    // gather value in v1
    for (uint i=0; i<v1.size(); ++i) map[i] = i;
    dolfin::Vector v1_local(v1.size());
    v1.gather(v1_local, map);
    for (uint i=0; i<v1.size(); ++i) {
      map[i] = floor(v1_local[i] + 0.5);
    }

    // modify map for collapsed spaces
    dofmap.clear();
    if (V1.dofmap()->is_view()) {
      for (boost::unordered_map<uint, uint>::iterator it = map1.begin(); it != map1.end(); ++it) {
        dofmap[it->second] = map[it->first];
      }
    }
    else {
      for (uint i=0; i<map.size(); ++i) {
        dofmap[i] = map[i];
      }
    }

    if (V2.dofmap()->is_view()) {
      for (boost::unordered_map<uint, uint>::iterator it = dofmap.begin(); it != dofmap.end(); ++it) {
        it->second = map2[it->second];
      }
    }

    // some debug output
    /*
    std::cout << "Collapse map:" << std::endl;
    for (boost::unordered_map<uint, uint>::iterator it = map1.begin(); it != map1.end(); ++it) {
      std::cout << it->first << " -> " << it->second << std::endl;
    }
    for (boost::unordered_map<uint, uint>::iterator it = map2.begin(); it != map2.end(); ++it) {
      std::cout << it->first << " -> " << it->second << std::endl;
    }
    */
  }
}
