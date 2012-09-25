#ifndef _MESHER_H_
#define _MESHER_H_

#include <vector>
#include <dolfin.h>
#include <gmsh/GModel.h>

namespace magnumfe {

  class Mesher {
    private:
      enum SampleType {
        NONE,
        MESHFILE,
        CUBOID
      };

      static const int vertex_data[12][2];
      static const int face_data[6][4];
      static const int face_point_data[6][4];

      GModel *model;
      GRegion *sample_region;
      std::vector<GFace*>   sample_faces;
      std::vector<GEdge*>   sample_edges;
      std::vector<GVertex*> sample_vertices;
      dolfin::Array<double> sample_size;
      SampleType            sample_type;

      void create_cuboid_geo(const dolfin::Array<double>& size, const dolfin::Array<int>& n, std::vector<GVertex*> &vertices, std::vector<GEdge*> &edges, std::vector<GFace*> &faces);

    public:
      double cuboid_mesh_size[3];

      Mesher();
      ~Mesher();

      void read_file(const std::string name);
      void write_file(const std::string name);

      void create_cuboid(const dolfin::Array<double>& size, const dolfin::Array<int>& n);
      void create_shell(int d, double margin = 0.0, const dolfin::Array<int>& n = dolfin::Array<int>(0), double shell_progression = 1.0);
      void mesh(dolfin::Mesh &mesh, double scale = 1.0);

      double get_sample_size(int i);
      double sample_volume();
      int num_sample_vertices();
  };

}
#endif
