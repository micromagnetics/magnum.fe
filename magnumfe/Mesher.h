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

#ifndef _MESHER_H_
#define _MESHER_H_

#include <vector>
#include <dolfin.h>
#include <gmsh/GModel.h>

namespace magnumfe {

  /// This class uses the GMSH library to generate dolfin meshes.
  /// Especially it it offers methods to surround arbitrary meshes
  /// with a cuboid shell with is need for the stray-field
  /// transformation method.
  /// 
  /// Creating a cuboid shell assumes that the coordinate origin is
  /// in the center of the sample mesh (The resulting mesh will be
  /// symmetric around the coordinate origin).
  ///
  /// *Example*
  ///
  ///     .. code-block:: c++
  ///
  ///         Mesher mesher();
  ///         mesher.read_file("sample.msh");
  ///         dolfin::Array<int> n(3);
  ///         // TODO put some values in n
  ///         mesher.create_shell(2, 0.1, n);
  ///         dolfin::Mesh mesh;
  ///         mesher.mesh(mesh);
  //

  class Mesher {
  public:

    /// The size of the mesh
    double cuboid_mesh_size[3];

    /// Create Mesher
    Mesher();

    /// Destroy Mesher
    ~Mesher();

    /// Reads in a mesh file (All filetypes of GMSH are supported).
    ///
    /// *Arguments*
    ///     name (std::string)
    ///         The filename.
    void read_file(const std::string name);

    /// Write the mesh in a file (All filetypes of GMSH are
    /// supported).
    ///
    /// *Arguments*
    ///     name (std::string)
    ///         The filename.
    void write_file(const std::string name);

    /// Creates a cuboid of given dimension around the coordinate
    /// origin.
    ///
    /// *Arguments*
    ///     size (dolfin::Array<double>)
    ///         Size of the cuboid.
    ///     n (dolfin::Array<int>)
    ///         Number of mesh points.
    void create_cuboid(const dolfin::Array<double>& size, const dolfin::Array<int>& n);

    /// Creates a cuboid shell around the current sample.
    ///
    /// *Arguments*
    ///     d (int)
    ///         Number of layers of the shell mesh.
    ///     margin (double)
    ///         Margin between sample and shell
    ///         (useful for meshes read from file).
    ///     n (dolfin::Array<int>)
    ///         Number of mesh points of shell.
    ///         (only needed if mesh was read from file)
    ///     shell_progression (int)
    ///         Defines coarsening of the shell layers.
    ///         1.0 means no coarsening.
    void create_shell(int d, double margin = 0.0, const dolfin::Array<int>& n = dolfin::Array<int>(0), double shell_progression = 1.0);

    /// Mesh the geometry and create dolfin mesh.
    ///
    /// *Arguments*
    ///     mesh (dolfin::Mesh)
    ///         The mesh.
    ///     scale (double)
    ///         Scaling of the mesh. Useful for very small meshes.
    void mesh(dolfin::Mesh &mesh, double scale = 1.0);

    /// Returns the sample size (inner-shell size) in a given
    /// dimension. Scaling used for mesh generation is not considered.
    ///
    /// *Arguments*
    ///     i (int)
    ///         The dimension (direction).
    /// *Returns*
    ///     double
    ///         The size.
    double get_sample_size(int i);

    /// Returns the volume of the sample. Scling used for mesh
    /// generation is not considered.
    ///
    /// *Returns*
    ///     double
    ///         The volume.
    double sample_volume();

    /// Returns the number vertices in the sample.
    ///
    /// *Returns*
    ///     int
    ///         The number of sample vertices.
    int num_sample_vertices();

  protected:
    // Type of sample (either cuboid of read from file)
    enum SampleType {
      NONE,     // not initialized
      MESHFILE, // read from mesh file
      CUBOID    // created via create_cuboid
    };

    // Lists of indices for the shell creation
    static const int vertex_data[12][2];
    static const int face_data[6][4];
    static const int face_point_data[6][4];

    // GMSH geometry and mesh containers
    GModel *model;
    GRegion *sample_region;
    std::vector<GFace*>   sample_faces;
    std::vector<GEdge*>   sample_edges;
    std::vector<GVertex*> sample_vertices;
    dolfin::Array<double> sample_size;
    SampleType            sample_type;

    // Helper method, which creates a GMSH cuboid
    void create_cuboid_geo(const dolfin::Array<double>& size,
                           const dolfin::Array<int>& n,
                           std::vector<GVertex*> &vertices,
                           std::vector<GEdge*> &edges,
                           std::vector<GFace*> &faces);

  };

}
#endif
