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
// Last modified by Claas Abert, 2014-06-10

#include "Mesher.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <dolfin.h>
#include <assert.h>

#include <gmsh/Gmsh.h>
#include <gmsh/GModel.h>
#include <gmsh/GEntity.h>
#include <gmsh/MVertex.h>
#include <gmsh/MTetrahedron.h>
#include <gmsh/MTriangle.h>
#include <gmsh/GRegionCompound.h>

using namespace magnumfe;

Mesher::Mesher():sample_size(3) {
  GmshInitialize();
  GmshSetOption("Mesh", "Algorithm3D", 4.); // Frontal Algorithm
  model = new GModel();
  model->setFactory("Gmsh"); 
  sample_type = NONE;
}

Mesher::~Mesher() {
  delete model;
  GmshFinalize();
}
//-----------------------------------------------------------------------------
void Mesher::create_cuboid_geo(const dolfin::Array<double>& size, const dolfin::Array<int>& n, std::vector<GVertex*> &vertices, std::vector<GEdge*> &edges, std::vector<GFace*> &faces) {
  const double dx = size[0];
  const double dy = size[1];
  const double dz = size[2];
  const double lc = 1.0; // TODO what about this value?

  vertices.push_back(model->addVertex(-dx, -dy, -dz, lc));
  vertices.push_back(model->addVertex(-dx, -dy,  dz, lc));
  vertices.push_back(model->addVertex(-dx,  dy, -dz, lc));
  vertices.push_back(model->addVertex(-dx,  dy,  dz, lc));
  vertices.push_back(model->addVertex( dx, -dy, -dz, lc));
  vertices.push_back(model->addVertex( dx, -dy,  dz, lc));
  vertices.push_back(model->addVertex( dx,  dy, -dz, lc));
  vertices.push_back(model->addVertex( dx,  dy,  dz, lc));

  const int dir[12] = {2, 1, 2, 1, 2, 1, 2, 1, 0, 0, 0, 0};
  for (int i=0; i<12; ++i) {
    GEdge *edge = model->addLine(vertices[vertex_data[i][0]], vertices[vertex_data[i][1]]);
    edge->resetMeshAttributes();
    edge->meshAttributes.method = 1; // Transfinite
    edge->meshAttributes.typeTransfinite = 1;
    edge->meshAttributes.coeffTransfinite = 1;
    edge->meshAttributes.nbPointsTransfinite = n[dir[i]] + 1;
    edges.push_back(edge);
  }

  for (int i=0; i<6; ++i) {
    std::vector<GEdge*> a;
    for (int j=0; j<4; ++j) {
      a.push_back(edges[face_data[i][j]]);
    }
    GFace *face = model->addPlanarFace(std::vector<std::vector<GEdge *> >(1, a));
    face->meshAttributes.method = 1; // Transfinite
    face->meshAttributes.transfiniteArrangement = -1; // TODO try 1 here
    faces.push_back(face);
  }
}
//-----------------------------------------------------------------------------
void Mesher::read_file(const std::string name) {
  assert(sample_type == NONE);
  model->load(name);

  // retrieve all physical groups
  std::map<int, std::vector<GEntity*> > groups[4];
  model->getPhysicalGroups(groups);

  // handle simple mesh with just one volume and one surface
  if (groups[2].size() == 1 and groups[3].size() == 1) {
    //assert(model->getNumRegions() == 1); // TODO allow multiple regions
    model->deletePhysicalGroups();

    // load entities of sample (edges, vertices, region)
    for (GModel::fiter fit = model->firstFace(); fit != model->lastFace(); fit++) {
      (*fit)->addPhysicalEntity(1);
      sample_faces.push_back(*fit);
    }

    (*(model->firstRegion()))->addPhysicalEntity(1);
  }
  // handle sophisticated meshes
  else {
    // automatically retrieve outer faces of sample for meshing of shell
    // TODO handle holes
    model->createTopologyFromMesh();

    std::vector<GRegion*> all_regions;
    for (GModel::riter rit = model->firstRegion(); rit != model->lastRegion(); rit++) {
      all_regions.push_back(*rit);
    }
    GRegionCompound compound(model, 0, all_regions);

    std::list<GFace*> faces = compound.faces();
    for (std::list<GFace*>::iterator fit = faces.begin(); fit != faces.end(); fit++) {
      sample_faces.push_back(*fit);
    }

    // LEGACY CODE (requires outer boundary of sample to be defined as face with id "1")
    // require face "1" to be defined
    //assert(groups[2].find(1) != groups[2].end());
  }

  // retrieve box size
  SBoundingBox3d bounds = model->bounds();
  sample_size[0] = std::max(fabs(bounds.min().x()), fabs(bounds.max().x()));
  sample_size[1] = std::max(fabs(bounds.min().y()), fabs(bounds.max().y()));
  sample_size[2] = std::max(fabs(bounds.min().z()), fabs(bounds.max().z()));

  sample_type = MESHFILE;
}
//-----------------------------------------------------------------------------
void Mesher::write_file(const std::string name) {
  model->save(name);
}
//-----------------------------------------------------------------------------
void Mesher::create_cuboid(const dolfin::Array<double>& size, const dolfin::Array<int>& n) {
  assert(sample_type == NONE);
  create_cuboid_geo(size, n, sample_vertices, sample_edges, sample_faces);

  for (int i=0; i < sample_faces.size(); ++i) {
    sample_faces[i]->addPhysicalEntity(1);
  }
  std::vector<std::vector<GFace *> > faces;
  faces.push_back(sample_faces);

  GRegion *sample_region = model->addVolume(faces);
  sample_region->meshAttributes.method = 1; // Transfinite
  sample_region->addPhysicalEntity(1);

  for (uint i=0; i<size.size(); ++i) sample_size[i] = size[i];

  sample_type = CUBOID;
}
//-----------------------------------------------------------------------------
void Mesher::create_shell(int d, double margin, const dolfin::Array<int>& n, double shell_progression) {
  // TODO split up in smaller functions
  assert(sample_type != NONE);

  std::vector<GVertex*> inner_vertices;
  std::vector<GEdge*>   inner_edges;
  std::vector<GFace*>   inner_faces;

  // copy n
  dolfin::Array<int> nn(3);

  switch(sample_type) {
    case CUBOID:
      inner_vertices = sample_vertices;
      inner_edges    = sample_edges;
      inner_faces    = sample_faces;

      // get n values from inner cuboid
      nn[0] = sample_edges[11]->meshAttributes.nbPointsTransfinite - 1;
      nn[1] = sample_edges[1]->meshAttributes.nbPointsTransfinite - 1;
      nn[2] = sample_edges[0]->meshAttributes.nbPointsTransfinite - 1;
      break;
    case MESHFILE:
      {
        // TODO chose clever defaults for margin and n
        SBoundingBox3d bounds = model->bounds();

        // add margin to size
        dolfin::Array<double> size(sample_size.size());

        for (uint i=0; i<size.size(); ++i) {
          size[i] = sample_size[i] + margin;
        }

        create_cuboid_geo(size, n, inner_vertices, inner_edges, inner_faces);
        for (uint i=0; i<nn.size(); ++i) nn[i] = n[i];
        break;
      }
    case NONE:
      assert(false);
  }

  std::vector<GVertex*> outer_vertices;
  std::vector<GEdge*>   outer_edges;
  std::vector<GFace*>   outer_faces;

  double shell_width = std::min(sample_size[0], std::min(sample_size[1], sample_size[2])) + margin;

  dolfin::Array<double> size(sample_size.size());
  for (uint i=0; i<size.size(); ++i) {
    size[i] = sample_size[i] + shell_width;
  }
  create_cuboid_geo(size, nn, outer_vertices, outer_edges, outer_faces);

  // setup connections between cuboids
  std::vector<GEdge*> connection_edges;
  for (int i=0; i<8; ++i) {
    GEdge *connection = model->addLine(inner_vertices[i], outer_vertices[i]);
    connection->meshAttributes.method = 1; // Transfinite
    connection->meshAttributes.typeTransfinite = 1;
    connection->meshAttributes.coeffTransfinite = shell_progression;
    connection->meshAttributes.nbPointsTransfinite = d + 1;
    connection_edges.push_back(connection);
  }

  std::vector<GFace*> connection_faces;
  for (int i=0; i<12; ++i) {
    std::vector<GEdge*> a;
    a.push_back(inner_edges[i]);
    a.push_back(connection_edges[vertex_data[i][0]]);
    a.push_back(outer_edges[i]);
    a.push_back(connection_edges[vertex_data[i][1]]);

    GFace* face = model->addPlanarFace(std::vector<std::vector<GEdge *> >(1, a));
    face->meshAttributes.method = 1; // Transfinite
    connection_faces.push_back(face);
  }

  if (sample_type == MESHFILE) {
    // add air region
    std::vector<std::vector<GFace *> > faces;
    faces.push_back(inner_faces);
    faces.push_back(sample_faces);
    GRegion* air_region = model->addVolume(faces);
    air_region->addPhysicalEntity(1000);
  }

  // shell
  std::vector<GRegion*> shell_regions;
  for (int i=0; i<6; ++i) {
    std::vector<GFace*> a;
    a.push_back(inner_faces[i]);
    a.push_back(outer_faces[i]);
    for (int j=0; j<4; ++j) {
      a.push_back(connection_faces[face_data[i][j]]);
    }
    GRegion* region = model->addVolume(std::vector<std::vector<GFace *> >(1, a));
    region->meshAttributes.method = 1; // Transfinite
    shell_regions.push_back(region);
  }

  for (int i=0; i<6; ++i) {
    // i/2+2 leads to 1001, 1002, 1003 for x, y, z transformation area
    shell_regions[i]->addPhysicalEntity(i/2+1001);
  }
}
//-----------------------------------------------------------------------------
double Mesher::get_sample_size(int i) {
  return sample_size[i];
}
//-----------------------------------------------------------------------------
void Mesher::mesh(dolfin::Mesh &mesh, double scale) {
  model->mesh(3);
  if (scale != 1.0) { model->scaleMesh(scale); }

  dolfin::MeshEditor editor;
  editor.open(mesh, "tetrahedron", 3, 3);

  // add vertices to mesh
  editor.init_vertices(model->getNumMeshVertices());
  std::vector<GEntity*> entities; 
  model->getEntities(entities);       
  dolfin::Point p;

  unsigned int index = 0;
  for(unsigned int i = 0; i < entities.size(); i++) {
    for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++) {
      MVertex *v = entities[i]->mesh_vertices[j];
      v->setIndex(index);

      p[0] = v->x();
      p[1] = v->y();
      p[2] = v->z();

      editor.add_vertex(index, p);
      ++index;
    }
  }

  // elements
  // TODO is there a better way to get the total number of tets?
  int tetCount = 0;
  for(GModel::riter rit = model->firstRegion(); rit != model->lastRegion(); ++rit) {
    tetCount += (*rit)->tetrahedra.size();
  }

  editor.init_cells(tetCount);

  // initialize subdomains
  //mesh.domains().init(2); //TODO can this really be omitted?
  mesh.domains().init(3);
  std::map<size_t, size_t> &facetmarkers  = mesh.domains().markers(2);
  std::map<size_t, size_t> &cellmarkers   = mesh.domains().markers(3);

  // add tets to mesh and create subdomains
  std::vector<size_t> v(4);
  index = 0;
  for(GModel::riter rit = model->firstRegion(); rit != model->lastRegion(); ++rit) {
    for(unsigned int i = 0; i < (*rit)->tetrahedra.size(); ++i) {
      MTetrahedron *tet = (*rit)->tetrahedra[i];

      v[0] = tet->getVertex(0)->getIndex();
      v[1] = tet->getVertex(1)->getIndex();
      v[2] = tet->getVertex(2)->getIndex();
      v[3] = tet->getVertex(3)->getIndex();

      editor.add_cell(index, v);

      // apply custom celldomains if tet does not belong to shell
      int physical = (*rit)->physicals[0]; // TODO what if no physicals are present?
      if (physical < 1000 && celldomains.size() > 0) {

        // get center of current tetrahendron
        dolfin::Array<double> tet_center(3);
        tet_center[0] = tet->circumcenter()[0];
        tet_center[1] = tet->circumcenter()[1];
        tet_center[2] = tet->circumcenter()[2];

        // set domain according to custom celldomains
        for (int j=0; j < celldomains.size(); ++j) {
          if (celldomains[j].first->inside(tet_center, false)) {
            physical = celldomains[j].second;
          }
        }
      }

      // set domain marker
      assert ((*rit)->physicals.size() == 1); // TODO move up?
      cellmarkers[index] = physical;

      ++index;
    }
  }
  editor.close();

  // set facet subdomains
  mesh.init(0,2);
  for(GModel::fiter fit = model->firstFace(); fit != model->lastFace(); ++fit) {

    if ((*fit)->getPhysicalEntities().size() == 0 && facetdomains.size() == 0) continue;

    for(unsigned int i = 0; i < (*fit)->triangles.size(); ++i) {
      MTriangle *tri = (*fit)->triangles[i];

      const int v0 = tri->getVertex(0)->getIndex();
      const int v1 = tri->getVertex(1)->getIndex();
      const int v2 = tri->getVertex(2)->getIndex();

      const size_t ft_size   = mesh.topology()(0,2).size(v0);
      const unsigned int *ft = mesh.topology()(0,2)(v0);

      int physical = -1; 
      if ((*fit)->getPhysicalEntities().size() > 0) physical = (*fit)->physicals[0];

      for (int j=0; j < facetdomains.size(); ++j) {
        bool inside = true;
        for (int k=0; k < 3; ++k) {
          dolfin::Array<double> point(3);
          point[0] = tri->getVertex(k)->x();
          point[1] = tri->getVertex(k)->y();
          point[2] = tri->getVertex(k)->z();

          if (!facetdomains[j].first->inside(point, false)) {
            inside = false;
            break;
          }
        }
        if (inside == true) physical = facetdomains[j].second;
      }

      if (physical < 0) continue;

      for(unsigned int j = 0; j < ft_size; ++j) {
        const unsigned int *fv = mesh.topology()(2,0)(ft[j]);

        bool match = true;
        if (std::find(fv, fv + 3, v1) == fv + 3) { match = false; }
        if (std::find(fv, fv + 3, v2) == fv + 3) { match = false; }

        if (match) {
          facetmarkers[ft[j]] = physical;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void Mesher::create_celldomain(std::shared_ptr<const dolfin::SubDomain> subdomain, size_t id) {
  celldomains.push_back(std::pair<std::shared_ptr<const dolfin::SubDomain>, size_t>(subdomain, id));
}
//-----------------------------------------------------------------------------
void Mesher::create_facetdomain(std::shared_ptr<const dolfin::SubDomain> subdomain, size_t id) {
  facetdomains.push_back(std::pair<std::shared_ptr<const dolfin::SubDomain>, size_t>(subdomain, id));
}
//-----------------------------------------------------------------------------
const int Mesher::vertex_data[12][2] = {
      {0,1},{1,3},{3,2},{2,0},
      {4,5},{5,7},{7,6},{6,4},
      {0,4},{1,5},{2,6},{3,7}
};
//-----------------------------------------------------------------------------
const int Mesher::face_data[6][4] = {
      {0, 1, 2, 3},{4, 5, 6, 7},
      {0, 9, 4, 8},{10, 6,11,2},
      {8, 7,10, 3},{9, 5,11, 1}
};
//-----------------------------------------------------------------------------
