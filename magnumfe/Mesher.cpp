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
    edge->meshAttributes.nbPointsTransfinite = n[dir[i]];
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
  assert(model->getNumRegions() == 1); // TODO allow multiple regions
  model->deletePhysicalGroups();

  // load entities of sample (edges, vertices, region)
  for (GModel::eiter eit = model->firstEdge(); eit != model->lastEdge(); eit++)
    sample_edges.push_back(*eit);
  for (GModel::fiter fit = model->firstFace(); fit != model->lastFace(); fit++)
    sample_faces.push_back(*fit);
  sample_region = *(model->firstRegion());
  sample_region->addPhysicalEntity(0);

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

  std::vector<std::vector<GFace *> > faces;
  faces.push_back(sample_faces);
  sample_region = model->addVolume(faces);
  sample_region->meshAttributes.method = 1; // Transfinite
  sample_region->addPhysicalEntity(0);
  for (uint i=0; i<size.size(); ++i) sample_size[i] = size[i];

  sample_type = CUBOID;
}
//-----------------------------------------------------------------------------
void Mesher::create_shell(int d, double margin, const dolfin::Array<int>& n, double shell_progression) {
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
      nn[0] = sample_edges[11]->meshAttributes.nbPointsTransfinite;
      nn[1] = sample_edges[1]->meshAttributes.nbPointsTransfinite;
      nn[2] = sample_edges[0]->meshAttributes.nbPointsTransfinite;
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
    air_region->addPhysicalEntity(1);
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
    // i/2+2 leads to 2, 3, 4 for x, y, z transformation area
    shell_regions[i]->addPhysicalEntity(i/2+2);
  }
}
//-----------------------------------------------------------------------------
int Mesher::num_sample_vertices() {
  int result = sample_region->getNumMeshVertices();

  std::vector<GFace*>::iterator fit;
  for (fit = sample_faces.begin(); fit != sample_faces.end(); fit++) {
    result += (*fit)->getNumMeshVertices();
  }

  std::vector<GEdge*>::iterator eit;
  for (eit = sample_edges.begin(); eit != sample_edges.end(); eit++)
    result += (*eit)->getNumMeshVertices();

  std::vector<GVertex*>::iterator vit;
  for (vit = sample_vertices.begin(); vit != sample_vertices.end(); vit++)
    result += (*vit)->getNumMeshVertices();

  return result;
}
//-----------------------------------------------------------------------------
double Mesher::sample_volume() {
  double volume = 0.0;
  std::vector<MTetrahedron*>::iterator it;
  for (it = sample_region->tetrahedra.begin(); it < sample_region->tetrahedra.end(); ++it)
    volume += std::fabs((*it)->getVolume());

  return volume;
}
//-----------------------------------------------------------------------------
double Mesher::get_sample_size(int i) {
  return sample_size[i];
}
//-----------------------------------------------------------------------------
void Mesher::mesh(dolfin::Mesh &mesh, double scale) {
  model->mesh(3);
  if (scale != 1.0)
    model->scaleMesh(scale);

  dolfin::MeshEditor editor;
  editor.open(mesh, "tetrahedron", 3, 3);

  // vertices
  editor.init_vertices(model->getNumMeshVertices());
  std::vector<GEntity*> entities; 
  model->getEntities(entities);       
  dolfin::Point p;

  unsigned int index = 0;
  for(unsigned int i = 0; i < entities.size(); i++) {
    for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++) {
      MVertex *v = entities[i]->mesh_vertices[j];
      //const unsigned int index = v->getNum() - 1;
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
  for(GModel::riter it = model->firstRegion(); it != model->lastRegion(); ++it) {
    tetCount += (*it)->tetrahedra.size();
  }
  editor.init_cells(tetCount);

  std::vector<size_t> v(4);
  index = 0;
  for(GModel::riter it = model->firstRegion(); it != model->lastRegion(); ++it) {
    for(unsigned int i = 0; i < (*it)->tetrahedra.size(); i++) {
      MTetrahedron *tet = (*it)->tetrahedra[i];

      v[0] = tet->getVertex(0)->getIndex();
      v[1] = tet->getVertex(1)->getIndex();
      v[2] = tet->getVertex(2)->getIndex();
      v[3] = tet->getVertex(3)->getIndex();

      editor.add_cell(index, v);
      ++index;
    }
  }

  editor.close();

  // physical regions
  mesh.domains().init(3);
  std::map<size_t, size_t> &subdomains = mesh.domains().markers(3);

  index = 0;
  for(GModel::riter it = model->firstRegion(); it != model->lastRegion(); ++it) {
    for(unsigned int i = 0; i < (*it)->tetrahedra.size(); i++) {
      assert ((*it)->physicals.size() == 1);
      subdomains[index] = (*it)->physicals[0];
      ++index;
    }
  }
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
