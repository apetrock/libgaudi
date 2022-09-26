/*
 *  obj_loader.h
 *  Manifold
 *
 *  Created by John Delaney on 3/23/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <sstream>
#include <vector>

#include <GL/gl.h>
#include <GL/glu.h> // Header File For The GLu32 Library

#include "m2Includes.h"

#include "octree.hpp"
#include <cmath>

#ifndef __OBJ_LOADER__
#define __OBJ_LOADER__

namespace m2 {

template <class SPACE> // unfortunately I can't make this any more generic, but
                       // that's OK
class merge_proc {
  M2_TYPEDEFS;

public:
  typedef Octree<T, vertex_type> Vtree;

  void merge_vertices(surf_ptr obj_in, vertex_ptr v1, vertex_ptr v2,
                      vector<list<face_vertex_ptr>> &corners) {
    std::cout << "merging vertices" << std::endl;
    list<face_vertex_ptr> &va1 = corners[v1->position_in_set()];
    list<face_vertex_ptr> &va2 = corners[v2->position_in_set()];
    size_t sz = va2.size();
    fvl_iterator itb = va2.begin();
    fvl_iterator ite = va2.end();
    for (itb; itb != ite; itb++) {
      v1->add_face_vertex(*itb);
      va1.push_back(*itb);
    }
    va2.clear();
    obj_in->remove_vertex(v2->position_in_set());
  }

  bool operator()(const Vtree &o, surf_ptr data) {
    return this->merge(o, data);
  }

  bool merge(const Vtree &o, surf_ptr data) {
    surf<SPACE> *obj = (surf<SPACE> *)data;

    if (o.is_leaf()) {
      vector<vertex_type *> opts = o.points();
      for (int i = 0; i < opts.size(); i++) {
        vertex_type *p1 = opts[i];
        for (int j = i + 1; j < opts.size(); j++) {
          vertex_type *p2 = opts[j];
          T ds = dist(p1->coordinate(), p2->coordinate());
          if (ds < 0.000001 && p1 && p2) {
            merge_vertices(obj, p2, p1, corners);
            j = opts.size();
          }
        }
      }
      return false;
    } else {
      return true;
    }
  }
  vector<list<face_vertex_ptr>> corners;
}; // merge_proc

template <typename SPACE> class obj_loader : public default_interface<SPACE> {
  M2_TYPEDEFS

public:
  typedef typename SPACE::coordinate_type VecT;
  obj_loader() {}
  ~obj_loader() {}
  m2::surf<SPACE> &operator()(string s) { return this->loadfile(s); }
  // 1. load file into memory
  // 2. instantiate manifold object
  // 3. load up vertices in the object
  // 4. load up faces, for every vertex in the face, make a face vertex and
  // insert that into the appropriate vertex
  // 5. for every vertex loop around all the face vertices and see if there are
  // prev/next faces that share an adjacent vertex 	 if so, make an edge
  typedef Octree<T, vertex_type> Vtree;
  Vtree *build_octree(surf_ptr in_) {
    std::cout << "building octree: " << std::endl;
    Vtree *tree_out = new Octree<T, vertex_type>();
    vertex_array &va = in_->get_vertices();
    T *bounds = tree_out->calculate_cubic_bounds(va);
    coordinate_type avg;

    avg[0] = bounds[0];
    avg[1] = bounds[1];
    avg[2] = bounds[2];
    T r = bounds[3];

    bool built = tree_out->build(va, 4, 10, avg, r, 0);

    in_->reset_flags();
    return tree_out;
  }

  void merge_vertices(surf_ptr in_, vector<list<face_vertex_ptr>> &corners) {
    Vtree *oct = this->build_octree(in_);
    merge_proc<SPACE> merge;
    merge.corners = corners;
    traverse<Vtree, merge_proc<SPACE>, surf_type>(*oct, merge, in_);
    // oct->traverse(merge,(void*) in_);
  }

  coordinate_type find_min(surf<T> *in) {
    in->pack();
    vector<vertex_ptr> &verts = in->get_vertices();
    coordinate_type out;
    coordinate_type avg;
    for (int i = 0; i < verts.size(); i++) {
      avg += verts[i]->coordinate();
    }
    avg /= (T)verts.size();
    out = avg;

    for (int i = 0; i < verts.size(); i++) {
      out[0] = verts[i]->coordinate()[0] < out[0] ? verts[i]->coordinate()[0]
                                                  : out[0];
      out[1] = verts[i]->coordinate()[1] < out[1] ? verts[i]->coordinate()[1]
                                                  : out[1];
      out[2] = verts[i]->coordinate()[2] < out[2] ? verts[i]->coordinate()[2]
                                                  : out[2];
    }
    return out;
  }

  coordinate_type find_max(surf<T> *in) {
    in->pack();
    vector<vertex_ptr> &verts = in->get_vertices();
    coordinate_type out;
    coordinate_type avg;
    for (int i = 0; i < verts.size(); i++) {
      avg += verts[i]->coordinate();
    }
    avg /= (T)verts.size();
    out = avg;

    for (int i = 0; i < verts.size(); i++) {
      coordinate_type c = verts[i]->coordinate();
      if (c[0] > out[0])
        out[0] = c[0];
      if (c[1] > out[1])
        out[1] = c[1];
      if (c[2] > out[2])
        out[2] = c[2];
    }
    return out;
  }

  void stitch_edges(surf_ptr obj, vertex_ptr vert,
                    list<face_vertex_ptr> &rFaceVertices) {
    fvl_iterator itb = rFaceVertices.begin();
    fvl_iterator ite = rFaceVertices.end();
    for (itb; itb != ite; itb++) {
      fvl_iterator jtb = rFaceVertices.begin();
      fvl_iterator jte = rFaceVertices.end();
      face_vertex_ptr fvn = *itb;
      if (fvn->edge())
        continue;
      for (jtb; jtb != jte; jtb++) {
        face_vertex_ptr fvp = *jtb;

        if (fvp->prev()->edge())
          continue;
        if (fvn->next()->vertex() == fvp->prev()->vertex()) {
          edge_ptr e = new edge_type(fvn, fvp->prev());
          obj->push_edge(e);
        }
      }
    }
  }

  surf_ref buildObj(vector<coordinate_type> &inputVerts,
                    vector<vector<int>> &inputFaces) {
    surf_ptr obj = new surf_type();
    int i = 0;
    for (auto &v : inputVerts) {
      vertex_ptr vert = obj->insert_vertex();
      this->coordinate(v, vert);
    }

    vector<list<face_vertex_ptr>> cornersOnVertices;

    cornersOnVertices.resize(inputVerts.size());
    for (int i = 0; i < inputFaces.size() - 1; i++) {

      vector<face_vertex_ptr> cornersOnPoly;

      for (int j = 0; j < inputFaces[i].size(); j++) {
        face_vertex_ptr fv = new face_vertex_type();
        int verti = inputFaces[i][j];
        obj->vertex(verti)->set_front(fv);
        cornersOnVertices[verti].push_back(fv);
        cornersOnPoly.push_back(fv);
      }

      for (int j = 0; j < cornersOnPoly.size(); j++) {
        int jn = (j + 1) % cornersOnPoly.size();
        face_vertex_ptr prev = cornersOnPoly[j];
        face_vertex_ptr next = cornersOnPoly[jn];
        prev->set_next(next);
      }

      face_ptr f = new face_type();
      f->set_front(cornersOnPoly[0]);
      f->update_all();
      obj->push_face(f);
    }

    // manage<SPACE> man;
    // this->merge_vertices(obj, cornersOnVertices);

    vector<vertex_ptr> &vertices = obj->get_vertices();

    for (long i = 0; i < vertices.size(); i++) {
      list<face_vertex_ptr> &corners = cornersOnVertices[i];
      if (corners.size() > 0) {
        m2::construct<SPACE> cons;
        this->stitch_edges(obj, vertices[i], corners);
      }
    }

    m2::construct<SPACE> cons;
    std::cout << "filling holes ... ";
    cons.fill_holes(obj, cornersOnVertices);

    // cons.incremental_fill_holes(obj);
#if 1
    for (long i = 0; i < vertices.size(); i++) {
      if (vertices[i]) {
        size_t k = vertices[i]->size();
        if (k == 0)
          obj->remove_vertex(vertices[i]->position_in_set());
      }
    }
#endif
    // obj->clean_up();
    obj->pack();
    obj->update_all();

    return *obj;
  }

  surf_ref loadfile(string s) {

    string line;
    string fname = s;
    ifstream myfile(fname.c_str());

    bool bNormalsIncluded = 0;
    bool bFacesIncluded = 0;
    vector<coordinate_type> inputVerts;
    vector<vector<int>> inputFaces;
    cout << "loading " << s << " from " << fname << " . . . " << endl;
    // Parse File and push_back vectors for FACE, VERTS, MAP, and NORMALS
    if (myfile.is_open()) {
      cout << ".obj file open" << endl;
      while (!myfile.eof()) {

        getline(myfile, line);

        istringstream ss(line);
        vector<string> vals;
        while (ss) {
          string s;
          if (!getline(ss, s, ' '))
            break;
          vals.push_back(s);
        }
        if (vals.size() == 0)
          continue;
        int in = line.find_first_of(" ");
        string a = line.substr(0, in);

        int vIdx, fIdx = 0;

        if (vals[0] == "v") {

          string fvs = line.substr(2);
          float vx, vy, vz;
          sscanf(fvs.c_str(), "%f %f %f", &vx, &vy, &vz);
          coordinate_type v(vx, vy, vz);

          inputVerts.push_back(v);

          vIdx++;
          // cout << mPoints[mPoints.size()-1] << endl;

        } else if (vals[0] == "f") {
          bFacesIncluded = 1;
          int faceSize = vals.size() - 1;
          vector<int> curFace;

          for (int i = 1; i < vals.size(); i++) {

            int index;
            std::istringstream iss(vals[i]);
            if (iss >> index)
              curFace.push_back(index - 1);
          }

          inputFaces.push_back(curFace);
          fIdx++;
        } else if (vals[0] == "vn") {
          bNormalsIncluded = 1;

          string fvs = line.substr(3);
          int fi = fvs.find_first_of(" ");
          string fv1 = fvs.substr(0, fi);
          string fvs2 = fvs.substr(fi + 1);
          fi = fvs2.find_first_of(" ");
          string fv2 = fvs2.substr(0, fi);
          string fv3 = fvs2.substr(fi + 1, 10);

          VecT normal;

          stringstream(fv1) >> normal[0];
          stringstream(fv2) >> normal[1];
          stringstream(fv3) >> normal[2];

          // mNormal.push_back(normal);
        }
      }
    }
    myfile.close();

    surf_ref obj = buildObj(inputVerts, inputFaces);

    std::cout << "done" << std::endl;
    // obj.clean_up();
    obj.pack();
    obj.update_all();
    return obj;
  }
};

template <typename SPACE>
void write_obj(m2::surf<SPACE> &obj, string filename) {
  M2_TYPEDEFS;
  cout << " Writing file " << filename.c_str() << " ...";
  flush(cout);
  ofstream fout(filename.c_str());
  vector<vertex_ptr> vertices = obj.get_vertices();
  vector<face_ptr> faces = obj.get_faces();
  // dump vertex positions
  for (unsigned int x = 0; x < vertices.size(); x++) {
    coordinate_type vertex = ci::get_coordinate<SPACE>(vertices[x]);
    fout << "v " << vertex[0] << " " << vertex[1] << " " << vertex[2] << endl;
  }

  for (unsigned int x = 0; x < faces.size(); x++) {
    face_ptr f = faces[x];

    if (f->size() < 3)
      continue;

    vertex_ptr vert0 = f->fbegin()->vertex();
    vertex_ptr vert1 = f->fbegin()->next()->vertex();
    vertex_ptr vert2 = f->fbegin()->prev()->vertex();
    // if(f->size() != 3) std::cout << "size: "<< f->size() << std::endl;

    fout << "f " << vert0->position_in_set() + 1 << " "
         << vert1->position_in_set() + 1 << " " << vert2->position_in_set() + 1
         << endl;
  }

  fout.close();
  cout << " done." << endl;
}

template <typename SPACE>
void write_pbrt(m2::surf<SPACE> &obj, const char *filename) {
  M2_TYPEDEFS;
  FILE *file = fopen(filename, "w");

  // read in the vertices
  fprintf(file, " Shape \"trianglemesh\"\n");
  fprintf(file, "\"point P\" [\n");
  vector<vertex_ptr> vertices = obj.get_vertices();
  vector<face_ptr> faces = obj.get_faces();

  for (unsigned int x = 0; x < vertices.size(); x++) {
    coordinate_type vertex = vertices[x]->coordinate();
    fprintf(file, "%f %f %f\n", vertex[0], vertex[1], vertex[2]);
  }
  fprintf(file, "]\n");

  fprintf(file, "\"integer indices\" [\n");
  for (unsigned int x = 0; x < faces.size(); x++) {
    face_ptr f = faces[x];
    vertex_ptr vert0 = f->fbegin()->vertex();
    vertex_ptr vert1 = f->fbegin()->next()->vertex();
    vertex_ptr vert2 = f->fbegin()->prev()->vertex();

    fprintf(file, "%i %i %i\n", vert0->position_in_set(),
            vert1->position_in_set(), vert2->position_in_set());
  }
  fprintf(file, "]\n");

  fclose(file);
}

} // namespace m2

#endif
