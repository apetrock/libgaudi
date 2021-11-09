#ifndef __MESH_HELPER__
#define __MESH_HELPER__

#include <GaudiGraphics/buffers.hpp>
#include <GaudiGraphics/viewer.hpp>
#include <GaudiMath/typedefs.h>
#include <manifold/m2Includes.h>
#include <nanogui/glutil.h>

#include <iostream>
#include <string>

namespace gg {

template <typename SPACE>
void fillBuffer(m2::control<SPACE> *mesh, gg::BufferObjectPtr obj,
                m2::colorRGB col = m2::colorRGB(1.0, 0.5, 0.5, 1.0)) {
  using namespace nanogui;
  M2_TYPEDEFS;
  std::vector<face_ptr> faces = mesh->get_faces();
  std::vector<vertex_ptr> verts = mesh->get_vertices();

  int numVerts = 0;
  int numIndices = 0;

  for (int i = 0; i < faces.size(); i++) {
    numIndices += faces[i]->size();
  }
  numVerts = verts.size();

  obj->fillBuffer([&](gg::BufferObject &o) -> void {
    o.allocateVerts(faces.size(), numVerts);
    auto &indices = o.indices();
    auto &positions = o.positions();
    auto &colors = o.colors();

    for (int i = 0; i < faces.size(); i++) {
      face_vertex_ptr fvb = faces[i]->fbegin();
      face_vertex_ptr fve = faces[i]->fend();
      bool it = true;
      int j = 0;

      while (it) {
        it = fvb != fve;
        indices.col(i)[j] = fvb->vertex()->position_in_set();
        // std::cout << fvb->vertex()->position_in_set() << " ";
        //                    std::cout << fvb->vertex()->position_in_set() <<
        //                    std::endl;
        j++;
        fvb = fvb->next();
        if (j == 3)
          break;
      }
      // std::cout << std::endl;
      // std::cout << indices.col(i) << std::endl;
    }

    for (int i = 0; i < verts.size(); i++) {
      for (int j = 0; j < 3; j++) {
        // std::cout << verts[i]->coordinate()[j] << std::endl;
        positions.col(i)[j] = verts[i]->coordinate()[j];
      }
      colors.col(i)[0] = col.r;
      colors.col(i)[1] = col.g;
      colors.col(i)[2] = col.b;
    }
    // std::cout << indices << std::endl;
  });
};

template <typename SPACE>
void insertPoints(std::vector<typename SPACE::vec3> &points,
                  std::vector<typename SPACE::vec3> &p_colors,
                  gg::PointBufferPtr &obj) {
  using namespace nanogui;
  M2_TYPEDEFS;

  int numVerts = points.size();

  obj->fillBuffer([&](gg::BufferObject &o) -> void {
    o.allocateVerts(numVerts, numVerts);
    auto &indices = o.indices();
    auto &positions = o.positions();
    auto &colors = o.colors();


    for (int i = 0; i < points.size(); i++) {
      indices.col(i)[0] = i;

      for (int j = 0; j < 3; j++) {
        positions.col(i)[j] = points[i][j];
      }

      colors.col(i)[0] = p_colors[i][0];
      colors.col(i)[1] = p_colors[i][1];
      colors.col(i)[2] = p_colors[i][2];
    }
  });
};

} // namespace gg

#endif