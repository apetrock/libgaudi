//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __M2HARMONIC_INTEGRATOR__
#define __M2HARMONIC_INTEGRATOR__

#include "gaudi/asawa/coordinate_interface.hpp"
#include "gaudi/asawa/manifold.hpp"
#include "gaudi/geometry_types.hpp"

#include "tree_code.hpp"
#include <vector>

namespace calder {

template <typename SPACE> class mesh_calculator {
public:
  M2_TYPEDEFS;

  using Dist_Integrator =
      Geometry_Integrator<SPACE, coordinate_type, triangle_type, T>;
  using DTree = typename Dist_Integrator::Tree;
  using DNode = typename Dist_Integrator::Node;

  vector<T> windingNumber(surf_ptr mesh,
                          const std::vector<coordinate_type> &evalPoints,
                          T regLength = 0.0001) {

    auto pre = [](const vector<triangle_type> &tris, DNode &node, DTree &tree,
                  coordinate_type &netCharge, coordinate_type &avgPoint,
                  coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      avgNormal = coordinate_type(0, 0, 0);
      netCharge = coordinate_type(0, 0, 0);

      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
        // netCharge += w;
        netWeight += w;
      }
      // std::cout << "w: " << netWeight << std::endl;
      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
      T kappa = 0.25 / M_PI / dist3;

      return kappa;
    };

    auto compute =
        [this, regLength, computeK](
            int i_c, const coordinate_type &q, const coordinate_type &pc,
            const coordinate_type &pe, const coordinate_type &N,
            const vector<triangle_type> &tris, DNode &node, DTree &tree) -> T {
      T out = 0;
      coordinate_type dp = pc - pe;
      T dist = va::norm(dp);
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          auto tri = tris[ii];
          auto c = tri.center();
          out += tri.solidAngle(pe);
        }
      } else {
        T kappa = computeK(dist, regLength);
        out += kappa * va::dot(N, dp);
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }

    vector<T> u(evalPoints.size(), T(0));
    Dist_Integrator integrator;
    integrator.integrate(normals, triangles, evalPoints, u, pre, compute);
    return u;
  }

} // namespace calder

#endif