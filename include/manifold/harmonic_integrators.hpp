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

#include "tree_code.hpp"

namespace m2 {
template <typename SPACE> class mesh_calculator {
public:
  M2_TYPEDEFS;

  using Dist_Integrator =
      m2::Geometry_Integrator<SPACE, coordinate_type, triangle_type, T>;
  using DTree = typename Dist_Integrator::Tree;
  using DNode = typename Dist_Integrator::Node;

  vector<T> calcDistanceFromMesh(surf_ptr mesh,
                                 std::vector<coordinate_type> evalPoints) {

    T regLength = 0.0001;
    auto pre = [](const vector<triangle_type> &tris, DNode &node, DTree &tree,
                  coordinate_type &netCharge,
                  coordinate_type &avgPoint) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      netCharge = coordinate_type(0, 0, 0);

      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        netCharge += w * tri.normal();
        netWeight += w;
      }
      // std::cout << "w: " << netWeight << std::endl;
      avgPoint /= netWeight;
      // netCharge /= netWeight;
      coordinate_type l0 = avgPoint;
      coordinate_type l1 = l0 + netCharge;
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
            const coordinate_type &q, const coordinate_type &pc,
            const coordinate_type &pe,
            const vector<triangle_type> &tris, DNode &node, DTree &tree) -> T {
      T out = 0;
      coordinate_type dp = pc - pe;
      T dist = m2::va::norm(dp);
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          auto tri = tris[ii];
          auto c = tri.center();
          out += tri.solidAngle(pe);
        }
      } else {
        T kappa = computeK(dist, regLength);
        out += kappa * m2::va::dot(q, dp);
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      std::vector<triangle_type> tris = faces[i]->get_tris();
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



  using Avg_Integrator =
      m2::Geometry_Integrator<SPACE, coordinate_type, triangle_type,
                              coordinate_type>;
  using ATree = typename Avg_Integrator::Tree;
  using ANode = typename Avg_Integrator::Node;

  vector<coordinate_type> harmonicAvg(surf_ptr mesh,
                                      std::vector<coordinate_type> faceVectors,
                                      std::vector<coordinate_type> evalPoints, T regLength = 0.5) {

    auto pre = [faceVectors](const vector<triangle_type> &tris, ANode &node,
                             ATree &tree, coordinate_type &netCharge,
                             coordinate_type &avgPoint) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      netCharge = coordinate_type(0, 0, 0);

      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        netCharge += w * faceVectors[ii];
        netWeight += w;
      }

      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
#if 0
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
#elif 0
      T d2 = dist * dist;
      T l2 = C * C;
      T kappa = 1.0 / pow(d2 + l2, 1.5);
#elif 1
      T d2 = dist * dist;
      T l2 = C * C;
      T kappa = (1.0 - exp(-d2 / l2)) / d2;
#else
      T d2 = dist * dist;
      T l2 = C * C;
      T kappa = 1.0 /  (d2 + l2);
#endif
      return kappa;
    };

    auto compute = [this, faceVectors, regLength, computeK](
                       const coordinate_type &q, const coordinate_type &pc,
                       const coordinate_type &pe,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> coordinate_type {
      coordinate_type out = coordinate_type(0, 0, 0);
      coordinate_type dp = pc - pe;
      T dist = m2::va::norm(dp);

      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          coordinate_type q = faceVectors[ii];
          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          dp = c - pe;
          T dist = m2::va::norm(dp);
          //out += w * computeK(dist, regLength) * va::cross(q,dp);
          out += w * computeK(dist, regLength) * q;
        }
      } else {
        //out += computeK(dist, regLength) * va::cross(q, dp);
        out += computeK(dist, regLength) * q;
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      std::vector<triangle_type> tris = faces[i]->get_tris();
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }

    vector<coordinate_type> u(evalPoints.size(), coordinate_type(0, 0, 0));
    Avg_Integrator integrator;
    integrator.integrate(faceVectors, triangles, evalPoints, u, pre, compute);

    return u;
  }
};
} // namespace m2

#endif