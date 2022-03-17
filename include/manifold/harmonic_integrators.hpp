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
                  coordinate_type &avgPoint, coordinate_type &avgNormal) -> void {
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
        //netCharge += w;
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

    auto compute = [this, regLength, computeK](
                       const coordinate_type &q, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, DNode &node,
                       DTree &tree) -> T {
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
        out += kappa * m2::va::dot(N, dp);
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      std::vector<triangle_type> tris = m2::ci::get_tris<SPACE>(faces[i]);
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

  template <typename Q>
  vector<Q> harmonicAvg(surf_ptr mesh, std::vector<Q> faceVectors,
                        std::vector<coordinate_type> evalPoints,
                        T regLength = 0.5) {
    using Avg_Integrator = m2::Geometry_Integrator<SPACE, Q, triangle_type, Q>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;

    auto pre = [faceVectors](const vector<triangle_type> &tris, ANode &node,
                             ATree &tree, Q &netCharge,
                             coordinate_type &avgPoint,
                             coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      netCharge = z::zero<Q>();

      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
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
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 1
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3 + l3);
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 0
      T dist2 = dist * dist;
      T dt = 0.5;
      T kappa = exp(-dist2 / 4.0 / dt);
      return kappa / pow(4.0 * M_PI * dt, 1.5);
#endif
      
    };

    auto compute = [this, faceVectors, regLength, computeK](
                       const Q &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> Q {
      Q out = z::zero<Q>();

    
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          Q qi = faceVectors[ii];
          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type dp = c - pe;
          T dist = m2::va::norm(dp);

          //if(dist <= std::numeric_limits<T>::epsilon())
          //  continue;
          // out += w * computeK(dist, regLength) * va::cross(q,dp);

          //T Ndp = m2::va::dot(tri.normal(), dp);
          
          T k = computeK(dist, regLength);
          out += w * k * qi;
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pc - pe;
        T dist = m2::va::norm(dp);
        //T Ndp = m2::va::dot(N, dp);
        T k = computeK(dist, regLength);
        out += k * wq;
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      if (!mesh->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<triangle_type> tris = m2::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }

    vector<Q> u(evalPoints.size(), z::zero<Q>());
    Avg_Integrator integrator;
    integrator.integrate(faceVectors, triangles, evalPoints, u, pre, compute);

    return u;
  }

  template <typename Q>
  vector<Q> harmonicNormal(surf_ptr mesh, std::vector<Q> faceVectors,
                        std::vector<coordinate_type> evalPoints,
                        T regLength = 0.5) {
    using Avg_Integrator = m2::Geometry_Integrator<SPACE, Q, triangle_type, Q>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;

    auto pre = [faceVectors](const vector<triangle_type> &tris, ANode &node,
                             ATree &tree, Q &netCharge,
                             coordinate_type &avgPoint,
                             coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      netCharge = z::zero<Q>();

      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
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
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 1
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3 + l3);
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 0
      T dist2 = dist * dist;
      T dt = 0.5;
      T kappa = exp(-dist2 / 4.0 / dt);
      return kappa / pow(4.0 * M_PI * dt, 1.5);
#endif

    };

    auto compute = [this, faceVectors, regLength, computeK](
                       const Q &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> Q {
      Q out = z::zero<Q>();
      
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          Q qi = faceVectors[ii];
          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type dp = c - pe;
          T dist = m2::va::norm(dp);

          // if(dist <= std::numeric_limits<T>::epsilon())
          //  continue;
          // out += w * computeK(dist, regLength) * va::cross(q,dp);

          T Ndp = m2::va::dot(tri.normal(), qi);

          T k = Ndp * computeK(dist, regLength);
          out += w * k * qi;
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pc - pe;
        T dist = m2::va::norm(dp);
        T Ndp = m2::va::dot(N, wq);
        T k = Ndp * computeK(dist, regLength);
        out += k * wq;
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      std::vector<triangle_type> tris = m2::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }

    vector<Q> u(evalPoints.size(), z::zero<Q>());
    Avg_Integrator integrator;
    integrator.integrate(faceVectors, triangles, evalPoints, u, pre, compute);

    return u;
  }
};
} // namespace m2

#endif