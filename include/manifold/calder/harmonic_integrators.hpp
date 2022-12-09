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

#include "manifold/asawa/coordinate_interface.hpp"
#include "manifold/asawa/geometry_types.hpp"
#include "manifold/asawa/m2.hpp"

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

  /*
    ////////////////////////////////////////////////////////////////////////////
    // Point covariance
    ////////////////////////////////////////////////////////////////////////////

    vector<mat43> covariance(std::vector<coordinate_type> evalPoints,
                             std::vector<T> regLength) {
      using Integrator = asawa::Simple_BarnesHutt<SPACE, real, mat3>;
      using Tree = typename Integrator::Tree;
      using Node = typename Integrator::Node;

      real dummy = 0.0;
      auto pre = [](const vector<real> &charges, const vector<T> &weights,
                    const vector<coordinate_type> &points, Node &node, Tree
  &tree, coordinate_type &avgPoint, real &netCharge) -> void { avgPoint =
  coordinate_type(0, 0, 0); T netWeight = 0;

        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          avgPoint += weights[ii] * points[ii];
          netWeight += weights[ii];
        }
        // netCharge /= netWeight;
        avgPoint /= netWeight;
      };

      auto computeK = [](T dist, T C) {
  #if 0
        T d3 = dist * dist * dist;
        T l3 = C * C * C;
        T kappa = (1.0 - exp(-d3 / l3)) / d3;
        return kappa / (4.0 * M_PI);
  #elif 1
        T dist3 = dist * dist * dist;
        T l3 = C * C * C;
        T kappa = 1.0 / (dist3 + l3);
        return kappa / pow(4.0 * M_PI, 1.5);
  #endif
        // return 1 / (dist*dist + C*C);
      };

      auto compute = [this, regLength, computeK](
                         int i_c, const real &ci, const real &cj,
                         const coordinate_type &pi, const coordinate_type &pj,
                         const vector<triangle_type> &tris, Node &node,
                         Tree &tree) -> mat3 {
        mat3 out = z::zero<mat3>();
        T reg = regLength[i_c];
        coordinate_type dp = pj - pi;

        T dist = va::norm(dp);
        T k = computeK(dist, reg);
        out = k * va::outer<real>(dp, dp);
        return out;
      };

      vector<face_ptr> &faces = mesh->get_faces();
      vector<triangle_type> triangles;
      for (int i = 0; i < faces.size(); i++) {
        if (!mesh->has_face(i))
          continue;
        if (faces[i]->size() < 3)
          continue;
        std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
        triangles.insert(triangles.end(), tris.begin(), tris.end());
      }

      vector<mat3> u0(evalPoints.size(), z::zero<mat3>());
      Integrator integrator;
      integrator.integrate(dummies, triangles, evalPoints, u0, pre, compute);

      vector<mat43> u1(evalPoints.size(), z::zero<mat43>());
      int i = 0;
      for (const auto &ui : u0) {
        Eigen::JacobiSVD<mat3> svd(ui, Eigen::ComputeFullU);
        const mat3 U = svd.matrixU();
        const coordinate_type S = svd.singularValues();
        u1[i].row(0) = U.row(0);
        u1[i].row(1) = U.row(1);
        u1[i].row(2) = U.row(2);
        u1[i++].row(3) = S;
      }
      return u1;
    }
  */

  ////////////////////////////////////////////////////////////////////////////
  // surface covariance
  ////////////////////////////////////////////////////////////////////////////

  vector<mat43> covariance(surf_ptr mesh,
                           const std::vector<coordinate_type> &evalPoints,
                           const std::vector<T> &regLength) {
    using Avg_Integrator =
        Geometry_Integrator<SPACE, real, triangle_type, mat3>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;
    std::vector<real> dummies(mesh->get_faces().size(), z::zero<real>());

    real dummy = 0.0;
    auto pre = [](const vector<triangle_type> &tris, ANode &node, ATree &tree,
                  real &netCharge, coordinate_type &avgPoint,
                  coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
        netWeight += w;
      }
      // netCharge /= netWeight;
      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
#if 0
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / (4.0 * M_PI);
#elif 1
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3 + l3);
      return kappa / pow(4.0 * M_PI, 1.5);
#endif
      // return 1 / (dist*dist + C*C);
    };

    auto compute = [this, regLength, computeK](
                       int i_c, const real &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> mat3 {
      mat3 out = z::zero<mat3>();
      T reg = regLength[i_c];
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];

          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type Ni = tri.normal();
          coordinate_type dp = pe - c;
          T dist = va::norm(dp);
          T k = computeK(dist, reg);
          // T dotN = va::dot(Ni, dp);
          out += w * k * va::outer<real>(dp, dp);
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pe - pc;
        T dist = va::norm(dp);
        // T Ndp = va::dot(N, dp);
        T k = computeK(dist, reg);
        // T dotN = va::dot(N, dp);
        out += k * va::outer<real>(dp, dp);
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      if (!mesh->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    vector<mat3> u0(evalPoints.size(), z::zero<mat3>());
    Avg_Integrator integrator;
    integrator.integrate(dummies, triangles, evalPoints, u0, pre, compute);

    vector<mat43> u1(evalPoints.size(), z::zero<mat43>());
    int i = 0;
    for (const auto &ui : u0) {
      Eigen::JacobiSVD<mat3> svd(ui, Eigen::ComputeFullU);
      const mat3 U = svd.matrixU();
      const coordinate_type S = svd.singularValues();
      u1[i].row(0) = U.row(0);
      u1[i].row(1) = U.row(1);
      u1[i].row(2) = U.row(2);
      u1[i++].row(3) = S;
    }
    return u1;
  }

  vector<mat43> covariance(surf_ptr mesh,
                           std::vector<coordinate_type> evalPoints,
                           T regLength = 0.5) {
    std::vector<T> regs(evalPoints.size(), regLength);
    return covariance(mesh, evalPoints, regs);
  }
  ////////////////////////////////////////////////////////////////////////////
  // covariance
  ////////////////////////////////////////////////////////////////////////////

  vector<mat43> curvature(surf_ptr mesh,
                          const std::vector<coordinate_type> &evalPoints,
                          std::vector<T> regLength) {
    using Avg_Integrator =
        Geometry_Integrator<SPACE, real, triangle_type, mat3>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;
    std::vector<real> dummies(mesh->get_faces().size(), z::zero<real>());

    real dummy = 0.0;
    auto pre = [](const vector<triangle_type> &tris, ANode &node, ATree &tree,
                  real &netCharge, coordinate_type &avgPoint,
                  coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
        netWeight += w;
      }
      // netCharge /= netWeight;
      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
#if 1
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / (4.0 * M_PI);
#elif 0
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3 + l3);
      return kappa / pow(4.0 * M_PI, 1.5);
#endif
      // return 1 / (dist*dist + C*C);
    };

    auto compute = [this, regLength, computeK](
                       int i_c, const real &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> mat3 {
      mat3 out = z::zero<mat3>();
      T reg = regLength[i_c];
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];

          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type Ni = tri.normal();
          coordinate_type dp = pe - c;
          coordinate_type cN = va::cross(Ni, dp);
          T dist = va::norm(dp);
          T k = computeK(dist, reg);
          // T dotN = va::dot(Ni, dp);
          // out += w * k * va::outer<real>(dp, cN);
          mat3 Ndp = va::outer<real>(Ni, dp);
          out += w * k * (Ndp + Ndp.transpose());
          // out += w * k * va::outer<real>(dp, dp) *
          // va::outer<real>(Ni, Ni);
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pe - pc;
        T dist = va::norm(dp);
        // T Ndp = va::dot(N, dp);
        T k = computeK(dist, reg);
        coordinate_type cN = va::cross(N, dp);
        // T dotN = va::dot(N, dp);

        // out += k * va::outer<real>(dp, cN);
        mat3 Ndp = va::outer<real>(N, dp);
        out += k * (Ndp + Ndp.transpose());
        // out += k * va::outer<real>(dp, dp) * va::outer<real>(N, N);
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      if (!mesh->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    vector<mat3> u0(evalPoints.size(), z::zero<mat3>());
    Avg_Integrator integrator;
    integrator.integrate(dummies, triangles, evalPoints, u0, pre, compute);

    vector<mat43> u1(evalPoints.size(), z::zero<mat43>());
    int i = 0;
    for (const auto &ui : u0) {
      Eigen::JacobiSVD<mat3> svd(ui, Eigen::ComputeFullU);
      const mat3 U = svd.matrixU();
      const coordinate_type S = svd.singularValues();
      u1[i].row(0) = U.row(0);
      u1[i].row(1) = U.row(1);
      u1[i].row(2) = U.row(2);
      u1[i++].row(3) = S;
    }
    return u1;
  }

  vector<mat43> curvature(surf_ptr mesh,
                          std::vector<coordinate_type> evalPoints,
                          T regLength = 0.5) {
    std::vector<T> regs(evalPoints.size(), regLength);
    return curvature(mesh, evalPoints, regs);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Sphere
  ////////////////////////////////////////////////////////////////////////////
  using ls_sphere_t = ls_sphere<real, coordinate_type>;

  vector<ls_sphere_t>
  least_squares_sphere(surf_ptr mesh,
                       const std::vector<coordinate_type> &evalPoints,
                       T regLength = 0.5) {
    using Avg_Integrator =
        Geometry_Integrator<SPACE, real, triangle_type, ls_sphere_t>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;

    if (evalPoints.empty())
      return std::vector<ls_sphere_t>();

    std::vector<real> dummies(mesh->get_faces().size(), z::zero<real>());

    real dummy = 0.0;
    auto pre = [](const vector<triangle_type> &tris, ANode &node, ATree &tree,
                  real &netCharge, coordinate_type &avgPoint,
                  coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
        netWeight += w;
      }
      // netCharge /= netWeight;
      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
#if 0
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / (4.0 * M_PI);
#elif 1
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3 + l3);
      return kappa / pow(4.0 * M_PI, 1.5);
#endif
      // return 1 / (dist*dist + C*C);
    };

    auto compute = [this, regLength, computeK](
                       int i_c, const real &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> ls_sphere_t {
      ls_sphere_t out;
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];

          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type Ni = tri.normal();
          coordinate_type dp = pe - c;
          T dist = va::norm(dp);
          T k = computeK(dist, regLength);
          // T dotN = va::dot(Ni, dp);
          out.add_vec(-w * k * dp, w * k);
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pe - pc;
        T dist = va::norm(dp);
        // T Ndp = va::dot(N, dp);
        T k = computeK(dist, regLength);
        // T dotN = va::dot(N, dp);
        out.add_vec(-k * dp, k);
      }

      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      if (!mesh->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    vector<ls_sphere_t> u0(evalPoints.size());
    Avg_Integrator integrator;
    integrator.integrate(dummies, triangles, evalPoints, u0, pre, compute);
    return u0;
  }

  ////////////////////////////////////////////////////////////////////////////
  // AVG
  ////////////////////////////////////////////////////////////////////////////

  template <typename Q>
  vector<Q> harmonicAvg(surf_ptr mesh, std::vector<Q> faceVectors,
                        const std::vector<coordinate_type> &evalPoints,
                        T regLength = 0.5) {
    using Avg_Integrator = Geometry_Integrator<SPACE, Q, triangle_type, Q>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;

    if (evalPoints.empty())
      return std::vector<Q>();

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
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3 + l3);
      return kappa / pow(4.0 * M_PI, 1.5);
    };

    coordinate_array vertex_normals =
        asawa::ci::get_vertex_normals<SPACE>(mesh);
    auto compute = [this, &vertex_normals, faceVectors, regLength, computeK](
                       int i_c, const Q &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> Q {
      Q out = z::zero<Q>();

      coordinate_type vn = vertex_normals[i_c];

      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          Q qi = faceVectors[ii];
          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type dp = c - pe;
          T dist = va::norm(dp);
          T k = computeK(dist, regLength);
          out += w * k * qi;
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pc - pe;
        T dist = va::norm(dp);
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
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
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

  ////////////////////////////////////////////////////////////////////////////
  // AVG
  ////////////////////////////////////////////////////////////////////////////

  template <typename Q>
  vector<Q> gaussAvg(surf_ptr mesh, std::vector<Q> faceVectors,
                     const std::vector<coordinate_type> &evalPoints,
                     T regLength = 0.5) {
    using Avg_Integrator = Geometry_Integrator<SPACE, Q, triangle_type, Q>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;

    if (evalPoints.empty())
      return std::vector<Q>();

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
      T dist2 = dist * dist;
      T kappa = exp(-dist2 / C);
      return kappa / C / pow(2.0 * M_PI, 0.5);
    };

    coordinate_array vertex_normals =
        asawa::ci::get_vertex_normals<SPACE>(mesh);
    auto compute = [this, &vertex_normals, faceVectors, regLength, computeK](
                       int i_c, const Q &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> Q {
      Q out = z::zero<Q>();

      coordinate_type vn = vertex_normals[i_c];

      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          Q qi = faceVectors[ii];
          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type dp = c - pe;
          T dist = va::norm(dp);
          T k = computeK(dist, regLength);
          out += w * k * qi;
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pc - pe;
        T dist = va::norm(dp);
        // T Ndp = va::dot(N, dp);
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
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
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

  ////////////////////////////////////////////////////////////////////////////
  // Divergence
  ////////////////////////////////////////////////////////////////////////////

  template <typename G, typename Q>
  vector<G> divergence(surf_ptr mesh, std::vector<Q> faceVectors,
                       std::vector<coordinate_type> evalPoints,
                       T regLength = 0.5) {
    using Avg_Integrator = Geometry_Integrator<SPACE, Q, triangle_type, G>;
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
      // netCharge /= netWeight;
      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
#if 1
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / (4.0 * M_PI);
      // return kappa;

#elif 0
      T d1 = dist;
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T den = (d3 + l3) * (d3 + l3);
      T kappa = -3.0 * d1 / den;
      return kappa / pow(4.0 * M_PI, 1.5);
#endif
    };

    auto compute = [this, &faceVectors, regLength, computeK](
                       int i_c, const Q &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> G {
      G out = z::zero<G>();
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          const Q &qi = faceVectors[ii];

          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type dp = pe - c;
          T dist = va::norm(dp);
          T k = computeK(dist, regLength);
          out += w * k * va::dot<G>(dp, qi);
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pe - pc;
        T dist = va::norm(dp);
        // T Ndp = va::dot(N, dp);
        T k = computeK(dist, regLength);

        out += k * va::dot<G>(dp, wq);
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
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }

    vector<G> u(evalPoints.size(), z::zero<G>());
    Avg_Integrator integrator;
    integrator.integrate(faceVectors, triangles, evalPoints, u, pre, compute);

    return u;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Grad1
  ////////////////////////////////////////////////////////////////////////////

  template <typename G, typename Q>
  vector<G> gradient(surf_ptr mesh, std::vector<Q> faceVectors,
                     std::vector<coordinate_type> evalPoints,
                     T regLength = 0.5) {
    using Avg_Integrator = Geometry_Integrator<SPACE, Q, triangle_type, G>;
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
      // netCharge /= netWeight;
      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
#if 1
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / (4.0 * M_PI);
      // return kappa;

#elif 0
      T d1 = dist;
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T den = (d3 + l3) * (d3 + l3);
      T kappa = -3.0 * d1 / den;
      return kappa / pow(4.0 * M_PI, 1.5);
#endif
    };

    auto compute = [this, &faceVectors, regLength, computeK](
                       int i_c, const Q &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> G {
      G out = z::zero<G>();
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          Q qi = faceVectors[ii];

          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type dp = pe - c;
          T dist = va::norm(dp);
          T k = computeK(dist, regLength);
          out += w * k * dp * qi;
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pe - pc;
        T dist = va::norm(dp);
        // T Ndp = va::dot(N, dp);
        T k = computeK(dist, regLength);

        out += k * dp * wq;
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
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }

    vector<G> u(evalPoints.size(), z::zero<G>());
    Avg_Integrator integrator;
    integrator.integrate(faceVectors, triangles, evalPoints, u, pre, compute);

    return u;
  }
  ////////////////////////////////////////////////////////////////////////////
  // Curl
  ////////////////////////////////////////////////////////////////////////////

  template <typename Q>
  vector<Q> curl(surf_ptr mesh, std::vector<Q> faceVectors,
                 std::vector<coordinate_type> evalPoints, T regLength = 0.5) {
    using Avg_Integrator = Geometry_Integrator<SPACE, Q, triangle_type, Q>;
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
#if 1
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 0
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
                       int i_c, const Q &wq, const coordinate_type &pc,
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
          T dist = va::norm(dp);
          T k = computeK(dist, regLength);

          out += w * k * va::cross(qi, dp);
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pc - pe;
        T dist = va::norm(dp);
        // T Ndp = va::dot(N, dp);
        T k = computeK(dist, regLength);
        out += k * va::cross(wq, dp);
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
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
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

  ////////////////////////////////////////////////////////////////////////////
  // messing around
  ////////////////////////////////////////////////////////////////////////////

  template <typename G, typename Q>
  vector<G> gravitation(surf_ptr mesh, std::vector<Q> faceVals,
                        std::vector<coordinate_type> evalPoints,
                        std::vector<Q> evalVals, T regLength = 0.5) {
    using Avg_Integrator = Geometry_Integrator<SPACE, Q, triangle_type, G>;
    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;
    coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(mesh);
    auto pre = [faceVals](const vector<triangle_type> &tris, ANode &node,
                          ATree &tree, Q &netCharge, coordinate_type &avgPoint,
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
        netCharge += w * faceVals[ii];
        netWeight += w;
      }
      // netCharge /= netWeight;
      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3 + l3);
      return kappa / pow(4.0 * M_PI, 1.5);
    };

    auto compute = [this, &faceVals, &evalVals, &normals, regLength, computeK](
                       int i_c, const Q &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> G {
      G out = z::zero<G>();
      Q qc = evalVals[i_c];
      coordinate_type Nv = normals[i_c];
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          Q qi = faceVals[ii];
          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type Nf = tri.normal();
          coordinate_type dp = c - pe;
          T dist = va::norm(dp);

          // if(dist <= std::numeric_limits<T>::epsilon())
          //   continue;
          //  out += w * computeK(dist, regLength) * va::cross(q,dp);
          T k = computeK(dist, regLength);
          // T Ndp = tri.solidAngle(pe);
          T NdN = 0.5 - 0.5 * va::dot(Nv, Nf);
          out += w * k * dp * NdN * qc * qi;
        }
      } else {
        // out += computeK(dist, regLength) * va::cross(q, dp);
        coordinate_type dp = pc - pe;
        T dist = va::norm(dp);
        T k = computeK(dist, regLength);
        T Ndp = va::dot(N, dp);
        T NdN = 0.5 - 0.5 * va::dot(Nv, N.normalized());
        out += k * dp * NdN * qc * wq;
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();

    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      if (!mesh->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }
    vector<G> u(evalPoints.size(), z::zero<G>());
    Avg_Integrator integrator;
    integrator.integrate(faceVals, triangles, evalPoints, u, pre, compute);

    return u;
  }
};

template <typename SPACE>
vector<typename SPACE::vec3>
calcAvgVerts(asawa::surf<SPACE> *mesh,
             const std::vector<typename SPACE::vec3> &vertVals,
             const std::vector<typename SPACE::coordinate_type> &evalPoints,
             typename SPACE::real regLength = 0.5) {
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  std::vector<vec3> faceVals = asawa::ci::verts_to_faces<SPACE>(vertVals, mesh);

  mesh_calculator<SPACE> calc;
  std::vector<vec3> u =
      // calc.template gaussAvg<vec3>(mesh, faceVals, evalPoints, regLength);
      calc.template harmonicAvg<vec3>(mesh, faceVals, evalPoints, regLength);

  return u;
}
} // namespace calder

#endif