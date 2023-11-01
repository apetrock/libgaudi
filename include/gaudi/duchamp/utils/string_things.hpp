#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/walk.hpp"
#include "gaudi/common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __WALKY_TALKY__
#define __WALKY_TALKY__

namespace gaudi {
namespace duchamp {

using namespace asawa;

vec3 align_walk(const vec3 x0, const vec3 &d0, const vec3 &N0, real li,
                const std::vector<vec3> walk, const std::vector<vec3> &normals,
                real eps = 1e-1, vec4 C = vec4(5.0, -8.0, 0.0, 0.0)) {
  // dumb little test to see if I can align the walk to neighboring lines...
  // we'll do this N^2 for fun
  if (walk.size() < 16)
    return d0;

  mat3 M = mat3::Zero();
  mat3 R = va::rejection_matrix(N0);
  real w = 0.0;
  for (int i = 1; i < walk.size() - 1; i++) {
    vec3 xi1 = walk[i];
    vec3 xi0 = walk[i - 1];
    vec3 Ni = normals[i];
    Eigen::Quaterniond q;
    q.setFromTwoVectors(Ni.normalized(), N0.normalized());

    // gg::geometry_logger::line(x0, x0 + 0.1 * N0, vec4(0.0, 0.0, 1.0, 1.0));
    // gg::geometry_logger::line(x0, x0 + 0.1 * Ni, vec4(0.0, 0.5, 1.0, 1.0));

    vec3 di = xi1 - xi0;
    di = q * di;
    real d = (x0 - 0.5 * (xi1 + xi0)).norm();

    real c2 = std::exp(-d * d / eps / eps);
    real s = va::sgn(d0.dot(di));
    real wi = c2;
    w += wi;
    // M += wi * (d0 * di.transpose() + di * d0.transpose());
    M += wi * di * di.transpose();
  }
  M /= w;

  Eigen::SelfAdjointEigenSolver<mat3> es(M);
  if (es.eigenvalues().minCoeff() < 1e-16) {
    std::cout << "SelfAdjointEigenSolver singular!" << std::endl;
    // matrix is singular
  }

  if (es.eigenvalues().hasNaN())
    return d0;

  real ssum = es.eigenvalues().sum();
  real l0 = es.eigenvalues()[0] / ssum;
  vec3 f0 = es.eigenvectors().col(0);
  real l1 = es.eigenvalues()[1] / ssum;
  vec3 f1 = es.eigenvectors().col(1);
  real l2 = es.eigenvalues()[2] / ssum;
  vec3 f2 = es.eigenvectors().col(2);

  vec3 dir = C[0] * d0;
  dir += C[1] * va::sgn(d0, f0) * l0 * f0;
  dir += C[2] * va::sgn(d0, f1) * l1 * f1;
  dir += C[3] * va::sgn(d0, f2) * l2 * f2;
  dir.normalize();
  return dir;
}

vec3 rotate_walk(const shell::shell &M, const index_t &ci, const vec3 &d0,
                 const vec3 &N0, real li, vec2 C) {
  // dumb little test to see if I can align the walk to neighboring lines...
  // we'll do this N^2 for fun

  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  real a = asawa::shell::angle(M, ci, x);
  return Eigen::AngleAxis<real>(C[0] * cos(C[1] * a) * li * M_PI, N0) * d0;
}

std::vector<vec3> walk(const shell::shell &M, const real &thet = 0.0,
                       const index_t i0 = 0, const index_t &N_steps = 4000,
                       real eps = 1.0e-8) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> &v = asawa::const_get_vec_data(M, 1);

  vec3 N = asawa::shell::edge_normal(M, i0, x);
  vec3 T = asawa::shell::edge_tangent(M, i0, x).normalized();
  vec3 B = N.cross(T).normalized();

  vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
  std::vector<index_t> corners;
  std::vector<real> S;
  std::vector<vec3> points;
  std::vector<vec3> normals;
  real l = 0.0;

  asawa::shell::walk(M, x, i0, dir, 0.5, N_steps, eps,
                     [&](const asawa::shell::shell &M,
                         const std::vector<vec3> &x, const index_t &corner,
                         const real &s, const real &accumulated_length,
                         vec3 &dir) {
                       asawa::shell::index_t ci = corner;
                       S.push_back(s);
                       corners.push_back(corner);
                       vec3 pt = asawa::shell::edge_vert(M, ci, s, x);
                       vec3 Ni = asawa::shell::edge_normal(M, corner, x);
                       real li = 0.0;
                       if (points.size() > 0)
                         li = (pt - points.back()).norm();
                       points.push_back(pt);
                       normals.push_back(Ni);
                       return true;
                     });
  std::cout << "walk resulted in: " << points.size() << " points" << std::endl;
  return points;
}

std::vector<vec3> silly_walk(const shell::shell &M, const real &thet = 0.0,
                             const index_t i0 = 0,
                             const index_t &N_steps = 4000,               //
                             bool rotate = false, vec2 cr = vec2::Zero(), //
                             bool align = false, vec4 ca = vec4::Zero(),  //
                             real eps = 1.0e-8) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> &v = asawa::const_get_vec_data(M, 1);

  vec3 N = asawa::shell::edge_normal(M, i0, x);
  vec3 T = asawa::shell::edge_tangent(M, i0, x).normalized();
  vec3 B = N.cross(T).normalized();

  vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
  std::vector<index_t> corners;
  std::vector<real> S;
  std::vector<vec3> points;
  std::vector<vec3> normals;
  real l = 0.0;

  asawa::shell::walk(
      M, x, i0, dir, 0.5, N_steps, 1e-8,
      [&](const asawa::shell::shell &M, const std::vector<vec3> &x,
          const index_t &corner, const real &s, const real &accumulated_length,
          vec3 &dir) {
        asawa::shell::index_t ci = corner;
        S.push_back(s);
        corners.push_back(corner);
        vec3 pt = asawa::shell::edge_vert(M, ci, s, x);
        vec3 Ni = asawa::shell::edge_normal(M, corner, x);
        real li = 0.0;
        if (points.size() > 0)
          li = (pt - points.back()).norm();
        if (rotate)
          dir = rotate_walk(M, corner, dir, Ni, li, cr);
        if (align)
          dir = align_walk(pt, dir, Ni, li, points, normals, eps, ca);

        points.push_back(pt);
        normals.push_back(Ni);
        return true;
      });
  std::cout << "walk resulted in: " << points.size() << " points" << std::endl;
  return points;
}
} // namespace duchamp
} // namespace gaudi
#endif