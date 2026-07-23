#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/vec_addendum.h"

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
#include "gaudi/geometry_logger.hpp"

#ifndef __WALKY_TALKY__
#define __WALKY_TALKY__

namespace gaudi {
namespace duchamp {

using namespace asawa;

// Knobs for silly_walk / rotate_walk / align_walk.
// Aggregate — use designated initializers:
//   silly_walk_config{.i0 = 0, .N_steps = 10000, .rotate = true,
//                     .twist_amp = 12.0, .twist_freq = -6.0};
struct silly_walk_config {
  // Starting half-edge / corner index on the shell.
  index_t i0 = 0;
  // Maximum number of geodesic steps before the walk stops.
  index_t N_steps = 4000;
  // Initial heading in the edge tangent–binormal plane:
  //   dir = cos(thet) * T + sin(thet) * B
  real thet = 0.0;

  // --- rotate_walk: spin the heading about the local edge normal -------------
  // Enabled when rotate == true. Each step applies
  //   AngleAxis(twist_amp * cos(twist_freq * a) * li * π, N)
  // where a is the local edge angle and li is the last step length.
  bool rotate = false;
  // Scale of the per-step twist about the edge normal.
  real twist_amp = 5.0;
  // Spatial frequency of twist modulation vs. local edge angle a.
  real twist_freq = -8.0;

  // --- align_walk: steer toward PCA of nearby prior walk tangents ------------
  // Enabled when align == true. Directions are analyzed in the local tangent
  // plane (2×2 PCA in {T,B}); heading becomes:
  //   dir ∝ ca[0]*d0 + ca[1]*±ℓ₊*f₊ + ca[2]*±ℓ₋*f₋
  // where f₊/f₋ are the major/minor in-plane PCA axes and ℓ± their
  // normalized eigenvalues.
  //   ca[0]  weight on the current direction d0
  //   ca[1]  weight on the major PCA axis (largest eigenvalue)
  //   ca[2]  weight on the minor PCA axis (smallest eigenvalue)
  bool align = false;
  vec3 ca = vec3(1.0, 0.5, 0.0);
  // Neighborhood radius for align_walk's Gaussian (NOT the walk geometric eps).
  real align_eps = 1.0e-1;

  // Geometric tolerance for the underlying shell walk only.
  real eps = 1.0e-8;

  std::vector<vec3> run(const shell::shell &M) const;
};

vec3 align_walk(const vec3 x0, const vec3 &d0, const vec3 &N0, real li,
                const std::vector<vec3> walk, const std::vector<vec3> &normals,
                real eps = 1e-1, vec3 C = vec3(1.0, 0.5, 0.0)) {
  // PCA of nearby walk tangents in the local tangent plane at N0.
  if (walk.size() < 16)
    return d0;

  const vec3 N = N0.normalized();
  vec3 T = va::reject(N, d0);
  if (T.norm() < 1e-12) {
    // d0 nearly normal — pick any tangent basis.
    T = va::reject(N, vec3(1.0, 0.0, 0.0));
    if (T.norm() < 1e-12)
      T = va::reject(N, vec3(0.0, 1.0, 0.0));
  }
  T.normalize();
  const vec3 B = N.cross(T).normalized();

  mat2 S = mat2::Zero();
  real w = 0.0;
  for (int i = 1; i < static_cast<int>(walk.size()) - 1; i++) {
    const vec3 xi1 = walk[i];
    const vec3 xi0 = walk[i - 1];
    const vec3 Ni = normals[i];
    Eigen::Quaterniond q;
    q.setFromTwoVectors(Ni.normalized(), N);

    vec3 di = q * (xi1 - xi0);
    di = va::reject(N, di);
    if (di.norm() < 1e-16)
      continue;

    const real d = (x0 - 0.5 * (xi1 + xi0)).norm();
    const real wi = std::exp(-d * d / eps / eps);
    w += wi;

    const vec2 uv(di.dot(T), di.dot(B));
    S += wi * (uv * uv.transpose());
  }
  if (w < 1e-16)
    return d0;
  S /= w;

  Eigen::SelfAdjointEigenSolver<mat2> es(S);
  if (es.info() != Eigen::Success || es.eigenvalues().hasNaN())
    return d0;

  // Eigenvalues ascending: col(0) = minor, col(1) = major.
  const real ssum = es.eigenvalues().sum();
  if (ssum < 1e-16)
    return d0;
  const real l_minus = es.eigenvalues()[0] / ssum;
  const real l_plus = es.eigenvalues()[1] / ssum;
  const vec2 e_minus = es.eigenvectors().col(0);
  const vec2 e_plus = es.eigenvectors().col(1);

  const vec3 f_minus = e_minus[0] * T + e_minus[1] * B;
  const vec3 f_plus = e_plus[0] * T + e_plus[1] * B;

  vec3 dir = C[0] * d0;
  dir += C[1] * va::sgn(d0, f_plus) * l_plus * f_plus;
  dir += C[2] * va::sgn(d0, f_minus) * l_minus * f_minus;
  dir = va::reject(N, dir);
  if (dir.norm() < 1e-16)
    return d0;
  dir.normalize();
  return dir;
}

vec3 rotate_walk(const shell::shell &M, shell::CornerId ci, const vec3 &d0,
                 const vec3 &N0, real li, vec2 C) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  real a = asawa::shell::angle(M, ci, x);
  return Eigen::AngleAxis<real>(C[0] * cos(C[1] * a) * li * M_PI, N0) * d0;
}

std::vector<vec3> walk(const shell::shell &M, const real &thet = 0.0,
                       shell::CornerId c0 = shell::corner_id(0),
                       const index_t &N_steps = 4000,
                       real eps = 1.0e-8) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> &v = asawa::const_get_vec_data(M, 1);

  vec3 N = asawa::shell::edge_normal(M, c0, x);
  vec3 T = asawa::shell::edge_tangent(M, c0, x).normalized();
  vec3 B = N.cross(T).normalized();

  vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
  std::vector<index_t> corners;
  std::vector<real> S;
  std::vector<vec3> points;
  std::vector<vec3> normals;
  real l = 0.0;

  asawa::shell::walk(M, x, c0, dir, 0.5, N_steps, eps,
                     [&](const asawa::shell::shell &M,
                         const std::vector<vec3> &x,
                         const asawa::shell::CornerId &ci,
                         const real &s, const real &accumulated_length,
                         vec3 &dir) {
                       S.push_back(s);
                       corners.push_back(ci);
                       vec3 pt = asawa::shell::edge_vert(M, ci, s, x);
                       vec3 Ni = asawa::shell::edge_normal(M, ci, x);
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
                             shell::CornerId c0 = shell::corner_id(0),
                             const index_t &N_steps = 4000,               //
                             bool rotate = false, vec2 cr = vec2::Zero(), //
                             bool align = false, vec3 ca = vec3::Zero(),  //
                             real eps = 1.0e-8, real align_eps = 1.0e-1) {
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> &v = asawa::const_get_vec_data(M, 1);

  vec3 N = asawa::shell::edge_normal(M, c0, x);
  vec3 T = asawa::shell::edge_tangent(M, c0, x).normalized();
  vec3 B = N.cross(T).normalized();

  vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
  std::vector<index_t> corners;
  std::vector<real> S;
  std::vector<vec3> points;
  std::vector<vec3> normals;
  real l = 0.0;

  asawa::shell::walk(
      M, x, c0, dir, 0.5, N_steps, eps,
      [&](const asawa::shell::shell &M, const std::vector<vec3> &x,
          const asawa::shell::CornerId &ci, const real &s,
          const real &accumulated_length, vec3 &dir) {
        S.push_back(s);
        corners.push_back(ci);
        vec3 pt = asawa::shell::edge_vert(M, ci, s, x);
        vec3 Ni = asawa::shell::edge_normal(M, ci, x);
        real li = 0.0;
        if (points.size() > 0)
          li = (pt - points.back()).norm();
        if (rotate)
          dir = rotate_walk(M, ci, dir, Ni, li, cr);
        if (align)
          dir = align_walk(pt, dir, Ni, li, points, normals, align_eps, ca);

        points.push_back(pt);
        normals.push_back(Ni);
        return true;
      });
  std::cout << "walk resulted in: " << points.size() << " points" << std::endl;
  return points;
}

inline std::vector<vec3> silly_walk_config::run(const shell::shell &M) const {
  return silly_walk(M, thet, shell::corner_id(i0), N_steps, rotate,
                    vec2(twist_amp, twist_freq), align, ca, eps, align_eps);
}

} // namespace duchamp
} // namespace gaudi
#endif
