#ifndef GAUDI_KUSAMA_TWISTED_LAPLACIAN_HPP
#define GAUDI_KUSAMA_TWISTED_LAPLACIAN_HPP

/// Cotan Laplacian on one complex value per vertex (interleaved `(Re, Im)`)
/// whose edges carry a unit twist: the quadratic form
///   sum over edges ij of w_ij |z_j - e^{i a_ij} C_ij z_i|^2,
/// `w_ij` the cotan weight, `a_ij` the edge's angle and `C_ij` the identity
/// or (when the edge's sign is -1) conjugation. Stripes twist by a phase
/// advance (`kusama::stripes`); a direction field would twist by the
/// transport angle between vertex tangent charts (Knöppel 2013). The mass
/// is the lumped vertex area on both components.

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"

#include <Eigen/Sparse>
#include <cmath>
#include <vector>

namespace gaudi {
namespace kusama {
namespace twisted_laplacian {

struct twist {
  real angle = 0;
  /// +1, or -1 to conjugate `z_i` across the edge.
  real sign = 1;
};

/// Assembles `A` and the diagonal mass `B`. `edge(i, j)` returns the twist
/// carrying `z_i` to `z_j` for one orientation of each edge; triangles only.
template <typename EdgeFn>
inline void assemble(const asawa::shell::shell &M, const std::vector<vec3> &x, EdgeFn edge,
                     matS &A, vecX &B) {
  const int nV = static_cast<int>(M.vert_count());
  B = vecX::Zero(2 * nV);
  std::vector<Eigen::Triplet<real>> t;
  t.reserve(static_cast<size_t>(M.face_count()) * 36);
  for (asawa::shell::FaceId f : M.get_face_range()) {
    if (M.fsize(f) != 3)
      continue;
    const auto tri = M.get_tri(f);
    const int v[3] = {static_cast<int>(tri[0]), static_cast<int>(tri[1]),
                      static_cast<int>(tri[2])};
    const real area = real(0.5) * (x[v[1]] - x[v[0]]).cross(x[v[2]] - x[v[0]]).norm();
    for (int c = 0; c < 3; ++c) {
      B[2 * v[c]] += area / real(3);
      B[2 * v[c] + 1] += area / real(3);
    }
    // Each corner's half cotangent weights the edge opposite it.
    for (int c = 0; c < 3; ++c) {
      const int k = v[c], i = v[(c + 1) % 3], j = v[(c + 2) % 3];
      const vec3 a = x[i] - x[k], b = x[j] - x[k];
      const real cr = a.cross(b).norm();
      if (!(cr > real(1e-20)))
        continue;
      const real w = real(0.5) * a.dot(b) / cr;
      const twist tw = edge(i, j);
      // w (|z_i|^2 + |z_j|^2 - 2 z_j . R C z_i), R the rotation by the angle.
      const real co = std::cos(tw.angle), si = std::sin(tw.angle);
      const real RC[2][2] = {{co, -si * tw.sign}, {si, co * tw.sign}};
      for (int d = 0; d < 2; ++d) {
        t.emplace_back(2 * i + d, 2 * i + d, w);
        t.emplace_back(2 * j + d, 2 * j + d, w);
      }
      for (int r = 0; r < 2; ++r)
        for (int q = 0; q < 2; ++q) {
          t.emplace_back(2 * j + r, 2 * i + q, -w * RC[r][q]);
          t.emplace_back(2 * i + q, 2 * j + r, -w * RC[r][q]);
        }
    }
  }
  A.resize(2 * nV, 2 * nV);
  A.setFromTriplets(t.begin(), t.end());
}

} // namespace twisted_laplacian
} // namespace kusama
} // namespace gaudi

#endif
