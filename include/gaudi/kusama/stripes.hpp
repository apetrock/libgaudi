#ifndef GAUDI_KUSAMA_STRIPES_HPP
#define GAUDI_KUSAMA_STRIPES_HPP

/// Stripe patterns (Knöppel, Crane, Pinkall, Schröder 2015,
/// `reference/knoeppel2015-stripe-patterns.pdf`), the core solve only: one
/// complex value `psi` per vertex whose phase advances at frequency `nu`
/// along a unit tangent field `X`, so its level lines run across `X` and
/// close up wherever the pattern has no zero. The energy (their Eq. 9,
/// extrinsic form of Alg. 3–4) is
///   sum over edges ij of w_ij |psi_j - e^{i omega_ij} conj^?(psi_i)|^2,
///   omega_ij = nu / 2 <x_j - x_i, X_i + s_ij X_j>,
/// with cotan weights `w_ij`, and `s_ij = sign <X_i, X_j>` conjugating across
/// edges where the field flips (a line field). The pattern is its smallest
/// generalized eigenvector against the lumped vertex mass, by inverse power
/// iteration (Alg. 6) on one factorization; warm starts from the last
/// pattern converge in a few backsolves. Texture coordinates and the
/// paper's rendering of zeros and branch triangles are not here: callers
/// read the phase gradient (`level_directions`).

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/kusama/twisted_laplacian.hpp"
#include "gaudi/sparse_solver.h"

#include <Eigen/Sparse>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <vector>

namespace gaudi {
namespace kusama {
namespace stripes {

/// Energy `A` and lumped mass `B` (diagonal) on interleaved `(Re, Im)` per
/// vertex, plus the edge data the phase read-out needs. `A` is shifted by
/// `shift * B`, which moves every eigenvalue by `shift` and keeps the factor
/// definite when the field is exactly integrable.
struct system {
  matS A;
  vecX B;
  real nu = 0;
  real shift = 0;
  std::vector<vec3> X;
};

/// Phase step and sign of edge `i -> j`.
inline real edge_omega(const std::vector<vec3> &x, const std::vector<vec3> &X, int i, int j,
                       real nu, real &s) {
  s = X[static_cast<size_t>(i)].dot(X[static_cast<size_t>(j)]) < real(0) ? real(-1) : real(1);
  return real(0.5) * nu *
         (x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)])
             .dot(X[static_cast<size_t>(i)] + s * X[static_cast<size_t>(j)]);
}

/// Builds the system for unit tangents `X` (zero entries contribute no
/// phase advance) at frequency `nu` (radians per unit length; the stripe
/// period is `2 pi / nu`). Triangles only.
inline system build(const asawa::shell::shell &M, const std::vector<vec3> &x,
                    const std::vector<vec3> &X, real nu) {
  const int nV = static_cast<int>(M.vert_count());
  system S;
  S.nu = nu;
  S.X = X;
  S.X.resize(static_cast<size_t>(nV), vec3::Zero());
  twisted_laplacian::assemble(
      M, x,
      [&](int i, int j) {
        twisted_laplacian::twist tw;
        tw.angle = edge_omega(x, S.X, i, j, nu, tw.sign);
        return tw;
      },
      S.A, S.B);
  real a_mean = 0, b_mean = 0;
  for (int i = 0; i < 2 * nV; ++i) {
    a_mean += S.A.coeff(i, i);
    b_mean += S.B[i];
  }
  S.shift = b_mean > real(0) ? real(1e-8) * a_mean / b_mean : real(0);
  for (int i = 0; i < 2 * nV; ++i)
    S.A.coeffRef(i, i) += S.shift * S.B[i];
  return S;
}

/// `psi^T (A - shift B) psi / psi^T B psi`: the stripe energy per unit mass.
inline real rayleigh(const system &S, const vecX &psi) {
  const real den = psi.dot(S.B.cwiseProduct(psi));
  return den > real(0) ? psi.dot(S.A * psi) / den - S.shift : real(0);
}

/// Inverse power iteration from `psi` (interleaved `(Re, Im)`; a zero or
/// wrongly sized start is replaced by a seeded random one), `iterations`
/// backsolves, `B`-normalized. Returns the Rayleigh quotient, NaN when the
/// factorization fails. Patterns of neighboring windings have nearly equal
/// energies, so a warm start keeps the pattern it starts from (it drifts
/// toward the minimum by the energy ratio per backsolve, under a percent on
/// a torus) while a cold start finds the minimum.
inline real solve(const system &S, vecX &psi, int iterations, unsigned seed = 1) {
  const int n = static_cast<int>(S.A.rows());
  if (psi.size() != n || !(psi.squaredNorm() > real(0)) || !psi.allFinite()) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> pick(real(-1), real(1));
    psi.resize(n);
    for (int i = 0; i < n; ++i)
      psi[i] = pick(rng);
  }
  matS A = S.A;
  m_solver solver(A);
  if (!solver.success())
    return std::numeric_limits<real>::quiet_NaN();
  for (int it = 0; it < iterations; ++it) {
    vecX b = S.B.cwiseProduct(psi);
    psi = solver.solve(b);
    const real nrm = std::sqrt(std::max(psi.dot(S.B.cwiseProduct(psi)), real(0)));
    if (!(nrm > real(0)) || !std::isfinite(nrm))
      return std::numeric_limits<real>::quiet_NaN();
    psi /= nrm;
  }
  return rayleigh(S, psi);
}

/// Phase change along edge `i -> j`: the target step `omega_ij` corrected by
/// the angle from the transported `psi_i` to `psi_j` (the paper's spinning
/// form), in `psi_j`'s sheet when `s_ij` is -1.
inline real edge_phase(const system &S, const std::vector<vec3> &x, const vecX &psi, int i, int j) {
  real s = 1;
  const real om = edge_omega(x, S.X, i, j, S.nu, s);
  const std::complex<real> zi(psi[2 * i], s * psi[2 * i + 1]);
  const std::complex<real> zj(psi[2 * j], psi[2 * j + 1]);
  return om - std::arg(std::polar(real(1), om) * zi / zj);
}

/// Per vertex the unit direction of the pattern's level lines, `N x grad
/// phase` summed over incident faces (area weighted), tangent to the vertex
/// normal. A face is skipped when it holds a zero of `psi` (its phase changes
/// around its boundary sum to a nonzero turn) or the field flips across one
/// of its edges; `skipped` counts those. A vertex with no usable face is zero.
inline std::vector<vec3> level_directions(const asawa::shell::shell &M,
                                          const std::vector<vec3> &x, const system &S,
                                          const vecX &psi, int *skipped = nullptr) {
  const int nV = static_cast<int>(M.vert_count());
  std::vector<vec3> sum(static_cast<size_t>(nV), vec3::Zero());
  int skip = 0;
  for (asawa::shell::FaceId f : M.get_face_range()) {
    if (M.fsize(f) != 3)
      continue;
    const auto tri = M.get_tri(f);
    const int v[3] = {static_cast<int>(tri[0]), static_cast<int>(tri[1]),
                      static_cast<int>(tri[2])};
    bool flips = false;
    for (int c = 0; c < 3; ++c)
      flips |= S.X[v[c]].dot(S.X[v[(c + 1) % 3]]) < real(0);
    const real dij = edge_phase(S, x, psi, v[0], v[1]);
    const real djk = edge_phase(S, x, psi, v[1], v[2]);
    const real dki = edge_phase(S, x, psi, v[2], v[0]);
    if (flips || std::abs(dij + djk + dki) > real(M_PI)) {
      ++skip;
      continue;
    }
    const vec3 cr = (x[v[1]] - x[v[0]]).cross(x[v[2]] - x[v[0]]);
    const real A2 = cr.norm();
    if (!(A2 > real(1e-20)))
      continue;
    const vec3 n = cr / A2;
    // Linear phase (0, dij, -dki) at (v0, v1, v2); grad lambda_c = n x e_c / A2
    // with e_c the edge opposite corner c, counterclockwise.
    const vec3 g = (dij * n.cross(x[v[0]] - x[v[2]]) - dki * n.cross(x[v[1]] - x[v[0]])) / A2;
    const vec3 along = real(0.5) * A2 * n.cross(g);
    for (int c = 0; c < 3; ++c)
      sum[static_cast<size_t>(v[c])] += along;
  }
  if (skipped)
    *skipped = skip;
  for (int i = 0; i < nV; ++i) {
    vec3 &d = sum[static_cast<size_t>(i)];
    vec3 nrm = asawa::shell::vert_normal(M, asawa::shell::vert_id(i), x);
    if (nrm.squaredNorm() > real(0))
      d -= nrm.normalized() * d.dot(nrm.normalized());
    d = d.squaredNorm() > real(1e-30) ? vec3(d.normalized()) : vec3::Zero();
  }
  return sum;
}

} // namespace stripes
} // namespace kusama
} // namespace gaudi

#endif
