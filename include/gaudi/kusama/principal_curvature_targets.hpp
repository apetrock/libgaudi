#ifndef GAUDI_KUSAMA_PRINCIPAL_CURVATURE_TARGETS_HPP
#define GAUDI_KUSAMA_PRINCIPAL_CURVATURE_TARGETS_HPP

/// Target principal curvatures for `principal_curvature_solver`.

#include "gaudi/common.h"

#include <cmath>
#include <limits>
#include <vector>

namespace gaudi {
namespace kusama {

/// Eigenvalues of the tangent tensor `(k11, k12, k22)`, larger first.
inline vec2 principal_curvature_sym2_eigs(const vec3 &c) {
  const real m = real(0.5) * (c[0] + c[2]);
  const real d = std::sqrt(real(0.25) * (c[0] - c[2]) * (c[0] - c[2]) +
                           c[1] * c[1]);
  return vec2(m + d, m - d);
}

/// Sign that puts a second operator's tensors in Eigensatz's convention
/// (convex positive): the sign of the mean-curvature correlation over the
/// vertices where both are defined.
inline real principal_curvature_sign_against(const std::vector<vec2> &k,
                                   const std::vector<vec3> &S,
                                   const std::vector<char> &ok) {
  real c = 0;
  for (size_t v = 0; v < k.size() && v < S.size() && v < ok.size(); ++v)
    if (ok[v] && k[v].allFinite() && S[v].allFinite())
      c += (k[v][0] + k[v][1]) * (S[v][0] + S[v][2]);
  return c < real(0) ? real(-1) : real(1);
}

/// `kappa_hat = kappa0 + (eig(H) - eig(S0))`: the change another operator's
/// shape block asks for, laid onto Eigensatz's own start-of-step curvature,
/// so the gap between the two discretizations is never a target. Vertices
/// without `ok` or a finite `kappa0` get NaN (no row).
inline std::vector<vec2>
principal_curvature_delta_targets(const std::vector<vec2> &kappa0,
                        const std::vector<vec3> &H, const std::vector<vec3> &S0,
                        const std::vector<char> &ok, real sign) {
  std::vector<vec2> kh(kappa0.size(),
                       vec2::Constant(std::numeric_limits<real>::quiet_NaN()));
  for (size_t v = 0; v < kappa0.size(); ++v) {
    if (v >= H.size() || v >= S0.size() || v >= ok.size() || !ok[v] ||
        !kappa0[v].allFinite() || !H[v].allFinite() || !S0[v].allFinite())
      continue;
    kh[v] = kappa0[v] + (principal_curvature_sym2_eigs(sign * H[v]) -
                         principal_curvature_sym2_eigs(sign * S0[v]));
  }
  return kh;
}

/// Analytic torus about +Z (outward normals): tube `1 / tube`, ring
/// `cos(phi) / (R + tube cos(phi))`, with `phi` read from each position.
inline std::vector<vec2> principal_curvature_torus_targets(const std::vector<vec3> &x,
                                                 real R, real tube) {
  std::vector<vec2> kh(x.size());
  for (size_t v = 0; v < x.size(); ++v) {
    const real phi = std::atan2(x[v][2], std::hypot(x[v][0], x[v][1]) - R);
    kh[v] = vec2(real(1) / tube,
                 std::cos(phi) / (R + tube * std::cos(phi)));
  }
  return kh;
}

} // namespace kusama
} // namespace gaudi

#endif
