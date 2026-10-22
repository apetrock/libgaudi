#ifndef GAUDI_HEPWORTH_BLOCK_VARIATIONAL_AUDIT_HPP
#define GAUDI_HEPWORTH_BLOCK_VARIATIONAL_AUDIT_HPP

/// Which constraints minimize an energy. A projection onto a set fixed during
/// the solve, `p(q) = argmin_{p in C} |A q - p|^2`, has energy
/// `E(q) = |A q - p(q)|^2` with gradient `2 A^T (A q - p(q))` (the envelope
/// theorem: `p`'s own dependence on `q` drops out at the minimizer). A
/// projection whose target moves with `q` (live charts, live goals, a polar
/// fit against the wrong spokes) breaks that identity, acts as a follower
/// force, and can add energy. The audit compares a central difference of `E`
/// along random directions with the envelope gradient, per constraint; for
/// Gauss–Newton rows it checks `residual_sq` against `2 add_gradient`.

#include "gaudi/hepworth/projection_constraint.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

struct variational_audit_row {
  std::string name;
  int count = 0;
  /// Constraints whose relative mismatch exceeds `flag_tol`.
  int flagged = 0;
  real median = 0;
  real max = 0;
  bool gauss_newton = false;
  /// Per-instance mismatch in constraint order, kept for small groups.
  std::vector<real> instances;
};

namespace detail {

/// `E(q)` and the envelope gradient's directional derivative along `d`.
inline void projection_energy(projection_constraint &c, const vecX &q, const vecX &d,
                              real &E, real *dE) {
  index_t id0 = 0;
  std::vector<trip> trips;
  c.fill_A(id0, trips);
  const index_t n = id0;
  if (n <= 0) {
    E = 0;
    if (dE)
      *dE = 0;
    return;
  }
  // Row-major: one constraint touches few rows of a system with many columns.
  Eigen::SparseMatrix<real, Eigen::RowMajor> A(n, q.size());
  A.setFromTriplets(trips.begin(), trips.end());
  vecX p = vecX::Zero(n);
  c.project(q, p);
  const vecX r = A * q - p;
  E = r.squaredNorm();
  if (dE)
    *dE = real(2) * (A * d).dot(r);
}

} // namespace detail

/// Relative mismatch `|FD - envelope| / max(|FD|, |envelope|, floor)` per
/// constraint, over `dirs` random directions of size `eps`, summarized by
/// constraint name. `floor` keeps constraints at rest (zero energy, zero
/// gradient) from reading as mismatched.
inline std::vector<variational_audit_row>
audit_variational(const std::vector<projection_constraint::ptr> &constraints, const vecX &q,
                  int dirs = 3, real eps = real(1e-6), real flag_tol = real(1e-3),
                  real floor = real(1e-12), unsigned seed = 1) {
  std::mt19937 rng(seed);
  std::normal_distribution<real> g(0, 1);
  std::map<std::string, std::vector<real>> by_name;
  std::map<std::string, bool> is_gn;
  std::map<std::string, std::vector<real>> order;
  for (const auto &cp : constraints) {
    if (!cp)
      continue;
    projection_constraint &c = *cp;
    const std::string name = c.name();
    // fill_A records the constraint's row offset; the solver's must survive.
    const index_t id0_keep = c._id0;
    real worst = 0;
    for (int t = 0; t < dirs; ++t) {
      vecX d(q.size());
      for (index_t i = 0; i < q.size(); ++i)
        d[i] = g(rng);
      d /= std::max(d.norm(), real(1e-30));
      real fd = 0, an = 0;
      if (c.gauss_newton()) {
        vecX grad = vecX::Zero(q.size());
        c.add_gradient(q, grad);
        an = real(2) * grad.dot(d);
        fd = (c.residual_sq(q + eps * d) - c.residual_sq(q - eps * d)) / (real(2) * eps);
      } else {
        real E0 = 0, Ep = 0, Em = 0;
        detail::projection_energy(c, q, d, E0, &an);
        detail::projection_energy(c, q + eps * d, d, Ep, nullptr);
        detail::projection_energy(c, q - eps * d, d, Em, nullptr);
        fd = (Ep - Em) / (real(2) * eps);
      }
      const real scale = std::max({std::abs(fd), std::abs(an), floor});
      worst = std::max(worst, std::abs(fd - an) / scale);
    }
    c._id0 = id0_keep;
    by_name[name].push_back(worst);
    order[name].push_back(worst);
    is_gn[name] = c.gauss_newton();
  }
  std::vector<variational_audit_row> out;
  for (auto &[name, v] : by_name) {
    std::sort(v.begin(), v.end());
    variational_audit_row row;
    row.name = name;
    row.count = static_cast<int>(v.size());
    row.median = v[v.size() / 2];
    row.max = v.back();
    row.flagged = static_cast<int>(std::count_if(v.begin(), v.end(), [&](real x) { return x > flag_tol; }));
    row.gauss_newton = is_gn[name];
    if (order[name].size() <= 4)
      row.instances = order[name];
    out.push_back(row);
  }
  std::sort(out.begin(), out.end(),
            [](const auto &a, const auto &b) { return a.median > b.median; });
  return out;
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
