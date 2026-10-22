#ifndef GAUDI_RESIDUAL_FOREST_TESTS_HPP
#define GAUDI_RESIDUAL_FOREST_TESTS_HPP

#include "gaudi/asawa/rod/residual_forest.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

namespace {

inline std::vector<vec3> make_helix(int n = 64) {
  std::vector<vec3> pts;
  pts.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    const real t = real(2) * real(3.14159265358979323846) * real(i) / real(n);
    pts.emplace_back(std::cos(t), std::sin(t), real(0.15) * t);
  }
  return pts;
}

inline real forest_max_abs_diff(const std::vector<vec3> &a,
                                const std::vector<vec3> &b) {
  GAUDI_ASSERT(a.size() == b.size());
  real m = 0;
  for (size_t i = 0; i < a.size(); ++i)
    m = std::max(m, (a[i] - b[i]).cwiseAbs().maxCoeff());
  return m;
}

inline real residual_l2(const std::vector<vec3> &d) {
  real s = 0;
  for (const auto &v : d)
    s += v.squaredNorm();
  return std::sqrt(s);
}

inline real forest_field_dot(const std::vector<vec3> &a,
                             const std::vector<vec3> &b) {
  const size_t n = std::min(a.size(), b.size());
  real s = 0;
  for (size_t i = 0; i < n; ++i)
    s += a[i].dot(b[i]);
  return s;
}

inline mat3 forest_rot_z(real theta) {
  const real c = std::cos(theta);
  const real s = std::sin(theta);
  mat3 R = mat3::Identity();
  R(0, 0) = c;
  R(0, 1) = -s;
  R(1, 0) = s;
  R(1, 1) = c;
  return R;
}

} // namespace

GAUDI_TEST(rod_forest_identity_recover) {
  auto R = asawa::rod::rod::create(make_helix(64), true);
  const std::vector<vec3> x0 = R->x();
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F->valid());
  GAUDI_ASSERT(F->n_levels() >= 2);

  const int N = F->n_levels() - 1;
  F->set(N, F->pickle(N));
  F->prolong();
  GAUDI_ASSERT(forest_max_abs_diff(R->x(), x0) < real(1e-9));
}

GAUDI_TEST(rod_forest_rigid_coarse) {
  auto R = asawa::rod::rod::create(make_helix(64), true);
  const std::vector<vec3> x0 = R->x();
  auto F = asawa::rod::residual_forest::create(R);
  const int N = F->n_levels() - 1;

  std::vector<real> delta_norms;
  for (int ell = 0; ell < N; ++ell)
    delta_norms.push_back(residual_l2(F->residual(ell)));

  const mat3 Rg = forest_rot_z(real(0.35));
  const vec3 t(real(1.2), real(-0.4), real(0.7));
  std::vector<vec3> xN = F->pickle(N);
  for (auto &p : xN)
    p = Rg * p + t;
  F->set(N, xN);
  F->prolong();

  for (size_t i = 0; i < x0.size(); ++i) {
    const vec3 expect = Rg * x0[i] + t;
    GAUDI_ASSERT((R->x()[i] - expect).norm() < real(1e-5));
  }
  for (int ell = 0; ell < N; ++ell)
    GAUDI_ASSERT(std::abs(residual_l2(F->residual(ell)) -
                          delta_norms[static_cast<size_t>(ell)]) < real(1e-9));
}

GAUDI_TEST(rod_forest_midlevel_reencode) {
  auto R = asawa::rod::rod::create(make_helix(64), true);
  auto F = asawa::rod::residual_forest::create(R);
  const int N = F->n_levels() - 1;
  GAUDI_ASSERT(N >= 2);
  const int M = N - 1;

  {
    std::vector<vec3> xN = F->pickle(N);
    for (auto &p : xN)
      p += vec3(real(0.05), 0, 0);
    F->set(N, xN);
  }

  const real dM0 = residual_l2(F->residual(M));
  std::vector<vec3> xM = F->pickle(M);
  for (size_t i = 0; i < xM.size(); ++i)
    xM[i] += vec3(0, real(0.02) * real(i % 3), 0);
  F->set(M, xM);
  const real dM1 = residual_l2(F->residual(M));
  GAUDI_ASSERT(std::abs(dM1 - dM0) > real(1e-8));

  F->prolong();
  const std::vector<vec3> x0_after_mid = R->x();

  {
    const mat3 Rg = forest_rot_z(real(0.2));
    std::vector<vec3> xN = F->pickle(N);
    for (auto &p : xN)
      p = Rg * p;
    F->set(N, xN);
  }
  F->prolong();

  GAUDI_ASSERT(forest_max_abs_diff(R->x(), x0_after_mid) > real(1e-4));
  GAUDI_ASSERT(std::abs(residual_l2(F->residual(M)) - dM1) < real(1e-6));
}

GAUDI_TEST(rod_forest_invalidate) {
  auto R = asawa::rod::rod::create(make_helix(32), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F->valid());
  F->invalidate();
  GAUDI_ASSERT(!F->valid());
  bool threw = false;
  try {
    (void)F->pickle(0);
  } catch (...) {
    threw = true;
  }
  GAUDI_ASSERT(threw);
}

GAUDI_TEST(rod_forest_pickle_topology) {
  auto R = asawa::rod::rod::create(make_helix(64), true);
  auto F = asawa::rod::residual_forest::create(R);
  for (int ell = 0; ell < F->n_levels(); ++ell) {
    auto Rc = F->pickle_rod(ell);
    GAUDI_ASSERT(Rc->corner_count() == F->live(ell).size());
    index_t linked = 0;
    for (size_t i = 0; i < Rc->corner_count(); ++i) {
      if (Rc->next(asawa::rod::corner_id(static_cast<int>(i))) >=
          asawa::rod::corner_id(0))
        ++linked;
    }
    GAUDI_ASSERT(linked == static_cast<index_t>(Rc->corner_count()));
  }
}

GAUDI_TEST(rod_forest_k_depth) {
  auto R = asawa::rod::rod::create(make_helix(64), true);
  auto F5 = asawa::rod::residual_forest::create(R, /*k=*/5);
  GAUDI_ASSERT(F5->k() == 5);
  GAUDI_ASSERT(F5->min_coarse() == 10);
  const int n_coarse5 =
      static_cast<int>(F5->live(F5->n_levels() - 1).size());
  GAUDI_ASSERT(n_coarse5 > 10);

  auto F21 = asawa::rod::residual_forest::create(R, /*k=*/21);
  GAUDI_ASSERT(F21->k() == 21);
  GAUDI_ASSERT(F21->n_levels() <= F5->n_levels());
  GAUDI_ASSERT(static_cast<int>(F21->live(F21->n_levels() - 1).size()) >
               F21->min_coarse());
}

GAUDI_TEST(rod_forest_restrict_prolong) {
  auto R = asawa::rod::rod::create(make_helix(64), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F->n_levels() >= 2);

  std::vector<vec3> x0 = F->pickle(0);
  for (auto &p : x0)
    p += vec3(real(0.01), 0, 0);
  F->set(0, x0);
  F->restrict_up(0);
  // Coarse positions solve Jᵀ J x_c = Jᵀ x_f (reference ApplyPinv).
  const std::vector<vec3> Jxc = F->prolong_J(0, F->live(1));
  const std::vector<vec3> JTJxc = F->restrict_JT(0, Jxc);
  const std::vector<vec3> JTx = F->restrict_JT(0, x0);
  GAUDI_ASSERT(JTJxc.size() == JTx.size());
  for (size_t j = 0; j < JTJxc.size(); ++j)
    GAUDI_ASSERT((JTJxc[j] - JTx[j]).norm() < real(1e-8));
}

GAUDI_TEST(rod_forest_restrict_after_split) {
  auto R = asawa::rod::rod::create(make_helix(48), true);
  // Insert corners mid-strand (pool order ≠ curve order after pack).
  for (int s = 0; s < 6; ++s) {
    const int n = static_cast<int>(R->corner_count());
    for (int i = 0; i < n; i += 7) {
      asawa::rod::CornerId ci = asawa::rod::corner_id(i);
      if (R->next(ci) < asawa::rod::corner_id(0))
        continue;
      asawa::rod::CornerId c10 = R->next(ci);
      asawa::rod::CornerId c11 = R->next(c10);
      asawa::rod::CornerId c01 = R->prev(ci);
      asawa::rod::CornerId cnew = R->insert_edge();
      R->link(ci, cnew);
      R->link(cnew, c10);
      for (auto &d : R->get_data())
        d->rod_split(*R, c01, ci, c10, c11, cnew);
    }
    R->pack();
  }
  R->_update_frames();

  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F->n_levels() >= 2);
  F->restrict_up(0);

  const auto &keep = F->keep(0);
  const auto &x0 = F->live(0);
  const auto &x1 = F->live(1);
  GAUDI_ASSERT(x1.size() == keep.size());
  const std::vector<vec3> Jxc = F->prolong_J(0, x1);
  const std::vector<vec3> JTJxc = F->restrict_JT(0, Jxc);
  const std::vector<vec3> JTx = F->restrict_JT(0, x0);
  GAUDI_ASSERT(JTJxc.size() == JTx.size());
  for (size_t j = 0; j < JTJxc.size(); ++j)
    GAUDI_ASSERT((JTJxc[j] - JTx[j]).norm() < real(1e-8));

  // Coarse edges must be local (thinned knot), not chordal diameter.
  real mean_e = 0;
  real max_e = 0;
  const int nc = static_cast<int>(x1.size());
  for (int j = 0; j < nc; ++j) {
    const int jn = (j + 1) % nc;
    const real e = (x1[static_cast<size_t>(jn)] - x1[static_cast<size_t>(j)]).norm();
    mean_e += e;
    max_e = std::max(max_e, e);
  }
  mean_e /= real(nc);
  const real lavg = R->lavg();
  GAUDI_ASSERT(mean_e < real(6) * lavg);
  GAUDI_ASSERT(max_e < real(20) * lavg);
}

GAUDI_TEST(rod_forest_J_adjoint) {
  auto R = asawa::rod::rod::create(make_helix(48), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F->valid() && F->n_levels() >= 2);
  GAUDI_ASSERT(F->j_idx_slot(0) >= 0);

  const size_t nc = F->level_rod(1)->corner_count();
  const size_t nf = R->corner_count();
  std::vector<vec3> e(nc), r(nf);
  for (size_t i = 0; i < nc; ++i) {
    const real t = real(i);
    e[i] = vec3(std::sin(t), std::cos(real(0.7) * t), real(0.1) * t);
  }
  for (size_t i = 0; i < nf; ++i) {
    const real t = real(i);
    r[i] = vec3(std::cos(real(1.1) * t), std::sin(t), real(-0.05) * t);
  }
  const std::vector<vec3> Je = F->prolong_J(0, e);
  const std::vector<vec3> JTr = F->restrict_JT(0, r);
  const real lhs = forest_field_dot(Je, r);
  const real rhs = forest_field_dot(e, JTr);
  GAUDI_ASSERT(std::abs(lhs - rhs) < real(1e-9) * (1 + std::abs(lhs)));
}

} // namespace test
} // namespace gaudi

#endif // GAUDI_RESIDUAL_FOREST_TESTS_HPP
