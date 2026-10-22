#ifndef GAUDI_GAUSS_NEWTON_TESTS_HPP
#define GAUDI_GAUSS_NEWTON_TESTS_HPP

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/kusama/principal_curvature_solver.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <random>
#include <vector>

namespace gaudi {
namespace test {
namespace {

/// GN pin: `r_i = w (x_i - t_i)`, the residual form of `pinned`.
class gn_pin : public hepworth::block::block_constraint {
public:
  DEFINE_CREATE_FUNC(gn_pin)
  gn_pin(std::vector<index_t> ids, std::vector<vec3> t, real w,
         std::vector<hepworth::sim_block::ptr> blocks)
      : block_constraint(ids, w, std::move(blocks)), _t(std::move(t)) {}
  bool gauss_newton() const override { return true; }
  void fill_J(const vecX &, index_t &id0, std::vector<trip> &trips) override {
    for (size_t k = 0; k < _ids.size(); ++k)
      for (int c = 0; c < 3; ++c)
        trips.push_back(trip(id0 + 3 * k + c,
                             _blocks[0]->get_offset_idx(_ids[k]) + c, _w));
    id0 += 3 * static_cast<index_t>(_ids.size());
  }
  void add_gradient(const vecX &q, vecX &g) override {
    for (size_t k = 0; k < _ids.size(); ++k) {
      const vec3 d = _blocks[0]->get_vec3(_ids[k], q) - _t[k];
      for (int c = 0; c < 3; ++c)
        g[_blocks[0]->get_offset_idx(_ids[k]) + c] += _w * _w * d[c];
    }
  }
  real residual_sq(const vecX &q) override {
    real e = 0;
    for (size_t k = 0; k < _ids.size(); ++k)
      e += _w * _w * (_blocks[0]->get_vec3(_ids[k], q) - _t[k]).squaredNorm();
    return e;
  }
  std::vector<vec3> _t;
};

/// GN sphere: `r_i = w (|x_i| - R)`, nonlinear.
class gn_sphere : public hepworth::block::block_constraint {
public:
  DEFINE_CREATE_FUNC(gn_sphere)
  gn_sphere(std::vector<index_t> ids, real R, real w,
            std::vector<hepworth::sim_block::ptr> blocks)
      : block_constraint(ids, w, std::move(blocks)), _R(R) {}
  bool gauss_newton() const override { return true; }
  void fill_J(const vecX &q, index_t &id0, std::vector<trip> &trips) override {
    for (size_t k = 0; k < _ids.size(); ++k) {
      const vec3 n = _blocks[0]->get_vec3(_ids[k], q).normalized();
      for (int c = 0; c < 3; ++c)
        trips.push_back(trip(id0 + k, _blocks[0]->get_offset_idx(_ids[k]) + c,
                             _w * n[c]));
    }
    id0 += static_cast<index_t>(_ids.size());
  }
  void add_gradient(const vecX &q, vecX &g) override {
    for (size_t k = 0; k < _ids.size(); ++k) {
      const vec3 x = _blocks[0]->get_vec3(_ids[k], q);
      const real r = _w * (x.norm() - _R);
      const vec3 n = x.normalized();
      for (int c = 0; c < 3; ++c)
        g[_blocks[0]->get_offset_idx(_ids[k]) + c] += _w * r * n[c];
    }
  }
  real residual_sq(const vecX &q) override {
    real e = 0;
    for (size_t k = 0; k < _ids.size(); ++k) {
      const real r = _w * (_blocks[0]->get_vec3(_ids[k], q).norm() - _R);
      e += r * r;
    }
    return e;
  }
  real _R;
};

inline std::vector<index_t> gn_all(size_t n) {
  std::vector<index_t> ids(n);
  for (size_t i = 0; i < n; ++i)
    ids[i] = static_cast<index_t>(i);
  return ids;
}

GAUDI_TEST(gn_linear_matches_projection) {
  std::mt19937 rng(11);
  std::uniform_real_distribution<real> u(-1, 1);
  const int n = 20;
  std::vector<vec3> x0(n), t(n);
  for (int i = 0; i < n; ++i) {
    x0[i] = vec3(u(rng), u(rng), u(rng));
    t[i] = vec3(u(rng), u(rng), u(rng));
  }
  auto run = [&](bool gn) {
    std::vector<vec3> x = x0, v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
    auto pos = hepworth::vec3_block::create(m, x, v, f);
    std::vector<hepworth::sim_block::ptr> blocks{pos};
    std::vector<hepworth::projection_constraint::ptr> cs;
    for (int i = 0; i < n; ++i) {
      if (gn)
        cs.push_back(gn_pin::create(std::vector<index_t>{i},
                                    std::vector<vec3>{t[i]}, real(3), blocks));
      else
        cs.push_back(hepworth::block::pinned::create(
            std::vector<index_t>{i}, t[i], real(3), blocks));
    }
    hepworth::block::projection_solver solver;
    solver.set_constraints(cs);
    solver.step(blocks, real(0.1), real(0.0), 1);
    return x;
  };
  const std::vector<vec3> xp = run(false);
  const std::vector<vec3> xg = run(true);
  real worst = 0;
  for (int i = 0; i < n; ++i)
    worst = std::max(worst, (xp[i] - xg[i]).norm());
  console_logger::info << "gn vs projection pin: max diff " << worst
                       << std::endl;
  GAUDI_EXPECT(worst < 1e-10);
}

GAUDI_TEST(gn_lands_on_sphere) {
  std::mt19937 rng(13);
  std::uniform_real_distribution<real> u(-1, 1);
  const int n = 60;
  std::vector<vec3> x0(n);
  for (int i = 0; i < n; ++i) {
    vec3 p(u(rng), u(rng), u(rng));
    x0[i] = p.normalized() * (0.5 + 0.5 * std::abs(u(rng)));
  }
  const real R = 1.3;
  for (const int every : {1, 4}) {
    std::vector<vec3> x = x0, v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
    auto pos = hepworth::vec3_block::create(m, x, v, f);
    std::vector<hepworth::sim_block::ptr> blocks{pos};
    std::vector<hepworth::projection_constraint::ptr> cs{
        gn_sphere::create(gn_all(n), R, real(100), blocks)};
    hepworth::block::projection_solver solver;
    solver.gn_refresh_every = every;
    solver.set_constraints(cs);
    int refreshes = 0;
    for (int s = 0; s < 4; ++s) {
      solver.step(blocks, real(0.1), real(0.5), 12);
      refreshes += solver.last_gn_refreshes;
    }
    real worst = 0;
    for (const vec3 &p : x)
      worst = std::max(worst, std::abs(p.norm() - R));
    console_logger::info << "gn sphere refresh_every " << every
                         << ": max | |x| - R | " << worst << ", refreshes "
                         << refreshes << std::endl;
    GAUDI_EXPECT(worst < 1e-3);
  }
}

} // namespace
} // namespace test
} // namespace gaudi

#endif
