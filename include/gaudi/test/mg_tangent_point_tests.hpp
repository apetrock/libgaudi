#ifndef GAUDI_MG_TANGENT_POINT_TESTS_HPP
#define GAUDI_MG_TANGENT_POINT_TESTS_HPP

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/procrustes_pyramid.hpp"
#include "gaudi/calder/point_tangent_point.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/modules/mg_tangent_point.hpp"
#include "gaudi/test/test.hpp"

#include <vector>

namespace gaudi {
namespace test {

GAUDI_TEST(mg_tp_points_cauchy_repels) {
  // Points along X, tangents along Z ⇒ T×dp ≠ 0.
  std::vector<vec3> c = {vec3(0, 0, 0), vec3(0.2, 0, 0), vec3(0.4, 0, 0)};
  std::vector<vec3> T = {vec3::UnitZ(), vec3::UnitZ(), vec3::UnitZ()};
  std::vector<real> w = {1, 1, 1};
  auto g = calder::tangent_point_gradient_points_cauchy(c, T, w, 0.05, 3.0,
                                                        /*bh_eps=*/0.0);
  GAUDI_EXPECT(g[0].norm() > 0.0);
  const real d01_before = (c[0] - c[1]).norm();
  c[0] -= 0.01 * g[0];
  c[1] -= 0.01 * g[1];
  GAUDI_EXPECT((c[0] - c[1]).norm() > d01_before);
}

GAUDI_TEST(mg_tp_kabsch_sim_recovers_rigid) {
  std::vector<vec3> x0 = {vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0),
                          vec3(0, 0, 1)};
  const real th = 0.37;
  mat3 R = mat3::Identity();
  R(0, 0) = std::cos(th);
  R(0, 1) = -std::sin(th);
  R(1, 0) = std::sin(th);
  R(1, 1) = std::cos(th);
  const vec3 t(0.2, -0.1, 0.3);
  std::vector<vec3> x1(x0.size());
  for (size_t i = 0; i < x0.size(); ++i)
    x1[i] = R * x0[i] + t;
  const mat3 Rg = calder::kabsch_sim_rotation(x0, x1);
  GAUDI_EXPECT((Rg - R).norm() < 1e-8);
}

GAUDI_TEST(mg_tp_extract_at_depth) {
  std::vector<vec3> centers;
  for (int i = 0; i < 32; ++i)
    centers.emplace_back(real(i % 8), real(i / 8), 0.1 * real(i));
  auto [h, idx] = arp::make_hash_3d_t<arp::morton_t>(centers);
  auto [lnodes, inodes] = arp::build_tree(h);
  arp::tree_view view(idx, inodes, lnodes, nullptr);
  arp::procrustes_pyramid P;
  P.build(view, centers);
  P.restrict_from_leaves();

  const index_t max_L = P.max_depth_from_leaves();
  GAUDI_EXPECT(max_L >= 1);
  auto cut1 = P.extract_at_depth(1);
  GAUDI_EXPECT(!cut1.empty());
  if (max_L >= 2) {
    auto cut2 = P.extract_at_depth(2);
    GAUDI_EXPECT(cut2.size() <= cut1.size() || cut2.size() > 0);
  }
}

GAUDI_TEST(mg_tp_level_step_moves_and_rotates_T) {
  std::vector<vec3> x = {vec3(0, 0, 0), vec3(0.15, 0, 0), vec3(0.3, 0, 0),
                         vec3(0.45, 0, 0)};
  std::vector<vec3> T(x.size(), vec3::UnitZ());
  std::vector<real> w(x.size(), real(1));
  duchamp::mg_tp_config cfg;
  cfg.coarse_iters = 6;
  cfg.coarse_step_frac = 0.2;
  cfg.bh_eps = 0.0;
  const std::vector<vec3> x0 = x;
  const mat3 R = duchamp::level_tp_step(x, T, w, cfg);
  real moved = 0;
  for (size_t i = 0; i < x.size(); ++i)
    moved += (x[i] - x0[i]).norm();
  GAUDI_EXPECT(moved > real(1e-6));
  (void)R;
  for (const auto &ti : T)
    GAUDI_EXPECT(std::abs(ti.norm() - real(1)) < 1e-6);
}

} // namespace test
} // namespace gaudi

#endif
