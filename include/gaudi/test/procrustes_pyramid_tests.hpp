#ifndef GAUDI_PROCRUSTES_PYRAMID_TESTS_HPP
#define GAUDI_PROCRUSTES_PYRAMID_TESTS_HPP

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/procrustes_pyramid.hpp"
#include "gaudi/arp/tree_view.hpp"
#include "gaudi/common.h"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

namespace {

struct PointTreeFixture {
  std::vector<vec3> centers;
  std::vector<index_t> indices;
  std::vector<arp::radix_tree_node> internal_nodes;
  std::vector<arp::radix_tree_node> leaf_nodes;

  explicit PointTreeFixture(std::vector<vec3> pts) : centers(std::move(pts)) {
    auto [h, idx] = arp::make_hash_3d_t<arp::morton_t>(centers);
    indices = std::move(idx);
    auto [lnodes, inodes] = arp::build_tree(h);
    leaf_nodes = std::move(lnodes);
    internal_nodes = std::move(inodes);
  }

  arp::tree_view view() const {
    return arp::tree_view(indices, internal_nodes, leaf_nodes, nullptr);
  }
};

inline PointTreeFixture make_grid_fixture() {
  std::vector<vec3> centers;
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y)
      for (int z = 0; z < 4; ++z)
        centers.emplace_back(real(x), real(y), real(z) + 0.1 * real(x * y));
  return PointTreeFixture(std::move(centers));
}

inline PointTreeFixture make_two_cluster_fixture() {
  std::vector<vec3> centers;
  for (int i = 0; i < 8; ++i)
    centers.emplace_back(real(i) * 0.1, real(i % 3) * 0.1, real(i % 2) * 0.1);
  for (int i = 0; i < 8; ++i)
    centers.emplace_back(10.0 + real(i) * 0.1, real(i % 3) * 0.1,
                         real(i % 2) * 0.1);
  return PointTreeFixture(std::move(centers));
}

inline real max_abs_diff(const std::vector<vec3> &a,
                         const std::vector<vec3> &b) {
  GAUDI_ASSERT(a.size() == b.size());
  real m = 0;
  for (size_t i = 0; i < a.size(); ++i)
    m = std::max(m, (a[i] - b[i]).cwiseAbs().maxCoeff());
  return m;
}

inline mat3 rot_z(real theta) {
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

GAUDI_TEST(procrustes_avg_commutes_with_R) {
  auto fix = make_grid_fixture();
  vec3 avg = vec3::Zero();
  for (const auto &p : fix.centers)
    avg += p;
  avg /= real(fix.centers.size());

  const mat3 R = rot_z(0.37);
  const vec3 t(1.2, -0.4, 0.8);
  vec3 avg_t = vec3::Zero();
  for (const auto &p : fix.centers)
    avg_t += R * p + t;
  avg_t /= real(fix.centers.size());

  const vec3 expect = R * avg + t;
  GAUDI_EXPECT((avg_t - expect).norm() < 1e-12);
}

GAUDI_TEST(procrustes_pyramid_roundtrip_identity) {
  auto fix = make_grid_fixture();
  arp::procrustes_pyramid P;
  P.build(fix.view(), fix.centers);
  GAUDI_ASSERT(P.root() != arp::UNULL);

  const auto before = P.leaf_positions_original();
  P.prolong_identity_from_root();
  const auto after = P.leaf_positions_original();
  GAUDI_EXPECT(max_abs_diff(before, after) < 1e-12);
  GAUDI_EXPECT(max_abs_diff(after, fix.centers) < 1e-12);
}

GAUDI_TEST(procrustes_pyramid_rigid_mutate_prolong_recovers) {
  auto fix = make_grid_fixture();
  arp::procrustes_pyramid P;
  P.build(fix.view(), fix.centers);
  const index_t root = P.root();
  GAUDI_ASSERT(root != arp::UNULL);

  const mat3 R = rot_z(0.51);
  const vec3 t(-0.3, 1.1, 0.2);
  const vec3 c0 = P.avg_rest()[static_cast<size_t>(root)];
  const std::vector<vec3> delta0 = P.leaf_delta();

  P.set_node_pose(root, R * c0 + t, R);
  P.prolong_from(root);

  const auto got = P.leaf_positions_original();
  for (size_t i = 0; i < fix.centers.size(); ++i) {
    const vec3 expect = R * fix.centers[i] + t;
    GAUDI_EXPECT((got[i] - expect).norm() < 1e-10);
  }

  // Re-encode: δ magnitudes in parent frame should match rest (rigid).
  P.encode_from_leaves();
  real max_dmag = 0;
  for (size_t i = 0; i < delta0.size(); ++i)
    max_dmag =
        std::max(max_dmag, std::abs(delta0[i].norm() - P.leaf_delta()[i].norm()));
  GAUDI_EXPECT(max_dmag < 1e-10);
}

GAUDI_TEST(procrustes_pyramid_midlevel_set_reencode) {
  auto fix = make_grid_fixture();
  arp::procrustes_pyramid P;
  P.build(fix.view(), fix.centers);
  const index_t root = P.root();
  GAUDI_ASSERT(root != arp::UNULL);
  GAUDI_ASSERT(P.node_count() > 2);

  // Pick a non-root internal node.
  index_t mid = arp::UNULL;
  for (size_t i = 0; i < P.node_count(); ++i) {
    if (fix.internal_nodes[i].parent == root) {
      mid = static_cast<index_t>(i);
      break;
    }
  }
  GAUDI_ASSERT(mid != arp::UNULL);

  const mat3 R0 = rot_z(0.2);
  const vec3 t0(0.5, 0.0, 0.0);
  P.set_node_pose(root, R0 * P.avg_rest()[static_cast<size_t>(root)] + t0, R0);
  P.prolong_from(root);

  const mat3 R1 = rot_z(-0.35);
  const vec3 c_mid = P.avg_live()[static_cast<size_t>(mid)];
  P.set_node_pose(mid, R1 * (c_mid - t0) + t0 + vec3(0.0, 0.25, 0.0), R0 * R1);
  P.prolong_from(mid);

  const auto after_edit = P.leaf_positions_original();
  P.encode_from_leaves();
  P.set_node_pose(root, P.avg_live()[static_cast<size_t>(root)],
                  P.R_live()[static_cast<size_t>(root)]);
  P.prolong_from(root);
  const auto recovered = P.leaf_positions_original();
  GAUDI_EXPECT(max_abs_diff(after_edit, recovered) < 1e-9);
}

GAUDI_TEST(procrustes_pyramid_two_cluster_smoke) {
  auto fix = make_two_cluster_fixture();
  arp::procrustes_pyramid P;
  P.build(fix.view(), fix.centers);
  GAUDI_ASSERT(P.root() != arp::UNULL);

  P.prolong_identity_from_root();
  GAUDI_EXPECT(max_abs_diff(P.leaf_positions_original(), fix.centers) < 1e-12);

  const mat3 R = rot_z(0.15);
  const vec3 t(0.0, 0.0, 2.0);
  const index_t root = P.root();
  P.set_node_pose(root, R * P.avg_rest()[static_cast<size_t>(root)] + t, R);
  P.prolong_from(root);
  const auto got = P.leaf_positions_original();
  for (size_t i = 0; i < fix.centers.size(); ++i)
    GAUDI_EXPECT((got[i] - (R * fix.centers[i] + t)).norm() < 1e-10);
}

GAUDI_TEST(procrustes_restrict_kabsch_recovers_rigid) {
  auto fix = make_grid_fixture();
  arp::procrustes_pyramid P;
  P.build(fix.view(), fix.centers);

  const mat3 R = rot_z(0.42);
  const vec3 t(0.7, -0.2, 0.3);
  std::vector<vec3> moved(fix.centers.size());
  for (size_t i = 0; i < fix.centers.size(); ++i)
    moved[i] = R * fix.centers[i] + t;

  P.set_leaf_positions_original(moved);
  P.restrict_from_leaves();

  const index_t root = P.root();
  GAUDI_ASSERT(root != arp::UNULL);
  // Root Kabsch should match R (up to sign ambiguity on near-isotropic clouds —
  // grid is not isotropic; expect close Frobenius).
  const mat3 &Rr = P.R_live()[static_cast<size_t>(root)];
  GAUDI_EXPECT((Rr - R).norm() < 1e-6);
  GAUDI_EXPECT((P.avg_live()[static_cast<size_t>(root)] - (R * P.avg_rest()[static_cast<size_t>(root)] + t))
                   .norm() < 1e-6);

  P.prolong_from(root);
  GAUDI_EXPECT(max_abs_diff(P.leaf_positions_original(), moved) < 1e-8);
}

GAUDI_TEST(procrustes_extract_cut_covers_leaves) {
  auto fix = make_grid_fixture();
  arp::procrustes_pyramid P;
  P.build(fix.view(), fix.centers);
  const auto counts = P.subtree_leaf_counts();
  const real total = counts[static_cast<size_t>(P.root())];

  // Coarse: allow large clusters
  auto cut_coarse = P.extract_cut(total * 0.5);
  GAUDI_EXPECT(!cut_coarse.empty());
  GAUDI_EXPECT(cut_coarse.size() < P.node_count());

  // Finer cut
  auto cut_fine = P.extract_cut(4.0);
  GAUDI_EXPECT(cut_fine.size() >= cut_coarse.size());

  real sum = 0;
  for (index_t n : cut_fine)
    sum += counts[static_cast<size_t>(n)];
  GAUDI_EXPECT(std::abs(sum - total) < 1e-9);
}

GAUDI_TEST(procrustes_coarse_translate_prolong_additive) {
  auto fix = make_grid_fixture();
  arp::procrustes_pyramid P;
  P.build(fix.view(), fix.centers);
  P.restrict_from_leaves();

  auto cut = P.extract_cut(16.0);
  GAUDI_ASSERT(!cut.empty());
  const vec3 dc(0.1, -0.05, 0.2);
  const auto before = P.leaf_positions_original();

  for (index_t n : cut) {
    const vec3 c0 = P.avg_live()[static_cast<size_t>(n)];
    const mat3 R0 = P.R_live()[static_cast<size_t>(n)];
    P.set_node_pose(n, c0 + dc, R0);
    P.prolong_from(n);
  }
  const auto after = P.leaf_positions_original();
  for (size_t i = 0; i < before.size(); ++i)
    GAUDI_EXPECT((after[i] - (before[i] + dc)).norm() < 1e-9);
}

} // namespace test
} // namespace gaudi

#endif
