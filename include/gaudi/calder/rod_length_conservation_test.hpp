#ifndef GAUDI_CALDER_ROD_LENGTH_CONSERVATION_TEST_HPP
#define GAUDI_CALDER_ROD_LENGTH_CONSERVATION_TEST_HPP

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/arp/datums.hpp"
#include "gaudi/calder/leaf_count.hpp"
#include "gaudi/calder/tree_code.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace gaudi {
namespace test {

namespace {

template <typename T>
T t2_scalar_datum(calder::fast_summation<arp::T2>::Node_Type node_type,
                  index_t j, const std::vector<calder::datum::ptr> &data) {
  const typename calder::datum_t<T>::ptr F =
      std::static_pointer_cast<typename calder::datum_t<T>>(data[0]);
  if (node_type == calder::fast_summation<arp::T2>::Node_Type::LEAF)
    return F->sorted_leaf_data()[j];
  return F->node_data()[j];
}

asawa::rod::rod::ptr make_circle_rod(int n_seg, real radius) {
  GAUDI_ASSERT(n_seg >= 8);
  GAUDI_ASSERT(radius > 1e-12);
  std::vector<vec3> verts;
  verts.reserve(static_cast<size_t>(n_seg));
  for (int i = 0; i < n_seg; ++i) {
    const real theta = 2.0 * M_PI * static_cast<real>(i) / static_cast<real>(n_seg);
    // Slight out-of-plane component so the T3-style AABB volume (tree_code) is
    // nonzero; a perfectly flat circle can have V=0 → sc=0 → only branch path.
    const real z = 0.0005 * std::sin(2.0 * theta);
    verts.push_back(vec3(radius * std::cos(theta), radius * std::sin(theta), z));
  }
  return asawa::rod::rod::create(verts, true);
}

real rod_total_edge_length(asawa::rod::rod &R) {
  real L = 0;
  for (asawa::rod::CornerId ci : R.get_vert_range()) {
    L += R.l0()[static_cast<int>(ci)];
  }
  return L;
}

real rod_recovered_scalar_length(asawa::rod::rod &R, const vec3 &pi, real eps) {
  const std::vector<vec3> &x = R.x();
  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<asawa::rod::CornerId> rverts = R.get_vert_range();
  std::vector<index_t> edge_ids(rverts.begin(), rverts.end());
  // Corner-id indexed lengths (same layout as integrate_over_rod / Darboux).
  const std::vector<real> &lc = R.l0();
  arp::T2::ptr tree = arp::T2::create(edge_verts, x, 12);
  calder::fast_summation<arp::T2> sum(*tree);
  sum.bind(calder::scalar_datum::create(edge_ids, lc));
  std::vector<vec3> pov = {pi};
  std::vector<real> u = sum.calc<real>(
      pov,
      [](const index_t &i, const index_t &j, const vec3 &pi,
         const std::vector<calder::datum::ptr> &data,
         calder::fast_summation<arp::T2>::Node_Type node_type,
         const arp::T2 &tree) -> real {
        (void)i;
        (void)pi;
        (void)tree;
        return t2_scalar_datum<real>(node_type, j, data);
      },
      [](const index_t &i, const index_t &j, const vec3 &pi,
         const std::vector<calder::datum::ptr> &data,
         calder::fast_summation<arp::T2>::Node_Type node_type,
         const arp::T2 &tree) -> real {
        (void)i;
        (void)pi;
        (void)tree;
        return t2_scalar_datum<real>(node_type, j, data);
      },
      eps, false);
  return u[0];
}

} // namespace

GAUDI_TEST(calder_rod_pyramid_length_matches_total_edge_length) {
  auto rod_ptr = make_circle_rod(64, 2.0);
  asawa::rod::rod &R = *rod_ptr;

  const real total_len = rod_total_edge_length(R);
  GAUDI_ASSERT(total_len > 1e-12);

  std::vector<index_t> ev = R.get_edge_vert_ids();
  GAUDI_ASSERT(ev.size() >= 2);
  const std::vector<vec3> &x = R.x();
  const vec3 on_vertex = x[ev[0]];
  const vec3 mid_edge = 0.5 * (x[ev[0]] + x[ev[1]]);
  const vec3 far(1000.0, 0.0, 0.0);

  const real eps_bh = 0.5;
  const real eps_tight = 1e-6;
  const real tol = 1e-4 * std::max(total_len, real(1.0));

  struct Row {
    const char *label;
    vec3 p;
    real eps;
  };
  const Row rows[] = {
      {"vertex (__x)", on_vertex, eps_bh},
      {"vertex (__x)", on_vertex, eps_tight},
      {"mid_edge (__x)", mid_edge, eps_bh},
      {"mid_edge (__x)", mid_edge, eps_tight},
      {"far off curve", far, eps_bh},
      {"far off curve", far, eps_tight},
  };

  std::cerr << "\n[calder_rod_pyramid_length]\n"
            << "  total_edge_length (sum over edges) = " << std::fixed
            << std::setprecision(8) << total_len << "\n"
            << "  edge_count = "
            << static_cast<int>(R.get_vert_range().size()) << "\n"
            << "  tol = " << std::scientific << tol << std::fixed << "\n\n"
            << "  " << std::setw(22) << "POV"
            << "  " << std::setw(10) << "eps"
            << "  " << std::setw(16) << "expected"
            << "  " << std::setw(16) << "recovered"
            << "  " << std::setw(14) << "|delta|"
            << "  " << std::setw(12) << "leaf_hits"
            << "\n";

  for (const Row &row : rows) {
    const real rec = rod_recovered_scalar_length(R, row.p, row.eps);
    const real delta = std::abs(rec - total_len);
    const real leaves =
        calder::rod_leaf_visit_count(R, row.p, row.eps);
    std::cerr << "  " << std::setw(22) << row.label << "  ";
    if (row.eps >= real(0.1))
      std::cerr << std::setw(10) << std::fixed << std::setprecision(2)
                << row.eps;
    else
      std::cerr << std::setw(10) << std::scientific << std::setprecision(1)
                << row.eps;
    std::cerr << std::fixed << std::setprecision(8) << "  " << std::setw(16)
              << total_len << "  " << std::setw(16) << rec << "  "
              << std::setw(14) << delta << "  " << std::setw(12)
              << static_cast<int>(leaves + real(0.5)) << "\n";
    GAUDI_EXPECT(delta < tol);
  }
  std::cerr << "  (leaf_hits = leaf callbacks from leaf->1, branch->0)\n\n";
}

GAUDI_TEST(calder_rod_leaf_visit_counts) {
  auto rod_ptr = make_circle_rod(64, 2.0);
  asawa::rod::rod &R = *rod_ptr;

  // T2::create uses __x vertex positions (see integrate_over_rod); queries must
  // use the same coordinates or opening-angle tests see the wrong geometry.
  const std::vector<vec3> &x = R.x();
  std::vector<index_t> ev = R.get_edge_vert_ids();
  const vec3 on_vertex = x[ev[0]];
  const vec3 mid_edge = 0.5 * (x[ev[0]] + x[ev[1]]);
  const vec3 far(1000.0, 0.0, 0.0);

  const real eps_bh = 0.5;
  const real eps_leaf = 1e-6;
  const real n_far_bh = calder::rod_leaf_visit_count(R, far, eps_bh);
  const real n_far_tight = calder::rod_leaf_visit_count(R, far, eps_leaf);
  const real n_vertex_bh =
      calder::rod_leaf_visit_count(R, on_vertex, eps_bh);
  const real n_vertex = calder::rod_leaf_visit_count(R, on_vertex, eps_leaf);
  const real n_mid_bh = calder::rod_leaf_visit_count(R, mid_edge, eps_bh);
  const real n_mid = calder::rod_leaf_visit_count(R, mid_edge, eps_leaf);

  std::cerr << "\n[calder_rod_leaf_visit_counts]\n"
            << "  (queries in __x space; tree built from R.x())\n"
            << "  " << std::setw(22) << "POV"
            << "  " << std::setw(12) << "eps=0.5"
            << "  " << std::setw(12) << "eps=1e-6"
            << "\n"
            << "  " << std::setw(22) << "far"
            << "  " << std::setw(12) << static_cast<int>(n_far_bh + real(0.5))
            << "  " << std::setw(12) << static_cast<int>(n_far_tight + real(0.5))
            << "\n"
            << "  " << std::setw(22) << "vertex (__x)"
            << "  " << std::setw(12) << static_cast<int>(n_vertex_bh + real(0.5))
            << "  " << std::setw(12) << static_cast<int>(n_vertex + real(0.5))
            << "\n"
            << "  " << std::setw(22) << "mid_edge (__x)"
            << "  " << std::setw(12) << static_cast<int>(n_mid_bh + real(0.5))
            << "  " << std::setw(12) << static_cast<int>(n_mid + real(0.5))
            << "\n\n";

  GAUDI_EXPECT(n_vertex > 0.5);
  GAUDI_EXPECT(n_mid > 0.5);
}

} // namespace test
} // namespace gaudi

#endif
