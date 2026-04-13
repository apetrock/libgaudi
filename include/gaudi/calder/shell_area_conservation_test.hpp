#ifndef GAUDI_CALDER_SHELL_AREA_CONSERVATION_TEST_HPP
#define GAUDI_CALDER_SHELL_AREA_CONSERVATION_TEST_HPP

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/arp/datums.hpp"
#include "gaudi/calder/leaf_count.hpp"
#include "gaudi/calder/tree_code.hpp"
#include "gaudi/test/bvh_tests.hpp"
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
T t3_scalar_datum(calder::fast_summation<arp::T3>::Node_Type node_type,
                  index_t j, const std::vector<calder::datum::ptr> &data) {
  const typename calder::datum_t<T>::ptr F =
      std::static_pointer_cast<typename calder::datum_t<T>>(data[0]);
  if (node_type == calder::fast_summation<arp::T3>::Node_Type::LEAF)
    return F->sorted_leaf_data()[j];
  return F->node_data()[j];
}

real shell_recovered_scalar_area(asawa::shell::shell &M,
                                  const std::vector<vec3> &x, const vec3 &pi,
                                  real eps) {
  std::vector<real> areas = asawa::shell::face_areas(M, x);
  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  auto face_ids = M.get_face_range();
  std::vector<index_t> face_ix(face_ids.begin(), face_ids.end());
  std::vector<real> Ac =
      asawa::shell::compress_to_range<real>(face_ix, areas);
  arp::T3::ptr tree = arp::T3::create(face_vert_ids, x, 24);
  calder::fast_summation<arp::T3> sum(*tree);
  sum.bind<real>(face_ix, Ac);
  std::vector<vec3> pov = {pi};
  std::vector<real> u = sum.calc<real>(
      pov,
      [](const index_t &i, const index_t &j, const vec3 &pi,
         const std::vector<calder::datum::ptr> &data,
         calder::fast_summation<arp::T3>::Node_Type node_type,
         const arp::T3 &tree) -> real {
        (void)i;
        (void)pi;
        (void)tree;
        return t3_scalar_datum<real>(node_type, j, data);
      },
      [](const index_t &i, const index_t &j, const vec3 &pi,
         const std::vector<calder::datum::ptr> &data,
         calder::fast_summation<arp::T3>::Node_Type node_type,
         const arp::T3 &tree) -> real {
        (void)i;
        (void)pi;
        (void)tree;
        return t3_scalar_datum<real>(node_type, j, data);
      },
      eps, false);
  return u[0];
}

} // namespace

GAUDI_TEST(calder_shell_pyramid_area_matches_surface_area) {
  MeshData mesh = load_sphere_mesh();
  asawa::shell::shell &M = *mesh.shell;
  std::vector<vec3> &x = mesh.vertices;

  vec3 c(0, 0, 0);
  for (const auto &v : x)
    c += v;
  c /= static_cast<real>(x.size());
  real r = 0;
  for (const auto &v : x)
    r = std::max(r, (v - c).norm());
  GAUDI_ASSERT(r > 1e-12);

  const real total_area = asawa::shell::surface_area(M, x);
  GAUDI_ASSERT(total_area > 1e-12);

  std::vector<index_t> fv = M.get_face_vert_ids();
  GAUDI_ASSERT(fv.size() >= 3);
  const vec3 x0 = x[fv[0]];
  const vec3 x1 = x[fv[1]];
  const vec3 x2 = x[fv[2]];
  const vec3 face_centroid = (x0 + x1 + x2) / 3.0;

  const vec3 p_center = c;
  const vec3 p_vertex = x[0];
  const vec3 p_face_c = face_centroid;

  const real eps_bh = 0.5;
  const real eps_tight = 1e-6;
  const real tol = 1e-4 * std::max(total_area, real(1.0));

  struct Row {
    const char *label;
    vec3 p;
    real eps;
  };
  const Row rows[] = {
      {"centroid (mesh avg)", p_center, eps_bh},
      {"centroid (mesh avg)", p_center, eps_tight},
      {"vertex x[0]", p_vertex, eps_bh},
      {"vertex x[0]", p_vertex, eps_tight},
      {"face0 centroid", p_face_c, eps_bh},
      {"face0 centroid", p_face_c, eps_tight},
  };

  std::cerr << "\n[calder_shell_pyramid_area]\n"
            << "  total_surface_area (sum of face areas) = " << std::fixed
            << std::setprecision(8) << total_area << "\n"
            << "  face_count = " << M.face_count() << "\n"
            << "  tol = " << std::scientific << tol << std::fixed << "\n\n"
            << "  " << std::setw(22) << "POV"
            << "  " << std::setw(10) << "eps"
            << "  " << std::setw(16) << "expected"
            << "  " << std::setw(16) << "recovered"
            << "  " << std::setw(14) << "|delta|"
            << "  " << std::setw(12) << "leaf_hits"
            << "\n";

  for (const Row &row : rows) {
    const real rec = shell_recovered_scalar_area(M, x, row.p, row.eps);
    const real delta = std::abs(rec - total_area);
    const real leaves =
        calder::shell_leaf_visit_count(M, x, row.p, row.eps);
    std::cerr << "  " << std::setw(22) << row.label << "  ";
    if (row.eps >= real(0.1))
      std::cerr << std::setw(10) << std::fixed << std::setprecision(2)
                << row.eps;
    else
      std::cerr << std::setw(10) << std::scientific << std::setprecision(1)
                << row.eps;
    std::cerr << std::fixed << std::setprecision(8) << "  " << std::setw(16)
              << total_area << "  " << std::setw(16) << rec << "  "
              << std::setw(14) << delta << "  " << std::setw(12)
              << static_cast<int>(leaves + real(0.5)) << "\n";
    GAUDI_EXPECT(delta < tol);
  }
  std::cerr << "  (leaf_hits = sum of fast_summation leaf callbacks with "
               "leaf->1, branch->0)\n\n";
}

GAUDI_TEST(calder_shell_leaf_visit_counts) {
  MeshData mesh = load_sphere_mesh();
  asawa::shell::shell &M = *mesh.shell;
  std::vector<vec3> &x = mesh.vertices;

  vec3 c(0, 0, 0);
  for (const auto &v : x)
    c += v;
  c /= static_cast<real>(x.size());

  std::vector<index_t> fv = M.get_face_vert_ids();
  const vec3 x0 = x[fv[0]];
  const vec3 x1 = x[fv[1]];
  const vec3 x2 = x[fv[2]];
  const vec3 face_centroid = (x0 + x1 + x2) / 3.0;

  // Opening angle 0.5 often closes whole subtrees before per-triangle leaves
  // (same as fast_winding). Use a tight eps to force leaf evaluations on-surface.
  const real eps_bh = 0.5;
  const real eps_leaf = 1e-6;
  const real n_center_bh =
      calder::shell_leaf_visit_count(M, x, c, eps_bh);
  const real n_center_tight =
      calder::shell_leaf_visit_count(M, x, c, eps_leaf);
  const real n_vertex_bh =
      calder::shell_leaf_visit_count(M, x, x[0], eps_bh);
  const real n_vertex = calder::shell_leaf_visit_count(M, x, x[0], eps_leaf);
  const real n_face_bh =
      calder::shell_leaf_visit_count(M, x, face_centroid, eps_bh);
  const real n_face = calder::shell_leaf_visit_count(M, x, face_centroid, eps_leaf);

  std::cerr << "\n[calder_shell_leaf_visit_counts]\n"
            << "  (same integrator as leaf_count.hpp: leaf=1, branch=0)\n"
            << "  " << std::setw(22) << "POV"
            << "  " << std::setw(12) << "eps=0.5"
            << "  " << std::setw(12) << "eps=1e-6"
            << "\n"
            << "  " << std::setw(22) << "centroid"
            << "  " << std::setw(12) << static_cast<int>(n_center_bh + real(0.5))
            << "  " << std::setw(12) << static_cast<int>(n_center_tight + real(0.5))
            << "\n"
            << "  " << std::setw(22) << "vertex x[0]"
            << "  " << std::setw(12) << static_cast<int>(n_vertex_bh + real(0.5))
            << "  " << std::setw(12) << static_cast<int>(n_vertex + real(0.5))
            << "\n"
            << "  " << std::setw(22) << "face0 centroid"
            << "  " << std::setw(12) << static_cast<int>(n_face_bh + real(0.5))
            << "  " << std::setw(12) << static_cast<int>(n_face + real(0.5))
            << "\n\n";

  GAUDI_EXPECT(n_vertex > 0.5);
  GAUDI_EXPECT(n_face > 0.5);
}

} // namespace test
} // namespace gaudi

#endif
