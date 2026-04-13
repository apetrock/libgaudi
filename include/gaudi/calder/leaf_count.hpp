#ifndef GAUDI_CALDER_LEAF_COUNT_HPP
#define GAUDI_CALDER_LEAF_COUNT_HPP

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/calder/tree_code.hpp"

namespace gaudi {
namespace calder {

/// Counts `fast_summation` leaf callbacks (leaf lambda returns 1, branch 0) for
/// a triangle mesh. Same tree/bind layout as scalar pyramid tests.
inline real shell_leaf_visit_count(asawa::shell::shell &M,
                                   const std::vector<vec3> &x, const vec3 &pi,
                                   real eps = 0.5) {
  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  auto face_ids = M.get_face_range();
  std::vector<index_t> face_ix(face_ids.begin(), face_ids.end());
  std::vector<real> ones(face_ix.size(), 1.0);
  arp::T3::ptr tree = arp::T3::create(face_vert_ids, x, 24);
  fast_summation<arp::T3> sum(*tree);
  sum.bind<real>(face_ix, ones);
  std::vector<vec3> pov = {pi};
  std::vector<real> u = sum.calc<real>(
      pov,
      [](const index_t &, const index_t &, const vec3 &,
         const std::vector<datum::ptr> &,
         typename fast_summation<arp::T3>::Node_Type,
         const arp::T3 &) -> real { return 1; },
      [](const index_t &, const index_t &, const vec3 &,
         const std::vector<datum::ptr> &,
         typename fast_summation<arp::T3>::Node_Type,
         const arp::T3 &) -> real { return 0; },
      eps, false);
  return u[0];
}

/// Same for polylines / closed rods (`T2` over edges). Query positions should
/// match the vertex space used to build the tree (`rod.xc()` in integrators).
/// Note: `fast_summation` uses a 3D AABB volume for the opening test; a
/// perfectly planar curve can have zero volume at some nodes (`sc==0`), which
/// forces only the branch path—tests may need a slight out-of-plane wiggle or a
/// very small `eps` to see leaf callbacks.
inline real rod_leaf_visit_count(asawa::rod::rod &R, const vec3 &pi,
                                 real eps = 0.5) {
  std::vector<vec3> xc = R.xc();
  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<asawa::rod::CornerId> rverts = R.get_vert_range();
  std::vector<index_t> edge_ids(rverts.begin(), rverts.end());
  std::vector<real> ones(edge_ids.size(), 1.0);
  arp::T2::ptr tree = arp::T2::create(edge_verts, xc, 12);
  fast_summation<arp::T2> sum(*tree);
  sum.bind<real>(edge_ids, ones);
  std::vector<vec3> pov = {pi};
  std::vector<real> u = sum.calc<real>(
      pov,
      [](const index_t &, const index_t &, const vec3 &,
         const std::vector<datum::ptr> &,
         typename fast_summation<arp::T2>::Node_Type,
         const arp::T2 &) -> real { return 1; },
      [](const index_t &, const index_t &, const vec3 &,
         const std::vector<datum::ptr> &,
         typename fast_summation<arp::T2>::Node_Type,
         const arp::T2 &) -> real { return 0; },
      eps, false);
  return u[0];
}

} // namespace calder
} // namespace gaudi

#endif
