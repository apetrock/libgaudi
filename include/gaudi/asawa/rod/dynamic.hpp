#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

#include "rod.hpp"

// #include "subdivide.hpp"

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/geometry_logger.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#ifndef __ASAWA_DYNAMIC_ROD__
#define __ASAWA_DYNAMIC_ROD__

namespace gaudi {
namespace asawa {
namespace rod {
// this is ugly, but we have to do it this way, with two lists because
// the callback on the data
real vert_line(const index_t &idT, //
               const std::vector<index_t> &t_inds,
               const vector<vec3> &t_x, //
               const index_t &idS,      //
               const std::vector<index_t> &s_inds, const vector<vec3> &s_x) {

  const vec3 &xA = t_x[t_inds[idT]];

  const vec3 &xB0 = s_x[s_inds[2 * idS + 0]];
  const vec3 &xB1 = s_x[s_inds[2 * idS + 1]];
  real d = va::distance_from_line(xB0, xB1, xA);
  return d;
};

// this is ugly, but we have to do it this way, with two lists because
// the callback on the data
real line_line(const index_t &idT, //
               const std::vector<index_t> &t_inds,
               const vector<vec3> &t_x, //
               const index_t &idS,      //
               const std::vector<index_t> &s_inds, const vector<vec3> &s_x) {

  if (t_inds[2 * idT + 0] == s_inds[2 * idS + 0])
    return std::numeric_limits<real>::max();
  if (t_inds[2 * idT + 1] == s_inds[2 * idS + 1])
    return std::numeric_limits<real>::max();

  if (t_inds[2 * idT + 1] == s_inds[2 * idS + 0])
    return std::numeric_limits<real>::max();
  if (t_inds[2 * idT + 0] == s_inds[2 * idS + 1])
    return std::numeric_limits<real>::max();
  //  std::cout << t_inds[2 * idT + 0] << " " << t_inds[2 * idT + 1] << "-"
  //            << t_inds[2 * idS + 0] << " " << t_inds[2 * idS + 1] <<
  //            std::endl;
  const vec3 &xA0 = t_x[t_inds[2 * idT + 0]];
  const vec3 &xA1 = t_x[t_inds[2 * idT + 1]];

  const vec3 &xB0 = s_x[s_inds[2 * idS + 0]];
  const vec3 &xB1 = s_x[s_inds[2 * idS + 1]];

  std::array<real, 3> d = va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
  real s = d[1];
  real t = d[2];

  vec3 xA = va::mix(s, xA0, xA1);
  vec3 xB = va::mix(t, xB0, xB1);
  vec3 dA = (xA1 - xA0).normalized();
  vec3 dB = (xB1 - xB0).normalized();

  vec3 xAB = (xB - xA).normalized();
  // geometry_logger::line(xA, xB, vec4(1.0, 1.0, 1.0, 1.0));
  /*
    real s = d[1];
    real t = d[2];
    vec3 xA = (1.0 - s) * xA0 + s * xA1;
    vec3 xB = (1.0 - t) * xB0 + t * xB1;

    vec3 dAB = (xA - xB).normalized();
  */
  /*
  if (abs(dA.dot(xAB)) > 0.95)
    return std::numeric_limits<real>::infinity();
  if (abs(dB.dot(xAB)) > 0.95)
    return std::numeric_limits<real>::infinity();
   */
  if (abs(dA.dot(xAB)) > 0.35)
    return std::numeric_limits<real>::max();
  if (abs(dB.dot(xAB)) > 0.35)
    return std::numeric_limits<real>::max();
  return d[0];
};

class dynamic {
public:
  typedef std::shared_ptr<dynamic> ptr;

  using data_type = std::vector<vec3>;
  using edge_index_type = std::vector<index_t>;
  using permutation_index_type = arp::bvh_tree<2>::permutation_index_type;
  using edge_view_type = adjacency_view<data_type, edge_index_type>;
  // Use permuted_simplex_view for tuple-based access
  using permuted_edge_view_type =
      permuted_simplex_view<2, data_type, edge_index_type>;

  template <typename A> using edge_slice = slice<2, A>;
  template <typename A> using point_slice = slice<1, A>;

  static ptr create(rod::ptr R, real Cc, real Cs, real Cm) {
    ptr Rd = std::make_shared<dynamic>(R, Cc, Cs, Cm);
    return Rd;
  }

  dynamic(rod::ptr R, real Cc, real Cs, real Cm)
      : _Cc(Cc), _Cs(Cs), _Cm(Cm), __R(R) {
    edge_verts_ = __R->get_edge_vert_ids();
    edge_set_.emplace(__R->__x, edge_verts_);
    bvh_tree = arp::bvh_tree<2>::create(*edge_set_);
  }

  void set_collapse_threshold(real Cc) { _Cc = Cc; }
  void set_stretch_threshold(real Cs) { _Cs = Cs; }
  void set_bridge_threshold(real Cm) { _Cm = Cm; }

  template <typename T>
  void set(const index_t &cnew, const T &xnew, std::vector<T> &x) {
    if (cnew >= x.size()) {
      x.push_back(xnew);
    } else {
      x[cnew] = xnew;
    }
  }

  template <typename T> T avg(index_t c0, index_t c1, const std::vector<T> &x) {
    T x0 = x[c0];
    T x1 = x[c1];
    T xnew = 0.5 * (x0 + x1);
    return xnew;
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<real> &x) {
    real xnew = avg<real>(c0, c1, x);
    set<real>(cnew, xnew, x);
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<vec3> &x) {
    vec3 xnew = avg<vec3>(c0, c1, x);
    set<vec3>(cnew, xnew, x);
  }

  void interp(index_t cp, index_t c0, index_t c1, index_t c2, index_t cnew,
              std::vector<vec3> &x) {

    if (cp < 0 || c2 < 0) {
      vec3 xnew = avg<vec3>(c0, c1, x);
      set<vec3>(cnew, xnew, x);
    } else {
      vec3 xp = x[cp];
      vec3 x0 = x[c0];
      vec3 x1 = x[c1];
      vec3 x2 = x[c2];
      vec3 xnew = va::catmull_rom(xp, x0, x1, x2, 0.5);
      set<vec3>(cnew, xnew, x);
    }
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<vec4> &x) {
    vec4 xnew = avg<vec4>(c0, c1, x);
    set<vec4>(cnew, xnew, x);
  }

  void interp(index_t c0, index_t c1, index_t cnew, std::vector<quat> &x) {
    quat x0 = x[c0];
    quat x1 = x[c1];
    quat xnew = va::slerp(x0, x1, 0.5);
    // quat xnew = va::exp_slerp(x0, x1, 0.5);
    set<quat>(cnew, xnew, x);
  }

  void interp_l(index_t c0, index_t c1, index_t cnew, std::vector<real> &z) {

    vec3 x0 = __R->__x[c0];
    vec3 x1 = __R->__x[c1];
    vec3 xn = __R->__x[cnew];
    real l10 = (x1 - x0).norm();
    real ln0 = (xn - x0).norm();
    real l1n = (x1 - xn).norm();
    real lnew = ln0 + l1n;

    real l0 = z[c0];
    // real r = 0.0;
    // real l = va::mix(r, l0, lnew);
    //  std::cout << " interp: " << (ln0 + l1n) / l10 << " " << ln0 / l10 << " "
    //            << l1n / l10 << std::endl;

    // set<real>(c0, ln0 / l10 * l, z);
    // set<real>(cnew, l1n / l10 * l, z);
    set<real>(c0, ln0 / l10 * l0, z);
    set<real>(cnew, l1n / l10 * l0, z);
  }

  void split_edge(CornerId c00) {
    CornerId c10 = __R->next(c00);
    CornerId c11 = __R->next(c10);
    CornerId c01 = __R->prev(c00);
    CornerId cnew = __R->insert_edge();
    __R->link(c00, cnew);
    __R->link(cnew, c10);

    interp(c01, c00, c10, c11, cnew, __R->__x);
    // interp(c00, c10, cnew, __R->__x);
    interp(c00, c10, cnew, __R->__v);
    interp(c00, c10, cnew, __R->__u);
    interp(c00, c10, cnew, __R->__o);
    interp(c00, c10, cnew, __R->__M);
    interp(c00, c10, cnew, __R->__J);
    interp_l(c00, c10, cnew, __R->__l0);
    interp(c00, c10, cnew, __R->__t);
  }

  void collapse_edge(CornerId c00) {
    CornerId c10 = __R->next(c00);
    CornerId c11 = __R->next(c10);
    CornerId c01 = __R->prev(c00);
    if (c10 < -1)
      return;
    __R->link(c01, c10);
    __R->set_next(c00, corner_id(-1));
    __R->set_prev(c00, corner_id(-1));

    interp(c01, c10, c10, __R->__x);
    // interp(c00, c10, cnew, __R->__x);

    interp(c01, c10, c10, __R->__v);
    interp(c01, c10, c10, __R->__u);
    interp(c01, c10, c10, __R->__o);
    interp(c01, c10, c10, __R->__M);
    interp(c01, c10, c10, __R->__J);
    interp_l(c01, c10, c10, __R->__l0);
    interp(c01, c10, c10, __R->__t);
  }

#if 1
  template <Vec3View PTYPE>
  vector<std::array<index_t, 2>> get_collisions(PTYPE edges_B, real tol) {
    int num_edges = edges_B.size() / 2;
    std::vector<std::vector<std::array<index_t, 2>>> per_edge(num_edges);

#pragma omp parallel for
    for (int k = 0; k < num_edges; k++) {
      const edge_slice<PTYPE> edge(edges_B, k);
      int kk = edges_B.get_index(k);
      std::vector<index_t> neighbors =
          bvh_tree->find_neighbors(edge, tol);
      if (neighbors.empty()) {
        // Emit a sentinel so every query edge has at least one result slot;
        // downstream wrappers filter on negative ids.
        per_edge[k].push_back({-1, -1});
      } else {
        for (index_t neighbor : neighbors) {
          per_edge[k].push_back({kk, neighbor});
        }
      }
    }

    std::vector<std::array<index_t, 2>> collected;
    for (auto &pairs : per_edge) {
      collected.insert(collected.end(), pairs.begin(), pairs.end());
    }
    return collected;
  }

  vector<std::array<index_t, 4>>
  get_collisions(const std::vector<index_t> &edge_verts_B,
                 const std::vector<vec3> &x_B, real tol) {
    edge_view_type edge_view(const_cast<std::vector<vec3> &>(x_B),
                             const_cast<std::vector<index_t> &>(edge_verts_B));
    std::vector<std::array<index_t, 2>> raw = get_collisions(edge_view, tol);
    std::vector<std::array<index_t, 4>> out;
    out.reserve(raw.size());
    for (const auto &collision : raw) {
      if (collision[0] < 0 || collision[1] < 0)
        continue;
      std::array<index_t, 2> rod_ids = get_edge_ids(collision[1]);
      out.push_back({edge_verts_B[2 * collision[0] + 0],
                     edge_verts_B[2 * collision[0] + 1], rod_ids[0],
                     rod_ids[1]});
    }
    return out;
  }
#endif

#if 1
  template <Vec3View PTYPE>
  vector<std::array<index_t, 2>> get_vert_collisions(PTYPE points_B, real tol) {
    int num_points = points_B.size();
    std::vector<std::vector<std::array<index_t, 2>>> per_point(num_points);

#pragma omp parallel for
    for (int k = 0; k < num_points; k++) {
      const point_slice<PTYPE> point(points_B, k);
      std::vector<index_t> neighbors =
          bvh_tree->find_neighbors(point, tol);
      if (neighbors.empty()) {
        // Emit a sentinel so every query point has at least one result slot;
        // downstream wrappers filter on negative ids.
        per_point[k].push_back({-1, -1});
      } else {
        for (index_t neighbor : neighbors) {
          per_point[k].push_back({k, neighbor});
        }
      }
    }

    std::vector<std::array<index_t, 2>> collected;
    for (auto &pairs : per_point) {
      collected.insert(collected.end(), pairs.begin(), pairs.end());
    }
    return collected;
  }

  vector<std::array<index_t, 3>>
  get_vert_collisions(const std::vector<index_t> &point_ids,
                      const std::vector<vec3> &x_B, real tol) {
    adjacency_view<std::vector<vec3>, std::vector<index_t>> point_view(
        const_cast<std::vector<vec3> &>(x_B),
        const_cast<std::vector<index_t> &>(point_ids));
    std::vector<std::array<index_t, 2>> raw =
        get_vert_collisions(point_view, tol);
    std::vector<std::array<index_t, 3>> out;
    out.reserve(raw.size());
    for (const auto &collision : raw) {
      if (collision[0] < 0 || collision[1] < 0)
        continue;
      std::array<index_t, 2> rod_ids = get_edge_ids(collision[1]);
      out.push_back({point_ids[collision[0]], rod_ids[0], rod_ids[1]});
    }
    return out;
  }
#endif

  std::array<vec3, 2> get_edge_verts(index_t edge_id) {
    rod &R = *__R;
    if (edge_verts_.empty()) {
      edge_verts_ = R.get_edge_vert_ids();
    }
    return {R.__x[edge_verts_[2 * edge_id + 0]],
            R.__x[edge_verts_[2 * edge_id + 1]]};
  }

  std::array<index_t, 2> get_edge_ids(const index_t &i) {
    if (!edge_set_) {
      edge_verts_ = __R->get_edge_vert_ids();
      edge_set_.emplace(__R->__x, edge_verts_);
    }
    return edge_set_->tuple_ids(i);
  }

  std::array<index_t, 4> get_collision_ids(const std::array<index_t, 2> &edge_ids) {
    if (edge_ids[0] < 0 || edge_ids[1] < 0) {
      return {-1, -1, -1, -1};
    }
    const std::array<index_t, 2> &id0 = get_edge_ids(edge_ids[0]);
    const std::array<index_t, 2> &id1 = get_edge_ids(edge_ids[1]);
    return {id0[0], id0[1], id1[0], id1[1]};
  }

#if 1
  vector<std::array<index_t, 2>>
  get_internal_collisions(const real &offset = 1.0) {
    rod &R = *__R;
    std::vector<vec3> &x = R.__x;
    const std::vector<index_t> &edge_verts = R.get_edge_vert_ids();
    edge_view_type edges(x, edge_verts);
    std::vector<std::array<index_t, 2>> raw =
        get_collisions(edges, 0.5 * offset * R._r);

    std::vector<std::array<index_t, 2>> collisions;
    collisions.reserve(raw.size());
    for (auto &c : raw) {
      if (c[0] == c[1])
        continue;
      if (c[0] > c[1])
        continue;
      collisions.push_back(c);
    }
    return collisions;
  }
#endif


  void step() {

    for (int i = 0; i < __R->corner_count(); i++) {
      CornerId ci = corner_id(i);
      CornerId jn = __R->next(ci);
      CornerId jp = __R->prev(ci);
      if (jn < 0 || jp < 0)
        continue;
      vec3 q0 = __R->__x[i];
      vec3 q1 = __R->__x[jn];
      real l = (q1 - q0).norm();
      if (l > _Cs) {
        split_edge(ci);
      }

      if (l < _Cc && jp > -1) {
        collapse_edge(ci);
      }
    }
    __R->pack();
    __R->update_mass();
    __R->update_lengths();
    edge_verts_ = __R->get_edge_vert_ids();
    edge_set_.emplace(__R->__x, edge_verts_);
    bvh_tree->update(*edge_set_);
  }

  rod::ptr __R;
  std::vector<index_t> edge_verts_;
  std::optional<arp::simplex_set<2>> edge_set_;
  arp::bvh_tree<2>::ptr bvh_tree;
  real _Cc, _Cs, _Cm; // collapse, stretch, bridge
};
} // namespace rod
} // namespace asawa
} // namespace gaudi

#endif