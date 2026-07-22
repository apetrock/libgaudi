#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/arp/arp.h"

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

//#include "GaudiGraphics/geometry_logger.h"

#include "../datums.hpp"

#include "datum_x.hpp"
#include "operations.hpp"
#include "shell.hpp"
// #include "subdivide.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/aabb.hpp"
#include "gaudi/arp/simplex_set.hpp"

#ifndef __ASAWA_DYNAMIC_SHELL__
#define __ASAWA_DYNAMIC_SHELL__
namespace gaudi {

namespace asawa {
namespace shell {

using corner1 = std::array<index_t, 1>;
using corner2 = std::array<index_t, 2>;
using corner4 = std::array<int, 4>;

using OpPredicateFcn = std::function<bool(shell &M, CornerId)>;
using MergePredicateFcn =
    std::function<bool(shell &M, CornerId, CornerId)>;

real dist_line_line(shell &M, CornerId cA0, CornerId cB0,
                    const std::vector<vec3> &x) {
  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);
  VertId vA0 = M.vert(cA0);
  VertId vA1 = M.vert(cA1);
  VertId vB0 = M.vert(cB0);
  VertId vB1 = M.vert(cB1);
  vec3 xA0 = x[vA0];
  vec3 xA1 = x[vA1];
  vec3 xB0 = x[vB0];
  vec3 xB1 = x[vB1];

  real d0 = 1.0 / 2.0 * ((xB0 - xA0).norm() + (xB1 - xA1).norm());
  return d0;
};

real dist_line_line_cen(shell &M, CornerId cA0, CornerId cB0,
                        const std::vector<vec3> &x) {
  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);
  VertId vA0 = M.vert(cA0);
  VertId vA1 = M.vert(cA1);
  VertId vB0 = M.vert(cB0);
  VertId vB1 = M.vert(cB1);
  vec3 xA0 = x[vA0];
  vec3 xA1 = x[vA1];
  vec3 xB0 = x[vB0];
  vec3 xB1 = x[vB1];

  real d0 = (0.5 * (xA1 + xA0) - 0.5 * (xB1 + xB0)).norm();
  return d0;
};

struct segment_proximity_result {
  real distance = std::numeric_limits<real>::max();
  real s = 0.0;
  real t = 0.0;
  vec3 xAB = vec3::Zero();
  bool valid = false;
};

// Full segment-segment proximity test with geometric filtering.
// Returns distance, parametric coords, and connecting vector.
// Rejects pairs that share vertices, have near-coincident endpoints,
// have closest-approach past segment ends, or whose connecting vector
// is too aligned with either edge tangent (colinearity filter).
segment_proximity_result
segment_segment_proximity(const vec3 &xA0, const vec3 &xA1, //
                          const vec3 &xB0, const vec3 &xB1, //
                          index_t vA0, index_t vA1,          //
                          index_t vB0, index_t vB1,          //
                          real tol,                           //
                          real param_margin = 0.2,            //
                          real colinear_threshold = 0.35,     //
                          real endpoint_eps_factor = 0.1) {
  segment_proximity_result result;

  // shared-vertex rejection
  if (vA0 == vB0 || vA0 == vB1 || vA1 == vB0 || vA1 == vB1)
    return result;

  // endpoint disjointness: reject if any endpoint pair is geometrically
  // near-coincident (catches separate topology with overlapping endpoints)
  real endpoint_eps = endpoint_eps_factor * tol;
  real ep2 = endpoint_eps * endpoint_eps;
  if ((xA0 - xB0).squaredNorm() < ep2 || (xA0 - xB1).squaredNorm() < ep2 ||
      (xA1 - xB0).squaredNorm() < ep2 || (xA1 - xB1).squaredNorm() < ep2)
    return result;

  std::array<real, 3> d = va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
  result.distance = d[0];
  result.s = d[1];
  result.t = d[2];

  if (result.distance > tol)
    return result;

  // parametric bounds: closest approach should be within or near the segments
  if (result.s < -param_margin || result.s > 1.0 + param_margin ||
      result.t < -param_margin || result.t > 1.0 + param_margin)
    return result;

  vec3 xA = va::mix(result.s, xA0, xA1);
  vec3 xB = va::mix(result.t, xB0, xB1);
  vec3 dv = xB - xA;
  real dv_norm = dv.norm();
  if (dv_norm < 1e-12)
    return result;

  result.xAB = dv / dv_norm;

  // colinearity filter: connecting vector should be roughly perpendicular
  // to both edge directions (rejects adjacent / colinear edges)
  vec3 dA = (xA1 - xA0).normalized();
  vec3 dB = (xB1 - xB0).normalized();
  if (std::abs(dA.dot(result.xAB)) > colinear_threshold)
    return result;
  if (std::abs(dB.dot(result.xAB)) > colinear_threshold)
    return result;

  result.valid = true;
  return result;
}

std::mt19937_64 rng;
std::uniform_real_distribution<real> unif(0.0, 1.0);

void debug_line(shell &M,           //
                const CornerId &cA0, //
                const vector<vec3> &x) {

  CornerId cA1 = M.other(cA0);
  VertId vA0 = M.vert(cA0);
  VertId vA1 = M.vert(cA1);

  const vec3 &ca0 = x[vA0];
  const vec3 &ca1 = x[vA1];

  vec4 c(unif(rng), unif(rng), unif(rng), 1.0);
  geometry_logger::line(ca0, ca1, c);
};

void debug_line_line(shell &M,           //
                     const CornerId &cA0, //
                     const CornerId &cB0, //
                     const vector<vec3> &x) {

  CornerId cA1 = M.other(cA0);
  CornerId cB1 = M.other(cB0);
  VertId vA0 = M.vert(cA0);
  VertId vA1 = M.vert(cA1);
  VertId vB0 = M.vert(cB0);
  VertId vB1 = M.vert(cB1);

  const vec3 &ca0 = x[vA0];
  const vec3 &ca1 = x[vA1];
  const vec3 &cb0 = x[vB0];
  const vec3 &cb1 = x[vB1];

  vec4 c(unif(rng), unif(rng), unif(rng), 1.0);
  geometry_logger::line(ca0, ca1, c);
  geometry_logger::line(cb0, cb1, c);

  geometry_logger::line(0.5 * (ca0 + ca1), 0.5 * (cb0 + cb1), c);
};

void debug_edge_normal(shell &M,          //
                       const CornerId &c0, //
                       const vector<vec3> &x) {

  vec3 N = edge_normal(M, c0, x);
  vec3 cen = edge_center(M, c0, x);
  vec4 col(unif(rng), unif(rng), unif(rng), 1.0);
  geometry_logger::line(cen, cen + 0.1 * N, col);
};

template <int OP, int C_ALLOC, int V_ALLOC, int F_ALLOC, int MSIZE>
std::vector<std::array<index_t, MSIZE>> //
op_edges(shell &M,                      //
         std::vector<index_t> &edges_to_op, std::vector<real> &S,
         const std::vector<vec3> &x,
         std::function<
             std::array<index_t, MSIZE>(index_t i,                         //
                                        index_t cs,                        //
                                        index_t vs,                        //
                                        index_t fs,                        //
                                        const std::vector<index_t> &edges, //
                                        shell &m)>
             func) {
  int STRIDE = OP < 2 ? 1 : 2;
  size_t cstart = M.corner_count();
  size_t estart = M.edge_count();
  size_t vstart = M.vert_count();
  size_t fstart = M.face_count();
  size_t Ne = edges_to_op.size();
  // Guard against size_t wrap / absurd counts that become std::length_error.
  constexpr size_t kMaxOpEdges = 1ull << 24;
  if (Ne > kMaxOpEdges) {
    std::cerr << "op_edges: absurd edge count Ne=" << Ne
              << " OP=" << OP << " C_ALLOC=" << C_ALLOC
              << " V_ALLOC=" << V_ALLOC << " F_ALLOC=" << F_ALLOC << std::endl;
    std::abort();
  }
  const size_t c_alloc_n = static_cast<size_t>(C_ALLOC) * Ne;
  const size_t v_alloc_n = static_cast<size_t>(V_ALLOC) * Ne;
  const size_t f_alloc_n = static_cast<size_t>(F_ALLOC) * Ne;
  M.inflate_edge_pairs(c_alloc_n);
  M.inflate_verts(v_alloc_n);
  M.inflate_faces(f_alloc_n);

  for (auto d : M.get_data()) {

    if (d->type() == EDGE)
      d->alloc(c_alloc_n);
    if (d->type() == VERTEX)
      d->alloc(v_alloc_n);
    if (d->type() == FACE)
      d->alloc(f_alloc_n);
  }

  for (index_t i = 0; i < edges_to_op.size(); i += STRIDE) {
    for (auto d : M.get_data()) {

      if (M.next(corner_id(edges_to_op[i])) < 0)
        continue;

      if (d->type() == EDGE) {
        if (OP == 0) {
          real s0 = S[i];
          real s1 = 1.0 - s0;
          index_t c0 = edges_to_op[i];
          index_t e0 = estart + 3 * i + 0;
          index_t e1 = estart + 3 * i + 1;
          index_t e2 = estart + 3 * i + 2;
          d->subdivide(M, c0 / 2, e0, e1, e2, s0, c0);
        }
        if (OP == 1) {
          index_t c0 = edges_to_op[i];
          d->collapse(M, c0);
        }
      }

      if (d->type() == FACE) {
        real s0 = S[i];
        real s1 = 1.0 - s0;

        index_t c0 = edges_to_op[i];
        CornerId c0i = corner_id(c0);
        CornerId c1 = M.other(c0i);
        FaceId f0 = M.face(c0i);
        FaceId f1 = M.face(c1);
        real a0 = face_area(M, f0, x);
        real a1 = face_area(M, f1, x);
        if (OP == 0) {
          index_t fs0 = fstart + 2 * i + 0;
          index_t fs1 = fstart + 2 * i + 1;
          d->subdivide(M, f0, f1, fs0, fs1, s0, c0);
        } else if (OP == 1) {
          d->collapse(M, c0);
        }
      }

      if (d->type() == VERTEX) {
        if (OP == 0) {
          // subd
          index_t c0 = edges_to_op[i];
          d->subdivide(M, vstart + i, -1, -1, -1, S[i], c0);
        } else if (OP == 1) {
          // collapse
          index_t c0 = edges_to_op[i];
          d->collapse(M, c0);
        } else if (OP == 2) {
          // merge
          index_t cA0 = edges_to_op[i + 0];
          CornerId cA0i = corner_id(cA0);
          CornerId cA1 = M.other(cA0i);
          index_t cB0 = edges_to_op[i + 1];
          CornerId cB0i = corner_id(cB0);
          CornerId cB1 = M.other(cB0i);
          index_t vs0 = vstart + 2 * i + 0;
          index_t vs1 = vstart + 2 * i + 1;
          d->merge(M,                        //
                   M.vert(cA0i), M.vert(cA1), //
                   vs0, vs1,                 //
                   M.vert(cA0i), M.vert(cA1), //
                   M.vert(cB0i), M.vert(cB1));
        }
      }
    }
  }

  // #pragma omp parallel for shared(post_edges)
  //  before parallelize, need to preallocate;
  std::vector<std::array<index_t, MSIZE>> collection(edges_to_op.size());

  for (index_t i = 0; i < edges_to_op.size(); i += STRIDE) {
    index_t ic0 = edges_to_op[i];
    // #pragma omp critical
    if (M.next(corner_id(ic0)) < 0) {
      collection[i].fill(-1);
      continue;
    }

    auto p = func(i, cstart, vstart, fstart, edges_to_op, M);
    collection[i] = p;
  }
  return collection;
}

// template <int STRIDE, int C_ALLOC, int V_ALLOC, int F_ALLOC, int MSIZE>
auto subdivide_op = op_edges<0, 3, 1, 2, 1>;
auto collapse_op = op_edges<1, 0, 0, 0, 1>;
auto merge_op = op_edges<2, 0, 2, 0, 4>;

void subdivide_edges(shell &M) {
  std::vector<index_t> edges_to_divide;
  const std::vector<vec3> &x = get_vec_data(M, 0);
  edges_to_divide.push_back(3);
  edges_to_divide.push_back(5);
  edges_to_divide.push_back(8);
  edges_to_divide.push_back(13);
  edges_to_divide.push_back(21);
  std::vector<real> S(edges_to_divide.size(), 0.5);
  subdivide_op(M, edges_to_divide, S, x,
               [](index_t i,  //
                  index_t cs, //
                  index_t vs, //
                  index_t fs, //
                  const std::vector<index_t> &edges, shell &m) -> corner1 {
                 return {subdivide_edge(m, corner_id(edges[i]))};
               });
}

void collapse_edges(shell &M) {
  std::vector<index_t> edges_to_divide;
  const std::vector<vec3> &x = get_vec_data(M, 0);
  edges_to_divide.push_back(4);
  // edges_to_divide.push_back(6);
  // edges_to_divide.push_back(9);
  //  edges_to_divide.push_back(14);
  edges_to_divide.push_back(22);

  std::vector<real> S(edges_to_divide.size(), 0.5);
  collapse_op(M, edges_to_divide, S, x,
              [](index_t i,  //
                 index_t cs, //
                 index_t vs, //
                 index_t fs, //
                 const std::vector<index_t> &edges,
                 shell &m) -> corner1 {
                return {collapse_edge(m, corner_id(edges[i]))};
              });
}

real length(CornerId c0, CornerId c1, const shell &M,
            const std::vector<vec3> &data) {
  vec3 v0 = data[M.vert(c0)];
  vec3 v1 = data[M.vert(c1)];
  return (v1 - v0).norm();
}

template <typename T, typename comp> class shell_data_comp {
public:
  shell_data_comp(const std::vector<T> &data, real eps, const shell &M)
      : _M(M), _data(data), _eps(eps) {}
  bool operator()(index_t c0, index_t c1) const {
    real dv = length(corner_id(c0), corner_id(c1), _M, _data);
    // std::cout << "comp: " << dv << " " << _eps << std::endl;
    return comp{}(dv, _eps);
  }

  std::vector<CornerId> get_edges() const {

    auto edges = _M.get_edge_range();
    std::vector<real> lengths = edge_lengths(_M, _data);
    sort(edges.begin(), edges.end(),
         [this, &lengths](CornerId ca, CornerId cb) -> bool {
           real dva = lengths[ca / 2];
           real dvb = lengths[cb / 2];
           return comp{}(dva, dvb);
         });

    return edges;
  }

  const std::vector<T> &_data;
  const shell &_M;
  real _eps;
};

struct face_dedup {
  std::vector<bool> face_flags;

  explicit face_dedup(const shell &M)
      : face_flags(static_cast<size_t>(M.face_count()), false) {}

  bool is_blocked(const shell &M, CornerId c0, CornerId c1) const {
    FaceId f0 = M.face(c0);
    FaceId f1 = M.face(c1);
    if (f0 < 0 || f1 < 0)
      return true;
    if (face_flags[static_cast<size_t>(f0)])
      return true;
    if (face_flags[static_cast<size_t>(f1)])
      return true;
    return false;
  }

  void mark(shell &M, CornerId c0, CornerId c1) {
    FaceId f0 = M.face(c0);
    FaceId f1 = M.face(c1);
    face_flags[static_cast<size_t>(f0)] = true;
    face_flags[static_cast<size_t>(f1)] = true;
  }
};

struct one_ring_dedup {
  std::vector<bool> edge_flags;

  explicit one_ring_dedup(const shell &M)
      : edge_flags(static_cast<size_t>(M.corner_count()) / 2, false) {}

  bool is_blocked(const shell & /*M*/, CornerId c0, CornerId /*c1*/) const {
    return edge_flags[static_cast<size_t>(c0) / 2];
  }

  void mark(shell &M, CornerId c0, CornerId c1) {
    VertId v0 = M.vert(c0);
    VertId v1 = M.vert(c1);
    if (v0 < 0 || v1 < 0)
      return;
    edge_flags[static_cast<size_t>(c0) / 2] = true;
    flag_one_ring(M, v0);
    flag_one_ring(M, v1);
  }

private:
  void flag_one_ring(shell &M, VertId v) {
    M.for_each_vertex(v, [this](CornerId ci, shell &m) {
      CornerId ca = ci;
      CornerId cb = m.next(ci);
      CornerId cc = m.next(cb);
      edge_flags[static_cast<size_t>(ca) / 2] = true;
      edge_flags[static_cast<size_t>(cb) / 2] = true;
      edge_flags[static_cast<size_t>(cc) / 2] = true;
    });
  }
};

template <typename comparator, typename DedupPolicy = one_ring_dedup>
std::vector<CornerId> gather_edges(shell &M, const comparator &comp) {
  auto edges = comp.get_edges();
  DedupPolicy policy(M);

  std::vector<CornerId> edges_out;
  edges_out.reserve(edges.size());

  for (int i = 0; i < static_cast<int>(edges.size()); i++) {
    CornerId c0 = edges[static_cast<size_t>(i)];
    CornerId c1 = M.other(c0);

    if (policy.is_blocked(M, c0, c1))
      continue;

    if (comp(c0, c1)) {
      edges_out.push_back(c0);
      policy.mark(M, c0, c1);
    }
  }

  return edges_out;
}
#if 0 
void smoothMesh(shell &M, real C, int N) {

  // return;
  vertex_array &vertices = this->_surf->get_vertices();
  int i = 0;
  // asawa::area_laplacian_0<SPACE, coordinate_type> M(this->_surf);
  asawa::laplacian3<SPACE> M(this->_surf);

  std::cout << "ugly smoothing " << std::flush;
  coordinate_array coords = asawa::ci::get_coordinates<SPACE>(this->_surf);
  coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(this->_surf);

  std::vector<real> sm =
      ci::get<SPACE, real>(_surf, SPACE::vertex_index::SMOOTH);
  for (int k = 0; k < N; k++) {
    std::cout << "." << std::flush;
    M.build();
    coords = M.smooth(coords, C, C + 3e-5);
  }

  asawa::ci::set_coordinates<SPACE>(coords, this->_surf);
  std::cout << "done!" << std::endl;
}
#endif

class dynamic {
public:
  typedef std::shared_ptr<dynamic> ptr;

  template <typename A> using edge_slice = slice<2, A>;
  template <typename A> using point_slice = slice<1, A>;

  static ptr create(shell::ptr M, real Cc, real Cs, real Cm) {
    return std::make_shared<dynamic>(M, Cc, Cs, Cm);
  }

  dynamic(shell::ptr M, real Cc, real Cs, real Cm) : __M(M) {
    _Cc = Cc;
    _Cs = Cs;
    _Cm = Cm;

    std::vector<vec3> velocities(__M->vert_count(), vec3::Zero());
    datum_t<vec3>::ptr vdata =
        datum_t<vec3>::create(prim_type::VERTEX, velocities);
    __vdatum_id = __M->insert_datum(vdata);
  };

  void set_collapse_threshold(real Cc) { _Cc = Cc; }
  void set_stretch_threshold(real Cs) { _Cs = Cs; }
  void set_bridge_threshold(real Cm) { _Cm = Cm; }

  /// Rebuild edge and face BVH trees from the current vertex positions (call
  /// once per frame or before collision queries after geometry changes).
  void update_trees() {
    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    std::vector<index_t> edge_verts = __M->get_edge_vert_ids();
    std::vector<index_t> face_verts = __M->get_face_vert_ids(true);
    // simplex_set is only a temporary adapter; bvh_tree copies positions +
    // adjacency into its own members.
    arp::simplex_set<2> edge_set(x, edge_verts);
    arp::simplex_set<3> face_set(x, face_verts);
    if (!edge_tree_)
      edge_tree_ = arp::bvh_tree<2>::create(edge_set);
    else
      edge_tree_->update(edge_set);
    if (!face_tree_)
      face_tree_ = arp::bvh_tree<3>::create(face_set);
    else
      face_tree_->update(face_set);
  }

  /// True if the shell still has at least one live edge (not collapsed to nothing).
  bool has_active_edges() const {
    for (int i = 0; i < static_cast<int>(__M->corner_count()); i += 2) {
      if (__M->next(corner_id(i)) >= 0)
        return true;
    }
    return false;
  }

  void delete_degenerates(shell &M) {
    vec3_datum::ptr x_datum = static_pointer_cast<vec3_datum>(M.get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    for (int i = 0; i < static_cast<int>(M.corner_count()); i++) {
      CornerId ci = corner_id(i);
      if (M.next(ci) < 0)
        continue;
      CornerId c0 = ci;
      CornerId c1 = M.other(ci);
      if (M.vert(c0) != M.vert(c1))
        continue;
      for (auto d : M.get_data()) {
        d->collapse(M, c0);
      }
      collapse_edge(M, c0, true);
    }

    for (int i = 0; i < static_cast<int>(M.face_count()); i++) {
      FaceId fi = face_id(i);
      if (M.fbegin(fi) < 0)
        continue;

      real a0 = face_area(M, fi, x);
      if (a0 < 1e-10) {
      }

      if (M.fsize(fi) > 2)
        continue;
      CornerId c0 = M.fbegin(fi);
      CornerId c1 = M.other(c0);
      for (auto d : M.get_data()) {
        d->collapse(M, c0);
      }
      merge_face(M, c0, c1);
    }

    for (int i = 0; i < static_cast<int>(M.vert_count()); i++) {
      VertId vi = vert_id(i);
      if (M.vbegin(vi) < 0)
        continue;
      /*
            if (M.vsize(i) > 16) {
              vec3 N = vert_normal(M, i, x);
              vec4 cola(0.0, 1.0, 1.0, 0.0);
              logger::line(x[i], x[i] + 0.1 * N, cola);
              M.for_each_vertex(i, [&x, &i, cola](index_t cid, shell &m) {
                index_t j = m.vert(m.next(cid));
                logger::line(x[i], x[j], cola);
              });

            }
      */
      if (M.vsize(vi) > 3)
        continue;

      vec3 N = vert_normal(M, vi, x);
      remove_vertex(M, vi);
    }
  }

  index_t align_edges(shell &M, index_t cA0, index_t cB0,
                      const std::vector<vec3> &x) {

    CornerId cB0i = corner_id(cB0);
    CornerId cB1 = M.other(cB0i);
    real d0 = dist_line_line(M, corner_id(cA0), cB0i, x);
    real d1 = dist_line_line(M, corner_id(cA0), cB1, x);

    if (d1 < d0)
      return cB1;

    return cB0;
  }

  void
  trim_edge_edge_collected(shell &M, const std::vector<vec3> &x,
                           std::vector<std::array<index_t, 2>> &collected) {
    std::vector<bool> flags(M.corner_count() / 2, false);
    real tol = _Cm;
    collected.erase(
        std::remove_if(
            collected.begin(), collected.end(),
            [tol, &M, &x, &flags](const auto &p) {
              if (p[0] < 0 || p[1] < 0)
                return true;

              CornerId cA0 = corner_id(p[0]);
              CornerId cB0 = corner_id(p[1]);
              CornerId cA1 = M.other(cA0);
              CornerId cB1 = M.other(cB0);

              const vec3 &xA0 = x[M.vert(cA0)];
              const vec3 &xA1 = x[M.vert(cA1)];
              const vec3 &xB0 = x[M.vert(cB0)];
              const vec3 &xB1 = x[M.vert(cB1)];

              // proper segment-segment distance
              std::array<real, 3> d =
                  va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
              if (d[0] > tol)
                return true;

              // connecting vector between closest points
              vec3 xA = va::mix(d[1], xA0, xA1);
              vec3 xB = va::mix(d[2], xB0, xB1);
              vec3 dv = xB - xA;
              real dv_norm = dv.norm();
              if (dv_norm < 1e-12)
                return true;
              vec3 xAB = dv / dv_norm;

              // normals must oppose each other
              vec3 NA = edge_normal(M, cA0, x);
              vec3 NB = edge_normal(M, cB0, x);
              if (va::dot(NA, NB) > -0.0)
                return true;

              // connecting vector must go outward from A and inward to B
              if (va::dot(NA, xAB) < 0.0)
                return true;
              if (va::dot(NB, xAB) > 0.0)
                return true;

              // dedup: one merge per edge
              if (flags[p[0] / 2] || flags[p[1] / 2])
                return true;

              flags[p[0] / 2] = true;
              flags[p[1] / 2] = true;

              return false;
            }),
        collected.end());
  }

  vector<std::array<index_t, 2>>
  get_edge_edge_collisions(const std::vector<index_t> &edge_verts_t, //
                           const std::vector<index_t> &edge_map_t,
                           const std::vector<vec3> &x_t, real tol) {
    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x_m = x_datum->data();

    std::vector<index_t> edge_verts_m = __M->get_edge_vert_ids();
    std::vector<index_t> edge_map_m = __M->get_edge_map();

    update_trees();

    adjacency_view<std::vector<vec3>, std::vector<index_t>> edge_view(
        const_cast<std::vector<vec3> &>(x_t),
        const_cast<std::vector<index_t> &>(edge_verts_t));

    std::vector<std::array<index_t, 2>> collected(edge_verts_t.size() / 2,
                                                  {-1, -1});
#pragma omp parallel for
    for (int i = 0; i < edge_verts_t.size(); i += 2) {
      index_t e0 = i / 2;
      index_t c0 = edge_map_t[e0];
      collected[i / 2] = {c0, -1};

      const edge_slice<decltype(edge_view)> edge(edge_view, e0);
      std::vector<index_t> nbrs = edge_tree_->find_neighbors(edge, tol);

      index_t best_e = -1;
      real best_d = std::numeric_limits<real>::max();
      for (index_t e1 : nbrs) {
        if (e0 >= e1)
          continue;

        index_t vT0 = edge_verts_t[2 * e0 + 0];
        index_t vT1 = edge_verts_t[2 * e0 + 1];
        index_t vS0 = edge_verts_m[2 * e1 + 0];
        index_t vS1 = edge_verts_m[2 * e1 + 1];

        auto pr = segment_segment_proximity(
            x_t[vT0], x_t[vT1], x_m[vS0], x_m[vS1], //
            vT0, vT1, vS0, vS1, tol);

        if (pr.valid && pr.distance < best_d) {
          best_d = pr.distance;
          best_e = e1;
        }
      }
      if (best_e >= 0 && c0 >= 0) {
        index_t c1 = edge_map_m[best_e];
        collected[i / 2] = {c0, c1};
      }
    }
    return collected;
  }

  vector<std::array<index_t, 2>>
  get_pnt_tri_collisions(const std::vector<index_t> &verts_t,
                         const std::vector<index_t> &verts_map_t,
                         const std::vector<vec3> &x_t, shell &M, real tol) {

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x_m = x_datum->data();
    std::vector<index_t> face_verts_m = __M->get_face_vert_ids(true);
    std::vector<index_t> face_map_m = __M->get_face_map(true);

    (void)M;
    update_trees();

    adjacency_view<std::vector<vec3>, std::vector<index_t>> point_view(
        const_cast<std::vector<vec3> &>(x_t),
        const_cast<std::vector<index_t> &>(verts_t));

    std::vector<std::array<index_t, 2>> collected(verts_t.size(), {-1, -1});
#pragma omp parallel for
    for (int i = 0; i < verts_t.size(); i++) {
      index_t iv = verts_map_t[i];
      collected[i] = {iv, -1};

      const point_slice<decltype(point_view)> point(point_view, i);
      std::vector<index_t> nbrs = face_tree_->find_neighbors(point, tol);

      index_t best_f = -1;
      real best_d = std::numeric_limits<real>::max();
      for (index_t fi : nbrs) {
        real d = arp::pnt_tri_min(i, verts_t, x_t, fi, face_verts_m, x_m);
        if (d < tol && d < best_d) {
          best_d = d;
          best_f = fi;
        }
      }
      if (best_f >= 0) {
        collected[i] = {iv, face_map_m[best_f]};
      }
    }

    return collected;
  }

  vector<std::array<index_t, 2>> get_internal_edge_edge_collisions(real tol) {
    shell &M = *__M;
    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    std::vector<index_t> edge_verts = __M->get_edge_vert_ids();
    std::vector<index_t> edge_map = __M->get_edge_map();

    std::vector<std::array<index_t, 2>> collected =
        get_edge_edge_collisions(edge_verts, edge_map, x, tol);

    for (auto &c : collected) {
      //std::cout << "aligning " << c[0] << " " << c[1] << std::endl;
      if (c[0] > -1 && c[1] > -1) {
        //std::cout << "aligning " << c[0] << " " << c[1] << std::endl;
        c[1] = align_edges(M, c[0], c[1], x);
      }
    }
    return collected;
  }

  vector<std::array<index_t, 2>> get_internal_pnt_tri_collisions(real tol) {
    shell &M = *__M;
    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    auto verts_typed = __M->get_vert_range();
    std::vector<index_t> verts(verts_typed.begin(), verts_typed.end());
    std::vector<index_t> verts_map = __M->get_vert_map();

    return get_pnt_tri_collisions(verts, verts_map, x, M, tol);
  }

  void merge_edges() {
    // edge e = c / 2;

    shell &M = *__M;
    if (!has_active_edges()) {
      std::cout << "merge: skip (no active edges)" << std::endl;
      return;
    }

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    // edge_tree_->debug();
    real tol = 0.25 * this->_Cm * this->_Cm;
    auto collected = get_internal_edge_edge_collisions(tol);
    trim_edge_edge_collected(M, x, collected);

    if (_merge_pred)
      collected.erase(std::remove_if(collected.begin(), collected.end(),
                                     [this](auto c) { return _merge_pred(*__M, corner_id(c[0]), corner_id(c[1])); }),
                      collected.end());

    for (const auto &cp : collected) {
      if (cp[0] < 0 || cp[1] < 0)
        continue;
      vec3 cenA = edge_center(M, corner_id(cp[0]), x);
      vec3 cenB = edge_center(M, corner_id(cp[1]), x);
      geometry_logger::line(cenA, cenB, vec4(0.0, 1.0, 0.1, 1.0));
    }

    std::vector<index_t> f_collect(2 * collected.size());
    for (int i = 0; i < collected.size(); i++) {
      f_collect[2 * i + 0] = collected[i][0];
      f_collect[2 * i + 1] = collected[i][1];
    }
    std::cout << "merge: pairs=" << collected.size()
              << " f_collect=" << f_collect.size()
              << " verts=" << M.vert_count()
              << " corners=" << M.corner_count() << std::endl;
    std::vector<real> S(f_collect.size(), 0.5);
    try {
               merge_op(*__M, f_collect, S, x,
             [&x, tol](index_t i,                         //
                       index_t cs,                        //
                       index_t vs,                        //
                       index_t fs,                        //
                       const std::vector<index_t> &edges, //
                       shell &M) -> corner4 {
               CornerId c0A = corner_id(edges[i + 0]);
               CornerId c0B = corner_id(edges[i + 1]);

               if (M.next(c0A) < 0 || M.next(c0B) < 0) {
                 return {-1, -1, -1, -1};
               }

               return merge_edge(M, c0A, c0B, vert_id(vs + 2 * i + 0),
                                 vert_id(vs + 2 * i + 1));
             });
    } catch (const std::exception &e) {
      std::cerr << "merge_op threw: " << e.what()
                << " pairs=" << collected.size()
                << " f_collect=" << f_collect.size()
                << " verts=" << M.vert_count()
                << " corners=" << M.corner_count() << std::endl;
      throw;
    }
  }

  void break_cycles() {
    // edge e = c / 2;
    shell &M = *__M;

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    auto edges = M.get_edge_range();
    std::vector<std::array<index_t, 2>> collected;

    for (int i = 0; i < static_cast<int>(edges.size()); i++) {
      CornerId c0 = corner_id(edges[static_cast<size_t>(i)]);
      CornerId c1 = M.other(c0);
      VertId v1 = M.vert(c1);
      if (count_cycle(M, c0) < 2)
        continue;

      std::array<index_t, 2> pair = {0, 0};
      int j = 0;
      M.for_each_vertex(M.vert(c0), [v1, &pair, &j](CornerId ci, shell &M) {
        VertId vi = M.vert(M.next(ci));
        if (vi == v1 && j < 2) {
          pair[j++] = ci;
        }
      });
      collected.push_back(pair);
    }

    if (_merge_pred)
      collected.erase(std::remove_if(collected.begin(), collected.end(),
                                     [this](auto c) { return _merge_pred(*__M, corner_id(c[0]), corner_id(c[1])); }),
                      collected.end());

    std::vector<index_t> f_collect(2 * collected.size());
    for (int i = 0; i < collected.size(); i++) {
      f_collect[2 * i + 0] = collected[i][0];
      f_collect[2 * i + 1] = collected[i][1];
    }

    std::vector<real> S(f_collect.size(), 0.5);
    merge_op(*__M, f_collect, S, x,
             [&x](index_t i,                         //
                  index_t cs,                        //
                  index_t vs,                        //
                  index_t fs,                        //
                  const std::vector<index_t> &edges, //
                  shell &M) -> corner4 {
               CornerId c0A = corner_id(edges[i + 0]);
               CornerId c0B = corner_id(edges[i + 1]);

               return merge_edge(M, c0A, c0B, vert_id(vs + 2 * i + 0),
                                 vert_id(vs + 2 * i + 1));
             });
  }

  bool skip_flip(shell &M, CornerId corner) {
    CornerId c0 = corner;
    CornerId c1 = M.other(c0);
    VertId v0 = M.vert(c0);
    VertId v1 = M.vert(c1);

    if (M.vsize(v0) < 3)
      return true;
    if (M.vsize(v1) < 3)
      return true;

    if (M.vert(M.prev(c0)) == M.prev(c1))
      return true;

    if (M.vert(M.prev(c0)) == M.vert(M.prev(c1)))
      return true;
    return false;
  }

  void flip_edges() {

    auto edges = __M->get_edge_range();
    for (int i = 0; i < edges.size(); i++) {
      int card = rand() % edges.size();
      std::swap(edges[i], edges[card]);
    }

    if (_flip_pred)
      edges.erase(std::remove_if(edges.begin(), edges.end(),
                     [this](CornerId c) { return _flip_pred(*__M, c); }), edges.end());

    for (int i = 0; i < static_cast<int>(edges.size()); i++) {
      CornerId c0 = corner_id(edges[static_cast<size_t>(i)]);
      CornerId c1 = __M->prev(c0);

      if (skip_flip(*__M, c0))
        continue;

      CornerId c2 = __M->other(c0);
      CornerId c3 = __M->prev(c2);

      vec3_datum::ptr coord_datum =
          static_pointer_cast<vec3_datum>(__M->get_datum(0));
      const std::vector<vec3> &data = coord_datum->data();
      vec3 v0 = data[__M->vert(c0)];
      vec3 v1 = data[__M->vert(c1)];
      vec3 v2 = data[__M->vert(c2)];
      vec3 v3 = data[__M->vert(c3)];

#if 0 // volume guard — set to 0 to disable
      {
        vec3 e01 = v1 - v0;
        vec3 e02 = v2 - v0;
        vec3 e03 = v3 - v0;
        real vol = std::abs(e03.dot(e01.cross(e02))) / 6.0;
        real l_avg = 0.25 * (e01.norm() + (v2 - v1).norm() +
                             (v3 - v2).norm() + e03.norm());
        real l3 = l_avg * l_avg * l_avg;
        if (l3 > 1e-20 && vol / l3 > _max_flip_vol_ratio)
          continue;
      }
#endif

      real m01 = 1.0 / (v0 - v1).norm();
      real m12 = 1.0 / (v1 - v2).norm();
      real m23 = 1.0 / (v2 - v3).norm();
      real m30 = 1.0 / (v3 - v0).norm();

      real cos0 = (v1 - v0).dot(v3 - v0) * m01 * m30;
      real cos1 = (v0 - v1).dot(v2 - v1) * m01 * m12;
      real cos2 = (v1 - v2).dot(v3 - v2) * m12 * m23;
      real cos3 = (v0 - v3).dot(v2 - v3) * m30 * m23;
      // half angle cos^2(2a) = 0.5*(1+cos(a))
      real cSame = acos(cos1) + acos(cos3); // corresponds to flipped edge
      real cFlip = acos(cos0) + acos(cos2); // corresponds to flipped edge
                                            // surface angles

      // current normals
      vec3 N00 = va::calculate_normal(v1, v0, v2);
      vec3 N01 = va::calculate_normal(v3, v2, v0);
      // new normals
      vec3 N10 = va::calculate_normal(v0, v3, v1);
      vec3 N11 = va::calculate_normal(v2, v1, v3);

      /*
      real cosN0 = va::dot(N00, N01);
      real tSame = M_PI - acos(cosN0);
      real cosN1 = va::dot(N10, N11);
      real tFlip = M_PI - acos(cosN1);
      real nFlip = tFlip;
      */
      real cosN0 = va::norm(vec3(N01 - N00));
      real sinN0 = va::norm(vec3(N01 + N00));
      real tSame = atan2(sinN0, cosN0);
      real cosN1 = va::norm(vec3(N11 - N10));
      real sinN1 = va::norm(vec3(N11 + N10));
      real tFlip = atan2(sinN1, cosN1);

#if 0 // dihedral guard — set to 0 to disable
      if (tFlip > _max_flip_dihedral)
        continue;
#endif

      real dt = tFlip - tSame;
      // std::cout << tFlip << " " << tSame << " " << dt << std::endl;
      real eFlip = cFlip * cFlip + 10.0 * dt * dt;
      real eSame = cSame * cSame;
      // real eFlip = cFlip * cFlip + tFlip * tFlip;
      // real eSame = cSame * cSame + tSame * tSame;

      // std::cout << eFlip << " " << eSame << std::endl;
      //  if (false) {
      if (eFlip < 1.0 * eSame) {
        for (auto d : __M->get_data()) {
          if (d->type() == EDGE) {
            d->flip(*__M, c0);
          }
        }

        flip_edge(*__M, c0);
      }
    }
  }

  void subdivide_edges() {
    using comp_great = shell_data_comp<vec3, std::greater<real>>;
    const std::vector<vec3> &x = get_vec_data(*__M, 0);
    auto cmp = comp_great(x, _Cs, *__M);

    auto edges_typed = gather_edges<comp_great>(*__M, cmp);
    std::vector<index_t> edges_to_divide(edges_typed.begin(), edges_typed.end());

    std::vector<real> S(edges_to_divide.size(), 0.5);
    subdivide_op(*__M, edges_to_divide, S, x,
                 [](index_t i,  //
                    index_t cs, //
                    index_t vs, //
                    index_t fs, //
                    const std::vector<index_t> &edges, shell &m) -> corner1 {
                   return {subdivide_edge(m, corner_id(edges[i]),    //
                                          vert_id(vs + i),         //
                                          corner_id(cs + 6 * i + 0), //
                                          corner_id(cs + 6 * i + 2), //
                                          corner_id(cs + 6 * i + 4), //
                                          face_id(fs + 2 * i + 0), //
                                          face_id(fs + 2 * i + 1)  //
                                          )};
                 });
  }

  void collapse_edges() {

    using comp = shell_data_comp<vec3, std::less<real>>;
    const std::vector<vec3> &x = get_vec_data(*__M, 0);
    auto cmp = comp(x, _Cc, *__M);

    auto edges_typed = gather_edges<comp>(*__M, cmp);
    std::vector<index_t> edges_to_divide(edges_typed.begin(), edges_typed.end());

    if (_collapse_pred)
      edges_to_divide.erase(std::remove_if(edges_to_divide.begin(), edges_to_divide.end(),
                     [this](index_t c) { return _collapse_pred(*__M, corner_id(c)); }), edges_to_divide.end());

    std::vector<real> S(edges_to_divide.size(), 0.5);
    collapse_op(*__M, edges_to_divide, S, x,
                [this](index_t i,  //
                       index_t cs, //
                       index_t vs, //
                       index_t fs, //
                       const std::vector<index_t> &edges, shell &m) -> corner1 {
                  if (_collapse_pred && _collapse_pred(m, corner_id(edges[i])))
                    return {-1};
                  return {collapse_edge(m, corner_id(edges[i]))};
                });
  }

  void test_nan() {
    std::vector<vec3> &x = get_vec_data(*__M, 0);
    for (int i = 0; i < x.size(); i++) {
      if (x[i].hasNaN()) {
        std::cout << "nan" << std::endl;
        std::cout << i << std::endl;
        std::cout << x[i] << std::endl;
        exit(0);
      }
    }
  }

  void update_positions(real dt, const std::vector<vec3> &dx) {

    std::vector<vec3> &x = get_vec_data(*__M, 0);
    std::vector<vec3> &_dx = get_vec_data(*__M, __vdatum_id);

    for (int i = 0; i < _dx.size(); i++) {
      _dx[i] = dx[i];
      x[i] += dt * dx[i];
    }
  }

  void step(bool merge_edges_ = true, bool break_cycles_ = true) {
    for (int k = 0; k < 1; k++) {
      if (!has_active_edges()) {
        std::cout << "dynamic::step: mesh has no active edges, skipping"
                  << std::endl;
        return;
      }

      std::cout << "subd" << std::endl;
      subdivide_edges();
      delete_degenerates(*__M);

      if (!has_active_edges()) {
        std::cout << "dynamic::step: collapsed to nothing after subd"
                  << std::endl;
        return;
      }

      if (break_cycles_) {
        std::cout << "break, ";
        break_cycles();
        std::cout << "degenerates ";
        delete_degenerates(*__M);
      }
      std::cout << std::endl;

      if (!has_active_edges()) {
        std::cout << "dynamic::step: collapsed to nothing after break"
                  << std::endl;
        return;
      }

      std::cout << "collapse" << std::endl;
      collapse_edges();
      delete_degenerates(*__M);

      if (!has_active_edges()) {
        std::cout << "dynamic::step: collapsed to nothing after collapse"
                  << std::endl;
        return;
      }

      if (merge_edges_) {
        std::cout << "merge" << std::endl;
        merge_edges();
        delete_degenerates(*__M);
      }

      if (!has_active_edges()) {
        std::cout << "dynamic::step: collapsed to nothing after merge"
                  << std::endl;
        return;
      }

      std::cout << "flip" << std::endl;
      flip_edges();
      delete_degenerates(*__M);
    }

    if (has_active_edges())
      pack(*__M);
  }

  void step(real dt, const std::vector<vec3> &dx) {
    update_positions(dt, dx);
    test_nan();
    step();
  }

  void set_flip_pred(OpPredicateFcn f) { _flip_pred = f; }
  void set_merge_pred(MergePredicateFcn f) { _merge_pred = f; }
  void set_collapse_pred(OpPredicateFcn f) { _collapse_pred = f; }

  shell::ptr __M;
  index_t __vdatum_id;
  real _Cc, _Cs, _Cm; // collapse, stretch, bridge
  real _max_flip_vol_ratio = 0.5;
  real _max_flip_dihedral = 1.0;

  OpPredicateFcn _flip_pred;
  MergePredicateFcn _merge_pred;
  OpPredicateFcn _collapse_pred;

  arp::bvh_tree<2>::ptr edge_tree_;
  arp::bvh_tree<3>::ptr face_tree_;
};

} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif