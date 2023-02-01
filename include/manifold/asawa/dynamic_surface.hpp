
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "manifold/arp/arp.h"

#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "datums.hpp"

#include "datum_x.hpp"

#include "m2_refactor.hpp"
#include "manifold/bontecou/laplacian_refactor.hpp"

#include "primitive_operations.hpp"
#include "subdivide.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>
#include <zlib.h>

#ifndef __ASAWA_DYNAMIC_SURFACE__
#define __ASAWA_DYNAMIC_SURFACE__
namespace asawa {
typedef double real;
typedef int index_t;
using vec4 = Eigen::Matrix<real, 4, 1>;
using corner1 = std::array<index_t, 1>;
using corner2 = std::array<index_t, 2>;
using corner4 = std::array<index_t, 4>;

real dist_line_line(manifold &M, index_t cA0, index_t cB0,
                    const std::vector<vec3> &x) {
  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);
  index_t vA0 = M.vert(cA0);
  index_t vA1 = M.vert(cA1);
  index_t vB0 = M.vert(cB0);
  index_t vB1 = M.vert(cB1);
  vec3 xA0 = x[vA0];
  vec3 xA1 = x[vA1];
  vec3 xB0 = x[vB0];
  vec3 xB1 = x[vB1];

  real d0 = 1.0 / 2.0 * ((xB0 - xA0).norm() + (xB1 - xA1).norm());
  return d0;
};

real dist_line_line_cen(manifold &M, index_t cA0, index_t cB0,
                        const std::vector<vec3> &x) {
  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);
  index_t vA0 = M.vert(cA0);
  index_t vA1 = M.vert(cA1);
  index_t vB0 = M.vert(cB0);
  index_t vB1 = M.vert(cB1);
  vec3 xA0 = x[vA0];
  vec3 xA1 = x[vA1];
  vec3 xB0 = x[vB0];
  vec3 xB1 = x[vB1];

  real d0 = (0.5 * (xA1 + xA0) - 0.5 * (xB1 + xA0)).norm();
  return d0;
};

// this is ugly, but we have to do it this way, with two lists because
// the callback on the data
real line_line_min(const index_t &idT, //
                   const std::vector<index_t> &t_inds,
                   const vector<vec3> &t_x, //
                   const index_t &idS,      //
                   const std::vector<index_t> &s_inds,
                   const vector<vec3> &s_x) {
  index_t vT0 = t_inds[2 * idT + 0];
  index_t vT1 = t_inds[2 * idT + 1];
  index_t vS0 = s_inds[2 * idS + 0];
  index_t vS1 = s_inds[2 * idS + 1];
  if (idT == idS)
    return std::numeric_limits<real>::infinity();

  if (vT0 == vS0)
    return std::numeric_limits<real>::infinity();
  if (vT1 == vS1)
    return std::numeric_limits<real>::infinity();
  if (vT0 == vS1)
    return std::numeric_limits<real>::infinity();
  if (vT1 == vS0)
    return std::numeric_limits<real>::infinity();

  const vec3 &xA0 = t_x[t_inds[2 * idT + 0]];
  const vec3 &xA1 = t_x[t_inds[2 * idT + 1]];
  const vec3 &xB0 = s_x[s_inds[2 * idS + 0]];
  const vec3 &xB1 = s_x[s_inds[2 * idS + 1]];
#if 0
  real d0 = 1.0 / 2.0 * ((xB0 - xA0).norm() + (xB1 - xA1).norm());
  real d1 = 1.0 / 2.0 * ((xB0 - xA1).norm() + (xB1 - xA0).norm());
  return min(d0, d1);
#else
  real d0 = (0.5 * (xA1 + xA0) - 0.5 * (xB1 + xA0)).norm();
  return d0;
#endif
};

std::mt19937_64 rng;
std::uniform_real_distribution<real> unif(0.0, 1.0);

void debug_line(manifold &M,        //
                const index_t &cA0, //
                const vector<vec3> &x) {

  index_t cA1 = M.other(cA0);
  index_t vA0 = M.vert(cA0);
  index_t vA1 = M.vert(cA1);

  const vec3 &ca0 = x[vA0];
  const vec3 &ca1 = x[vA1];

  vec4 c(unif(rng), unif(rng), unif(rng), 1.0);
  gg::geometry_logger::line(ca0, ca1, c);
};

void debug_line_line(manifold &M,        //
                     const index_t &cA0, //
                     const index_t &cB0, //
                     const vector<vec3> &x) {

  index_t cA1 = M.other(cA0);
  index_t cB1 = M.other(cB0);
  index_t vA0 = M.vert(cA0);
  index_t vA1 = M.vert(cA1);
  index_t vB0 = M.vert(cB0);
  index_t vB1 = M.vert(cB1);

  const vec3 &ca0 = x[vA0];
  const vec3 &ca1 = x[vA1];
  const vec3 &cb0 = x[vB0];
  const vec3 &cb1 = x[vB1];

  vec4 c(unif(rng), unif(rng), unif(rng), 1.0);
  gg::geometry_logger::line(ca0, ca1, c);
  gg::geometry_logger::line(cb0, cb1, c);

  gg::geometry_logger::line(0.5 * (ca0 + ca1), 0.5 * (cb0 + cb1), c);
};

void debug_edge_normal(manifold &M,       //
                       const index_t &c0, //
                       const vector<vec3> &x) {

  vec3 N = edge_normal(M, c0, x);
  vec3 cen = edge_center(M, c0, x);
  vec4 col(unif(rng), unif(rng), unif(rng), 1.0);
  gg::geometry_logger::line(cen, cen + 0.1 * N, col);
};

template <int OP, int C_ALLOC, int V_ALLOC, int F_ALLOC, int MSIZE>
void op_edges(manifold &M, //
              std::vector<index_t> &edges_to_op,
              std::function<std::array<index_t, MSIZE>(
                  index_t i,                         //
                  index_t cs,                        //
                  index_t vs,                        //
                  index_t fs,                        //
                  const std::vector<index_t> &edges, //
                  manifold &m)>
                  func) {
  int STRIDE = OP < 2 ? 1 : 2;
  size_t cstart = M.corner_count();
  size_t vstart = M.vert_count();
  size_t fstart = M.face_count();
  size_t Ne = edges_to_op.size();
  M.inflate_edge_pairs(C_ALLOC * Ne);
  M.inflate_verts(V_ALLOC * Ne);
  M.inflate_faces(F_ALLOC * Ne);

  for (auto d : M.get_data()) {
    if (d->type() != VERTEX)
      continue;
    d->alloc(V_ALLOC * edges_to_op.size());
  }

  for (index_t i = 0; i < edges_to_op.size(); i += STRIDE) {
    for (auto d : M.get_data()) {
      if (d->type() != VERTEX)
        continue;
      if (M.next(edges_to_op[i]) < 0)
        continue;

      if (OP == 0) {
        // subd
        index_t c0 = edges_to_op[i];
        index_t c1 = M.other(c0);
        d->calc(vstart + i, M, M.vert(c0), M.vert(c1));
      } else if (OP == 1) {
        // collapse
        index_t c0 = edges_to_op[i];
        index_t c1 = M.other(c0);
        d->calc(M.vert(c1), M, M.vert(c0), M.vert(c1));
      } else if (OP == 2) {
        // merge
        index_t cA0 = edges_to_op[i + 0];
        index_t cA1 = M.other(cA0);
        index_t cB0 = edges_to_op[i + 1];
        index_t cB1 = M.other(cB0);
        d->calc(M.vert(cA0), M, M.vert(cA0), M.vert(cB0));
        d->calc(M.vert(cA1), M, M.vert(cA1), M.vert(cB1));
        // d->calc(M.vert(cB0), M, M.vert(cA0), M.vert(cB0));
        // d->calc(M.vert(cB1), M, M.vert(cA1), M.vert(cB1));

        d->calc(vstart + 2 * i + 0, M, M.vert(cA0), M.vert(cB0));
        d->calc(vstart + 2 * i + 1, M, M.vert(cA1), M.vert(cB1));
      }
    }
  }

  //#pragma omp parallel for shared(post_edges)
  // before parallelize, need to preallocate;
  for (index_t i = 0; i < edges_to_op.size(); i += STRIDE) {
    index_t ic0 = edges_to_op[i];
    //#pragma omp critical
    if (M.next(ic0) < 0)
      continue;
    auto p = func(i, cstart, vstart, fstart, edges_to_op, M);
  }
}

// template <int STRIDE, int C_ALLOC, int V_ALLOC, int F_ALLOC, int MSIZE>
auto subdivide_op = op_edges<0, 3, 1, 2, 1>;
auto collapse_op = op_edges<1, 0, 0, 0, 1>;
auto merge_op = op_edges<2, 0, 2, 0, 4>;

void subdivide_edges(manifold &M) {
  std::vector<index_t> edges_to_divide;
  edges_to_divide.push_back(3);
  edges_to_divide.push_back(5);
  edges_to_divide.push_back(8);
  edges_to_divide.push_back(13);
  edges_to_divide.push_back(21);
  subdivide_op(M, edges_to_divide,
               [](index_t i,  //
                  index_t cs, //
                  index_t vs, //
                  index_t fs, //
                  const std::vector<index_t> &edges, manifold &m) -> corner1 {
                 return {subdivide_edge(m, edges[i])};
               });
}

void collapse_edges(manifold &M) {
  std::vector<index_t> edges_to_divide;
  edges_to_divide.push_back(4);
  // edges_to_divide.push_back(6);
  // edges_to_divide.push_back(9);
  //  edges_to_divide.push_back(14);
  edges_to_divide.push_back(22);
  collapse_op(M, edges_to_divide,
              [](index_t i,  //
                 index_t cs, //
                 index_t vs, //
                 index_t fs, //
                 const std::vector<index_t> &edges, manifold &m) -> corner1 {
                return {collapse_edge(m, edges[i])};
              });
}

real length(index_t c0, index_t c1, const manifold &M,
            const std::vector<vec3> &data) {
  vec3 v0 = data[M.vert(c0)];
  vec3 v1 = data[M.vert(c1)];
  return (v1 - v0).norm();
}

template <typename T, typename comp> class manifold_data_comp {
public:
  manifold_data_comp(const std::vector<T> &data, real eps, const manifold &M)
      : _M(M), _data(data), _eps(eps) {}
  bool operator()(index_t c0, index_t c1) const {
    real dv = length(c0, c1, _M, _data);
    // std::cout << "comp: " << dv << " " << _eps << std::endl;
    return comp{}(dv, _eps);
  }

  std::vector<index_t> get_edges() const {

    std::vector<index_t> edges = _M.get_edge_range();
    std::vector<real> lengths(_M.corner_count() / 2, -1);

    sort(edges.begin(), edges.end(),
         [this, &lengths](const index_t &ca, const index_t &cb) -> bool {
           real dva = lengths[ca / 2];
           real dvb = lengths[cb / 2];
           if (dva < 0) {
             dva = length(ca, _M.other(ca), _M, _data);
             lengths[ca / 2] = dva;
           }
           if (dvb < 0) {
             dvb = length(cb, _M.other(cb), _M, _data);
             lengths[cb / 2] = dvb;
           }

           return comp{}(dva, dvb);
         });

    return edges;
  }

  const std::vector<T> &_data;
  const manifold &_M;
  real _eps;
};

template <typename comparator>
std::vector<index_t> gather_edges(manifold &M, const comparator &comp) {
  std::vector<index_t> edges = comp.get_edges();
  std::vector<bool> face_flags(M.face_count(), false);
  std::vector<index_t> edges_out;
  edges_out.reserve(edges.size());

  for (int i = 0; i < edges.size(); i++) {
    index_t c0 = edges[i];
    index_t c1 = M.other(c0);

    index_t f0 = M.face(c0);
    index_t f1 = M.face(c1);

    if (f0 < 0 || f1 < 0)
      continue;
    if (face_flags[f0])
      continue;
    if (face_flags[f1])
      continue;

    if (comp(c0, c1)) {
      edges_out.push_back(c0);
      face_flags[f0] = true;
      face_flags[f1] = true;
    }
  } // namespace asawa

  return edges_out;
}
#if 0 
void smoothMesh(manifold &M, real C, int N) {

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

class dynamic_surface {
public:
  typedef std::shared_ptr<dynamic_surface> ptr;

  static ptr create(manifold::ptr M, real Cc, real Cs, real Cm) {
    return std::make_shared<dynamic_surface>(M, Cc, Cs, Cm);
  }

  dynamic_surface(manifold::ptr M, real Cc, real Cs, real Cm) : __M(M) {
    _Cc = Cc;
    _Cs = Cs;
    _Cm = Cm;

    std::vector<vec3> velocities(__M->vert_count(), vec3::Zero());
    datum_t<vec3>::ptr vdata =
        datum_t<vec3>::create(prim_type::VERTEX, velocities);
    __vdatum_id = __M->insert_datum(vdata);
  };

  void delete_degenerates(manifold &M) {
    vec3_datum::ptr x_datum = static_pointer_cast<vec3_datum>(M.get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    for (int i = 0; i < M.corner_count(); i++) {
      if (M.next(i) < 0)
        continue;
      index_t c0 = i;
      index_t c1 = M.other(i);
      if (M.vert(c0) != M.vert(c1))
        continue;

      collapse_edge(M, c0, true);
    }

    for (int i = 0; i < M.face_count(); i++) {
      if (M.fbegin(i) < 0)
        continue;

      real a0 = face_area(M, i, x);
      if (a0 < 1e-10) {
        // std::cout << "deg face" << std::endl;
        //         M.fprintv(i);
        //  M.cprint(M.fbegin(i));
      }

      if (M.fsize(i) > 2)
        continue;
      index_t c0 = M.fbegin(i);
      index_t c1 = M.other(c0);

      merge_face(M, c0, c1);
    }

    for (int i = 0; i < M.vert_count(); i++) {
      if (M.vbegin(i) < 0)
        continue;

      if (M.vsize(i) > 16) {
        vec3 N = vert_normal(M, i, x);
        vec4 cola(0.0, 1.0, 1.0, 0.0);
        gg::geometry_logger::line(x[i], x[i] + 0.1 * N, cola);
        M.for_each_vertex(i, [&x, &i, cola](index_t cid, manifold &m) {
          index_t j = m.vert(m.next(cid));
          gg::geometry_logger::line(x[i], x[j], cola);
        });
      }

      if (M.vsize(i) > 3)
        continue;

      vec3 N = vert_normal(M, i, x);
      vec4 cola(0.2, 0.5, 1.0, 0.0);
      gg::geometry_logger::line(x[i], x[i] + 0.1 * N, cola);
      remove_vertex(M, i);
    }
  }

  index_t align_edges(manifold &M, index_t cA0, index_t cB0,
                      const std::vector<vec3> &x) {

    index_t cB1 = M.other(cB0);
    real d0 = dist_line_line(M, cA0, cB0, x);
    real d1 = dist_line_line(M, cA0, cB1, x);

    if (d1 < d0)
      return cB1;

    return cB0;
  }

  void trim_collected(manifold &M, const std::vector<vec3> &x,
                      std::vector<std::array<index_t, 2>> &collected) {
    std::vector<bool> flags(M.corner_count() / 2, false);
    real tol = _Cm;
    collected.erase(std::remove_if(collected.begin(), collected.end(),
                                   [tol, &M, &x, &flags](const auto &p) {
                                     if (p[0] < 0)
                                       return true;
                                     if (p[1] < 0)
                                       return true;

                                     vec3 cenA = edge_center(M, p[0], x);
                                     vec3 cenB = edge_center(M, p[1], x);
                                     real dist = (cenA - cenB).norm();

                                     if (dist > tol) {
                                       return true;
                                     }

                                     vec3 NA = edge_normal(M, p[0], x);
                                     vec3 NB = edge_normal(M, p[1], x);
                                     real angle = va::dot(NA, NB);

                                     if (angle > -0.25) {
                                       return true;
                                     }

                                     if (flags[p[0] / 2])
                                       return true;
                                     if (flags[p[1] / 2])
                                       return true;

                                     flags[p[0] / 2] = true;
                                     flags[p[1] / 2] = true;

                                     return false;
                                   }),
                    collected.end());

    std::sort(collected.begin(), collected.end(),
              [&M, &x](const auto &pA, const auto &pB) {
                real dA = dist_line_line(M, pA[0], pA[1], x);
                real dB = dist_line_line(M, pB[0], pB[1], x);
                return dA < dB;
              });

    auto it = std::unique(collected.begin(), collected.end());

    collected.erase(it, collected.end());
  }

  vector<std::array<index_t, 2>> get_merge_pairs(manifold &M,
                                                 std::vector<vec3> &x) {

    std::vector<index_t> edge_verts = __M->get_edge_vert_ids();
    std::vector<index_t> edge_map = __M->get_edge_map();

    edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 64);
    // edge_tree->debug(edges, x);
    real tol = this->_Cm * this->_Cm;

    std::vector<std::array<index_t, 2>> collected(edge_verts.size() / 2);
    //#pragma omp parallel for
    for (int i = 0; i < edge_verts.size(); i += 2) {
      index_t e0 = i / 2;
      index_t e1 = arp::getNearest<2, 2>(e0, edge_verts, x, //
                                         *edge_tree,        //
                                         tol, &line_line_min);

      index_t c0 = edge_map[e0];
      // debug_line(M, c0, x);
      index_t c1 = -1;
      if (e1 > 0) {
        // std::cout << ii << " " << nearest << std::endl;
        c1 = edge_map[e1];
        c1 = align_edges(M, c0, c1, x);
        // debug_line_line(M, c0, c1, x);
        // debug_edge_normal(M, c0, x);
        // debug_edge_normal(M, c1, x);
      }

      if (c0 > c1)
        std::swap(c0, c1);

      collected[i / 2] = {c0, c1};
    }

    trim_collected(M, x, collected);
    return collected;
  }

  void merge_edges() {
    // edge e = c / 2;
    manifold &M = *__M;

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    auto collected = get_merge_pairs(M, x);

    std::vector<index_t> f_collect(2 * collected.size());
    for (int i = 0; i < collected.size(); i++) {
      f_collect[2 * i + 0] = collected[i][0];
      f_collect[2 * i + 1] = collected[i][1];
    }
    real tol = this->_Cm;
    merge_op(*__M, f_collect,
             [&x, tol](index_t i,                         //
                       index_t cs,                        //
                       index_t vs,                        //
                       index_t fs,                        //
                       const std::vector<index_t> &edges, //
                       manifold &M) -> corner4 {
               index_t c0A = edges[i + 0];
               index_t c0B = edges[i + 1];

               if (M.next(c0A) < 0 || M.next(c0B) < 0) {
                 return {-1, -1, -1, -1};
               }

               // debug_line_line(M, c0A, c0B, x);
               // debug_edge_normal(M, c0A, x);
               // debug_edge_normal(M, c0B, x);

               return merge_edge(M, c0A, c0B, vs + 2 * i + 0, vs + 2 * i + 1);
             });
  }

  void break_cycles() {
    // edge e = c / 2;
    manifold &M = *__M;

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    auto edges = M.get_edge_range();
    std::vector<std::array<index_t, 2>> collected;

    for (int i = 0; i < edges.size(); i++) {
      index_t c0 = edges[i + 0];
      index_t c1 = M.other(c0);
      index_t v1 = M.vert(c1);
      if (count_cycle(M, c0) < 2)
        continue;

      std::array<index_t, 2> pair = {0, 0};
      int j = 0;
      M.for_each_vertex(M.vert(c0), [v1, &pair, &j](index_t ci, manifold &M) {
        index_t vi = M.vert(M.next(ci));
        if (vi == v1 && j < 2) {
          pair[j++] = ci;
        }
      });
      collected.push_back(pair);
    }

    std::vector<index_t> f_collect(2 * collected.size());
    for (int i = 0; i < collected.size(); i++) {
      f_collect[2 * i + 0] = collected[i][0];
      f_collect[2 * i + 1] = collected[i][1];
    }

    merge_op(*__M, f_collect,
             [&x](index_t i,                         //
                  index_t cs,                        //
                  index_t vs,                        //
                  index_t fs,                        //
                  const std::vector<index_t> &edges, //
                  manifold &M) -> corner4 {
               index_t c0A = edges[i + 0];
               index_t c0B = edges[i + 1];

               debug_line_line(M, c0A, c0B, x);
               debug_edge_normal(M, c0A, x);
               debug_edge_normal(M, c0B, x);
               // std::cout << "break_cycle" << std::endl;
               // M.cprint(c0A);
               return merge_edge(M, c0A, c0B, vs + 2 * i + 0, vs + 2 * i + 1);
             });
  }

  bool skip_flip(manifold &M, index_t corner) {
    index_t c0 = corner;
    index_t c1 = M.other(c0);
    index_t v0 = M.vert(c0);
    index_t v1 = M.vert(c1);

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
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::vector<index_t> edges = __M->get_edge_range();
    for (int i = 0; i < edges.size(); i++) {
      int card = rand() % edges.size();
      index_t et = edges[i];
      edges[i] = edges[card];
      edges[card] = et;
    }

    for (int i = 0; i < edges.size(); i++) {
      // std::cout << "A" << std::endl;

      index_t c0 = edges[i];
      index_t c1 = __M->prev(c0);

      if (skip_flip(*__M, c0))
        continue;

      index_t c2 = __M->other(c0);
      index_t c3 = __M->prev(c2);

      vec3_datum::ptr coord_datum =
          static_pointer_cast<vec3_datum>(__M->get_datum(0));
      const std::vector<vec3> &data = coord_datum->data();
      vec3 v0 = data[__M->vert(c0)];
      vec3 v1 = data[__M->vert(c1)];
      vec3 v2 = data[__M->vert(c2)];
      vec3 v3 = data[__M->vert(c3)];

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
      real dt = tFlip - tSame;
      // std::cout << tFlip << " " << tSame << " " << dt << std::endl;
      real eFlip = cFlip * cFlip + 150.0 * dt * dt;
      real eSame = cSame * cSame;
      // real eFlip = cFlip * cFlip + tFlip * tFlip;
      // real eSame = cSame * cSame + tSame * tSame;

      // std::cout << eFlip << " " << eSame << std::endl;
      //  if (false) {
      if (eFlip < 0.75 * eSame) {
        flip_edge(*__M, c0);
      }
    }
  }

  void subdivide_edges() {
    using comp_great = manifold_data_comp<vec3, std::greater<real>>;
    vec3_datum::ptr coord_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    const std::vector<vec3> &coords = coord_datum->data();
    auto cmp = comp_great(coords, _Cs, *__M);

    std::vector<index_t> edges_to_divide = gather_edges<comp_great>(*__M, cmp);

    subdivide_op(*__M, edges_to_divide,
                 [](index_t i,  //
                    index_t cs, //
                    index_t vs, //
                    index_t fs, //
                    const std::vector<index_t> &edges, manifold &m) -> corner1 {
                   return {subdivide_edge(m, edges[i],    //
                                          vs + i,         //
                                          cs + 6 * i + 0, //
                                          cs + 6 * i + 2, //
                                          cs + 6 * i + 4, //
                                          fs + 2 * i + 0, //
                                          fs + 2 * i + 1  //
                                          )};
                 });
  }

  void collapse_edges() {
    using comp_less = manifold_data_comp<vec3, std::less<real>>;
    vec3_datum::ptr coord_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    const std::vector<vec3> &coords = coord_datum->data();
    auto cmp = comp_less(coords, _Cc, *__M);

    std::vector<index_t> edges_to_divide = gather_edges<comp_less>(*__M, cmp);

    collapse_op(*__M, edges_to_divide,
                [](index_t i,  //
                   index_t cs, //
                   index_t vs, //
                   index_t fs, //
                   const std::vector<index_t> &edges, manifold &m) -> corner1 {
                  return {collapse_edge(m, edges[i])};
                });
  }

  void update_positions(const std::vector<vec3> &dx) {

    vec3_datum::ptr c_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &coords = c_datum->data();

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(__vdatum_id));
    std::vector<vec3> &_dx = v_datum->data();

    for (int i = 0; i < _dx.size(); i++) {
      _dx[i] = dx[i];
      coords[i] += dx[i];
    }
  }

  void step(const std::vector<vec3> &dx) {
    update_positions(dx);

    subdivide_edges();
    delete_degenerates(*__M);

    break_cycles();
    delete_degenerates(*__M);

    collapse_edges();
    delete_degenerates(*__M);

    merge_edges();
    delete_degenerates(*__M);

    flip_edges();
    delete_degenerates(*__M);

    pack(*__M);
  }

  manifold::ptr __M;
  index_t __vdatum_id;
  real _Cc, _Cs, _Cm; // collapse, stretch, bridge

  arp::aabb_tree<2>::ptr edge_tree;
  arp::aabb_tree<3>::ptr face_tree;
};

} // namespace asawa
#endif