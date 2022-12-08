#include "m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "manifold/hepworth/constraints.hpp"
#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "m2_refactor.hpp"
#include "primitive_operations.hpp"

#include <functional>
#include <iostream>
#include <vector>
#include <zlib.h>

#ifndef __ASAWA_DYNAMIC_SURFACE__
#define __ASAWA_DYNAMIC_SURFACE__
namespace asawa {
typedef double real;
typedef int index_t;

using vec3 = Eigen::Matrix<real, 3, 1>;
using vec4 = Eigen::Matrix<real, 4, 1>;

using real_datum = datum_t<real>;
using vec3_datum = datum_t<vec3>;

void op_edges(manifold &M, //
              std::vector<index_t> &edges_to_op,
              std::function<index_t(index_t cid, manifold &m)> func) {

  for (auto d : M.get_data()) {
    if (d->type() != VERTEX)
      continue;
    d->alloc(edges_to_op.size());
  }

  for (index_t i = 0; i < edges_to_op.size(); i++) {
    for (auto d : M.get_data()) {
      if (d->type() != VERTEX)
        continue;
      if (M.next(edges_to_op[i]) < 0)
        continue;
      d->calc(i, M, edges_to_op[i]);
    }
  }

  std::vector<index_t> post_edges(edges_to_op.size(), -1);

  //#pragma omp parallel for shared(post_edges)
  // before parallelize, need to preallocate;
  for (index_t i = 0; i < edges_to_op.size(); i++) {
    index_t ic0 = edges_to_op[i];
    //#pragma omp critical
    if (M.next(ic0) < 0)
      continue;
    post_edges[i] = func(ic0, M);
  }

  // verts.resize(M.vert_count());
  for (auto d : M.get_data()) {
    if (d->type() != VERTEX)
      continue;
    d->resize(M.vert_count());
  }

  for (index_t i = 0; i < post_edges.size(); i++) {
    for (auto d : M.get_data()) {
      if (d->type() != VERTEX)
        continue;

      if (M.next(post_edges[i]) < 0)
        continue;
      d->map(M.vert(post_edges[i]), i);
    }
  }
}

void subdivide_edges(manifold &M) {
  std::vector<index_t> edges_to_divide;
  edges_to_divide.push_back(3);
  edges_to_divide.push_back(5);
  edges_to_divide.push_back(8);
  edges_to_divide.push_back(13);
  edges_to_divide.push_back(21);
  op_edges(M, edges_to_divide,
           [](index_t cid, manifold &m) { return subdivide_edge(m, cid); });
}

void collapse_edges(manifold &M) {
  std::vector<index_t> edges_to_divide;
  edges_to_divide.push_back(4);
  // edges_to_divide.push_back(6);
  // edges_to_divide.push_back(9);
  //  edges_to_divide.push_back(14);
  edges_to_divide.push_back(22);
  op_edges(M, edges_to_divide,
           [](index_t cid, manifold &m) { return collapse_edge(m, cid); });
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
    sort(edges.begin(), edges.end(),
         [this](const index_t &ca, const index_t &cb) -> bool {
           real dva = length(ca, _M.other(ca), _M, _data);
           real dvb = length(cb, _M.other(cb), _M, _data);
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

  void flip_edges() {

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

      index_t c2 = __M->other(c0);
      index_t c3 = __M->prev(c2);

      vec3_datum::ptr coord_datum =
          static_pointer_cast<vec3_datum>(__M->get_datum(0));
      const std::vector<vec3> &data = coord_datum->data();

      vec3 v0 = data[__M->vert(c0)];
      vec3 v1 = data[__M->vert(c1)];
      vec3 v2 = data[__M->vert(c2)];
      vec3 v3 = data[__M->vert(c3)];
      /*
      real A0 = va::calculate_area<real>(v0, v1, v2);
      real A1 = va::calculate_area<real>(v2, v3, v0);

      if (A0 / A1 > 4.0) {
        flip_edge(*__M, c0);
        continue;
      } else if (A1 / A0 > 4) {
        flip_edge(*__M, c0);
        continue;
      }
      */
      // std::cout << A0 << "  " << A1 << " " << A0 / A1 << " " << A1 / A0
      //           << std::endl;
      //  interior angles
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
      real cosN0 = va::dot(N00, N01);
      real tSame = M_PI - acos(cosN0);
      real cosN1 = va::dot(N10, N11);
      real tFlip = M_PI - acos(cosN1);
      real nFlip = tFlip;

      real eFlip = cFlip * cFlip + 0.5 * tFlip * tFlip;
      real eSame = cSame * cSame + 0.0 * tSame * tSame;

      // if (false) {
      if (eFlip < 0.5 * eSame) {
        flip_edge(*__M, c0);
      }
    }
  }

  void split_edges() {
    using comp_great = manifold_data_comp<vec3, std::greater<real>>;
    vec3_datum::ptr coord_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    const std::vector<vec3> &coords = coord_datum->data();
    auto cmp = comp_great(coords, _Cs, *__M);

    std::vector<index_t> edges_to_divide = gather_edges<comp_great>(*__M, cmp);

    op_edges(*__M, edges_to_divide,
             [](index_t cid, manifold &m) { return subdivide_edge(m, cid); });
  }

  void collapse_edges() {
    using comp_less = manifold_data_comp<vec3, std::less<real>>;
    vec3_datum::ptr coord_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    const std::vector<vec3> &coords = coord_datum->data();
    auto cmp = comp_less(coords, _Cc, *__M);

    std::vector<index_t> edges_to_divide = gather_edges<comp_less>(*__M, cmp);

    op_edges(*__M, edges_to_divide,
             [](index_t cid, manifold &m) { return collapse_edge(m, cid); });
  }

  void step(const std::vector<vec3> &dx) {

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

    split_edges();
    collapse_edges();
    flip_edges();
  }

  manifold::ptr __M;
  index_t __vdatum_id;
  real _Cc, _Cs, _Cm; // collapse, stretch, bridge
};

} // namespace asawa
#endif