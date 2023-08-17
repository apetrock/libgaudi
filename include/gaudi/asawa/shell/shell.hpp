#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

#ifndef __ASAWA_SHELL__
#define __ASAWA_SHELL__

namespace gaudi {
namespace asawa {

class datum;

namespace shell {

typedef int index_t;

typedef std::shared_ptr<datum> datum_ptr;

class shell {
public:
  typedef std::shared_ptr<shell> ptr;

  static ptr create(const std::vector<index_t> &corners_next,
                    const std::vector<index_t> &corners_vert,
                    const std::vector<index_t> &corners_face) {

    return std::make_shared<shell>(corners_next, corners_vert, corners_face);
  }

  shell(const std::vector<index_t> &corners_next,
        const std::vector<index_t> &corners_vert,
        const std::vector<index_t> &corners_face)
      : __corners_next(corners_next), __corners_vert(corners_vert),
        __corners_face(corners_face) {
    this->update_head();
    this->update_prev();
  }

  void update_head() {
    auto update = [](std::vector<index_t> &p_beg,
                     const std::vector<index_t> &c_ptr) {
      int size = *std::max_element(c_ptr.begin(), c_ptr.end());

      p_beg = std::vector<index_t>(size + 1, -1);

      for (int i = 0; i < c_ptr.size(); i++) {
        // std::cout << i << " " << c_ptr[i] << " " << std::endl;
        p_beg[c_ptr[i]] = i;
      }
    };

    update(__face_begin, __corners_face);
    update(__vert_begin, __corners_vert);
  }

  void update_prev() {
    // can't use the iterator because the iterator needs prev.
    __corners_prev = std::vector<index_t>(__corners_next.size(), -1);
    for (int i = 0; i < __face_begin.size(); i++) {
      int j_end = __face_begin[i];
      int j0 = __corners_next[j_end];
      bool it = true;
      while (it) {
        it = j0 != j_end;
        __corners_prev[__corners_next[j0]] = j0;
        j0 = __corners_next[j0];
      }
    }
  }

  index_t insert_datum(datum_ptr datum) {
    __data.push_back(datum);
    return __data.size() - 1;
  }
  /*incomplete data, bummer*/
  datum_ptr &get_datum(index_t i) { return __data[i]; }
  const datum_ptr &const_get_datum(index_t i) const { return __data[i]; }

  std::vector<datum_ptr> &get_data() { return __data; }

  size_t vert_count() const { return __vert_begin.size(); }
  size_t face_count() const { return __face_begin.size(); }
  size_t corner_count() const { return __corners_next.size(); }

  index_t other(index_t id) const { return 2 * (id / 2) + (id + 1) % 2; };
  bool edge_equal(index_t c0, index_t c1) { return c0 / 2 == c1 / 2; }
  index_t vprev(index_t id) const { return this->next(this->other(id)); }
  index_t vnext(index_t id) const { return this->other(this->prev(id)); }

  index_t next(index_t id) const { return __corners_next[id]; }
  index_t prev(index_t id) const { return __corners_prev[id]; }
  index_t vert(index_t id) const { return __corners_vert[id]; }
  index_t face(index_t id) const { return __corners_face[id]; }

  index_t fbegin(index_t id) const { return __face_begin[id]; }
  index_t fend(index_t id) const { return prev(__face_begin[id]); }

  index_t vbegin(index_t id) const { return __vert_begin[id]; }

  index_t vend(index_t id) const { return vprev(__vert_begin[id]); }

  void set_vbegin(index_t id, index_t c) { __vert_begin[id] = c; }
  void set_fbegin(index_t id, index_t c) { __face_begin[id] = c; }

  void set_next(index_t id, index_t c) { __corners_next[id] = c; }
  void set_prev(index_t id, index_t c) { __corners_prev[id] = c; }

  void set_vert(index_t cid, index_t v) {
    __corners_vert[cid] = v;
    set_vbegin(v, cid);
  }

  void set_face(index_t cid, index_t f) {
    __corners_face[cid] = f;
    set_fbegin(f, cid);
  }

  index_t find_edge_from_verts(index_t v0, index_t v1) const {
    index_t c = -1;
    const_for_each_vertex(v0, [&](int ci, const shell &M) {
      index_t cn = M.other(ci);
      if (M.vert(cn) == v1)
        c = ci;
    });
    return c;
  }

  void link(index_t c0, index_t c1) {
    // std::cout << __PRETTY_FUNCTION__ << c0 << " " << c1 << std::endl;
    /*
    if (vert(c0) == vert(c1)) {
      print_stacktrace();
      __builtin_frame_address(1);
    }
    assert(vert(c0) != vert(c1));
    */
    set_next(c0, c1);
    set_prev(c1, c0);
  }

  void uber_assert() {
    for (int i = 0; i < corner_count(); i++) {
      index_t c0 = i;
      index_t c1 = other(c0);
      if (next(c0) < 0)
        continue;

      if (vert(c0) == vert(c1)) {
        std::cout << "    " << __PRETTY_FUNCTION__ << " c: " << c0 << " " << c1
                  << " v: " << vert(c0) << " " << vert(c1) //
                  << " vs: " << vsize(vert(c0)) << " " << vsize(vert(c1))
                  << std::endl;
      }
      // std::cout << c0 << " " << next(c0) << " " << vert(c0) << std::endl;
      if (vsize(vert(c0)) > 1)
        assert(vert(c0) != vert(c1));
    }

    for (int i = 0; i < vert_count(); i++) {
      if (vbegin(i) < 0)
        continue;

      if (next(vbegin(i)) < 0) {
        std::cout << "-" << __PRETTY_FUNCTION__ << " v: " << i << " "
                  << vbegin(i) << " " << std::endl
                  << std::flush;
        vprintv(i);
      }
      assert(next(vbegin(i)) >= 0);

      if (vert(vbegin(i)) != i) {
        std::cout << "-" << __PRETTY_FUNCTION__ << " v: " //
                  << i << " "                             //
                  << vbegin(i) << " "                     //
                  << vert(vbegin(i)) << " "               //
                  << std::endl
                  << std::flush;

        for_each_vertex(
            i, [](index_t cid, shell &m) { std::cout << m.vert(cid) << " "; });
        std::cout << std::endl;
      }
      assert(vert(vbegin(i)) == i);
    }

    for (int i = 0; i < face_count(); i++) {
      if (__face_begin[i] < 0)
        continue;
      if (fsize(i) != 3) {
        std::cout << "-" << __PRETTY_FUNCTION__ << " f: " << i << " "
                  << __face_begin[i] << " " << fsize(i) << std::endl;
        fprintv(i);
      }
      assert(fsize(i) < 4);
    }
  }

  void inflate_verts(size_t N) {
    __vert_begin.resize(__vert_begin.size() + N, -1);
  }

  void inflate_faces(size_t N) {
    __face_begin.resize(__face_begin.size() + N, -1);
  }

  void inflate_edge_pairs(size_t N) {
    size_t Ns = __corners_next.size() + 2 * N;
    __corners_next.resize(Ns, -1);
    __corners_prev.resize(Ns, -1);
    __corners_vert.resize(Ns, -1);
    __corners_face.resize(Ns, -1);
  }

  index_t insert_vertex() {
    // returns id of new vertex
    __vert_begin.push_back(-1);
    return __vert_begin.size() - 1;
  }

  index_t insert_face() {
    // returns id of new vertex
    __face_begin.push_back(-1);
    return __face_begin.size() - 1;
  }

  index_t insert_edge_pair() {
    // returns id of new vertex
    __corners_next.push_back(-1);
    __corners_next.push_back(-1);
    __corners_prev.push_back(-1);
    __corners_prev.push_back(-1);
    __corners_vert.push_back(-1);
    __corners_vert.push_back(-1);
    __corners_face.push_back(-1);
    __corners_face.push_back(-1);

    return __corners_next.size() - 2;
  }

  void swap_rows(index_t cA, index_t cB) {

    index_t cAn = __corners_next[cA];
    index_t cAp = __corners_prev[cA];
    index_t vA = __corners_vert[cA];
    index_t fA = __corners_face[cA];

    index_t cBn = __corners_next[cB];
    index_t cBp = __corners_prev[cB];
    index_t vB = __corners_vert[cB];
    index_t fB = __corners_face[cB];

    set_next(prev(cA), cB);
    set_prev(next(cA), cB);
    set_vbegin(vA, cB);
    set_fbegin(fA, cB);

    set_next(prev(cB), cA);
    set_prev(next(cB), cA);
    set_vbegin(vB, cA);
    set_fbegin(fB, cA);

    __corners_next[cA] = cBn;
    __corners_prev[cA] = cBp;
    __corners_vert[cA] = vB;
    __corners_face[cA] = fB;

    __corners_next[cB] = cAn;
    __corners_prev[cB] = cAp;
    __corners_vert[cB] = vA;
    __corners_face[cB] = fA;
  }

  void flip_edge(index_t c0) {
    index_t c1 = other(c0);
    swap_rows(c0, c1);
  }

  void remove_vertex(index_t i) {
    // std::cout << "    " << __PRETTY_FUNCTION__ << " " << i << std::endl;
    __vert_begin[i] = -1;
  }

  void remove_face(index_t i) {
    // std::cout << "    " << __PRETTY_FUNCTION__ << " " << i << std::endl;
    __face_begin[i] = -1;
  }

  void remove_edge_pair(index_t i0) {
    // returns id of new vertex
    index_t i1 = other(i0);
    __corners_next[i0] = -1;
    __corners_next[i1] = -1;

    __corners_prev[i0] = -1;
    __corners_prev[i1] = -1;

    __corners_vert[i0] = -1;
    __corners_vert[i1] = -1;

    __corners_face[i0] = -1;
    __corners_face[i1] = -1;
  }

  void for_each_face(index_t i,
                     std::function<void(index_t cid, shell &m)> func) {
    int j0 = this->fbegin(i);
    int j_end = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j0 != j_end;
      func(j0, *this);
      j0 = this->next(j0);
    }
  }

  void const_for_each_face(
      index_t i, std::function<void(index_t cid, const shell &m)> func) const {
    int j0 = this->fbegin(i);
    int j_end = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j0 != j_end;
      func(j0, *this);
      j0 = this->next(j0);
    }
  }

  void for_each_face_tri(
      index_t i,
      std::function<void(index_t c0, index_t c1, index_t c2, shell &m)> func) {
    int j1 = this->fbegin(i);
    int j2 = this->next(j1);
    int j0 = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j2 != j0;
      func(j0, j1, j2, *this);
      index_t jn = this->next(j2);
      j1 = j2;
      j2 = jn;
    }
  }

  void const_for_each_face_tri(
      index_t i,
      std::function<void(index_t c0, index_t c1, index_t c2, const shell &m)>
          func) const {
    int j1 = this->fbegin(i);
    int j2 = this->next(j1);
    int j0 = this->fend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j2 != j0;
      func(j0, j1, j2, *this);
      index_t jn = this->next(j2);
      j1 = j2;
      j2 = jn;
    }
  }

  void for_each_vertex(index_t i,
                       std::function<void(index_t cid, shell &m)> func) {

    int j0 = this->vbegin(i);
    int j_end = this->vend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j0 != j_end;
      func(j0, *this);
      j0 = this->vnext(j0);
    }
  }

  void const_for_each_vertex(
      index_t i, std::function<void(index_t cid, const shell &m)> func) const {

    int j0 = this->vbegin(i);
    int j_end = this->vend(i);
    bool it = true;
    int k = 0;
    while (it && k++ < 100) {
      it = j0 != j_end;
      func(j0, *this);
      j0 = this->vnext(j0);
    }
  }

  std::vector<index_t> get_edge_range() const {
    std::vector<index_t> range;
    range.reserve(corner_count() / 2);
    // replace this with some c++isms
    for (int i = 0; i < corner_count(); i += 2) {
      if (__corners_next[i] > -1)
        range.push_back(i);
    }
    return range;
  }

  std::vector<index_t> get_edge_range_2() const {
    std::vector<index_t> range;
    range.reserve(corner_count() / 2);
    // replace this with some c++isms
    index_t ii = 0;
    for (int i = 0; i < corner_count(); i += 2) {
      if (__corners_next[i] > -1)
        range.push_back(ii++);
    }
    return range;
  }

  std::vector<index_t> get_vert_range() const {
    std::vector<index_t> range;
    range.reserve(vert_count());
    // replace this with some c++isms
    for (int i = 0; i < vert_count(); i++) {
      if (__vert_begin[i] > -1)
        range.push_back(i);
    }
    return range;
  }

  std::vector<index_t> get_face_range(bool filter_tris = true) const {
    std::vector<index_t> range;
    range.reserve(face_count());
    // replace this with some c++isms
    for (int i = 0; i < face_count(); i++) {
      if (fbegin(i) < 0)
        continue;
      if (filter_tris && fsize(i) != 3)
        continue;
      range.push_back(i);
    }
    return range;
  }

  std::vector<index_t> get_edge_vert_ids() {
    std::vector<index_t> range;
    range.reserve(corner_count());
    // replace this with some c++isms
    for (int i = 0; i < corner_count(); i += 2) {
      if (__corners_next[i] < 0)
        continue;
      range.push_back(vert(i));
      range.push_back(vert(other(i)));
    }
    return range;
  }

  std::vector<index_t> get_face_vert_ids(bool filter_tris = true) {
    std::vector<int> faces;
    faces.reserve(3 * face_count());

    for (int i = 0; i < face_count(); i++) {
      if (fbegin(i) < 0)
        continue;

      if (filter_tris && fsize(i) != 3)
        continue;

      for_each_face(
          i, [&faces](int ci, shell &M) { faces.push_back(M.vert(ci)); });
    }
    return faces;
  }

  std::vector<index_t> get_one_ring(index_t iv) const {
    std::vector<index_t> range;
    this->const_for_each_vertex(iv, [&range](int ci, const shell &M) {
      range.push_back(M.vert(M.next(ci)));
    });

    return range;
  }

  std::array<index_t, 3> get_tri(index_t i) const {
    std::array<index_t, 3> tri = {-1, -1, -1};
    index_t j = 0;
    this->const_for_each_face(
        i, [&tri, &j](int ci, const shell &M) { tri[j++] = M.vert(ci); });
    return tri;
  }

  std::vector<index_t> get_range_map(const std::vector<index_t> indices,
                                     int stride, bool filter_tris = true) {
    std::vector<index_t> map;
    map.reserve(indices.size() / stride);
    // replace this with some c++isms
    for (int i = 0; i < indices.size(); i += stride) {
      if (indices[i] < 0)
        continue;

      if (filter_tris && fsize(i) != 3)
        continue;

      map.push_back(i);
    }
    return map;
  }

  std::vector<index_t> get_edge_map() {
    return get_range_map(__corners_next, 2, false);
  }
  std::vector<index_t> get_vert_map() {
    return get_range_map(__vert_begin, 1, false);
  }
  std::vector<index_t> get_face_map(bool filter_tris = true) {
    return get_range_map(__face_begin, 1, filter_tris);
  }

  void fupdate(index_t f) {
    for_each_face(f, [f](index_t cid, shell &m) { m.set_face(cid, f); });
  }

  void vupdate(index_t v) {
    for_each_vertex(v, [v](index_t cid, shell &m) { m.set_vert(cid, v); });
  }

  size_t fsize(index_t v) const {
    size_t s = 0;
    const_for_each_face(v, [&s](index_t cid, const shell &m) { s++; });
    return s;
  }

  size_t vsize(index_t v) const {
    size_t s = 0;
    const_for_each_vertex(v, [&s](index_t cid, const shell &m) { s++; });
    return s;
  }

  void cprint(index_t c0) {
    index_t c1 = other(c0);
    index_t c0p = prev(c0);
    index_t c1p = prev(c1);

    index_t v0 = vert(c0);
    index_t v1 = vert(c1);
    index_t v0p = vert(c0p);
    index_t v1p = vert(c1p);
    std::cout << "// =========== " << std::endl;
    std::cout << "//    " << c0p << std::endl;
    std::cout << "//" << c0 << "---" << c1 << std::endl;
    std::cout << "//    " << c1p << std::endl;
    std::cout << "//    " << v0p << std::endl;
    std::cout << "//" << v0 << "---" << v1 << std::endl;
    std::cout << "//    " << v1p << std::endl;
    std::cout << "// vs: " << vsize(v0) << " " << vsize(v1) << std::endl;
    std::cout << "// vps: " << vsize(v0) << " " << vsize(v1) << std::endl;
    vprint_graph_viz(v0);
    vprint_graph_viz(v1);
    vprint_graph_viz(v0p);
    vprint_graph_viz(v1p);
    vprintvs(v0);
    vprintvs(v1);
    vprintvs(v0p);
    vprintvs(v1p);

    std::cout << "// =========== " << std::endl;
  }

  void fprint(index_t f) {
    std::cout << "//f-" << f << ": ";
    for_each_face(f, [](index_t cid, shell &m) { std::cout << cid << " "; });
    std::cout << std::endl;
  }

  void fprintv(index_t f) {
    std::cout << "//fvs-" << f << ": " << std::flush;
    for_each_face(f, [this](index_t cid, shell &m) {
      std::cout << this->vert(cid) << " " << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprint(index_t v) {
    std::cout << "//v-" << v << ": " << std::flush;
    for_each_vertex(v, [](index_t cid, shell &m) {
      std::cout << cid << " " << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprintv(index_t v) {
    std::cout << "//vvs-" << v << ": " << std::flush;
    for_each_vertex(v, [this](index_t cid, shell &m) {
      std::cout << this->vert(this->next(cid)) << " " << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprintvs(index_t v) {
    std::cout << "//vvs-" << v << ": " << std::flush;
    for_each_vertex(v, [this](index_t cid, shell &m) {
      std::cout << this->vsize(this->vert(this->next(cid))) << " "
                << std::flush;
    });
    std::cout << std::endl << std::flush;
  }

  void vprintfs(index_t v) {
    std::cout << "//vfs-" << v << ": ";
    for_each_vertex(v, [this](index_t cid, shell &m) {
      std::cout << this->fsize(this->face(cid)) << " ";
    });
    std::cout << std::endl;
  }

  void vprint_graph_viz(index_t v) {
    std::cout << "//    ========" << std::endl;
    std::cout << "//    vg-" << v << std::endl;
    for_each_vertex(v, [this](index_t cid, shell &m) {
      // std::cout << "    " << cid << " -> " << this->next(cid) << std::endl;
      // std::cout << "    " << cid << " -> v" << this->vert(cid) <<
      // std::endl;
      std::cout << "    " << this->vert(cid) << " -> "
                << this->vert(this->next(cid)) << std::endl;
    });
    std::cout << "//    ========" << std::endl;
  }

  // accessors
  std::vector<index_t> &corners_next() { return __corners_next; }
  std::vector<index_t> &corners_prev() { return __corners_prev; }
  std::vector<index_t> &corners_vert() { return __corners_vert; }
  std::vector<index_t> &corners_face() { return __corners_face; }
  std::vector<index_t> &vert_begin() { return __vert_begin; }
  std::vector<index_t> &face_begin() { return __face_begin; }

  std::vector<index_t> __corners_next;
  std::vector<index_t> __corners_prev;
  std::vector<index_t> __corners_vert;
  std::vector<index_t> __corners_face;

  std::vector<index_t> __vert_begin;
  std::vector<index_t> __face_begin;
  std::vector<datum_ptr> __data;
};
} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif