

#include "geometry_types.hpp"
#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include <cstddef>
#include <memory.h>
#include <ostream>
#include <stdio.h>
#include <vector>
#include <zlib.h>

#ifndef __ASAWA_MANIFOLD_SURFACE__
#define __ASAWA_MANIFOLD_SURFACE__
namespace asawa {

typedef int index_t;

enum prim_type {
  VERTEX,
  FACE,
  EDGE,
  CORNER,
};

class datum;
typedef std::shared_ptr<datum> datum_ptr;

class manifold {
public:
  typedef std::shared_ptr<manifold> ptr;

  static ptr create(const std::vector<index_t> &corners_next,
                    const std::vector<index_t> &corners_vert,
                    const std::vector<index_t> &corners_face) {

    return std::make_shared<manifold>(corners_next, corners_vert, corners_face);
  }

  manifold(const std::vector<index_t> &corners_next,
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
  /*
  void pack() {
    // TODO: safer pack, iterating from mRecycle[i] to mRecycle[i+1]
    // if (mFaceRecycle.size() > 0) {
    face_array tFaces;

    int j = 0;
    for (int i = 0; i < mFaces.size(); i++) {
      if (mHasFace[i]) {
        tFaces.push_back(mFaces[i]);
        tFaces.back()->position_in_set() = j;
        j++;
      }
    }
    mHasFace = std::vector<bool>(tFaces.size(), true);
    swap(mFaces, tFaces);
  }
  */
  index_t insert_datum(datum_ptr datum) {
    __data.push_back(datum);
    return __data.size() - 1;
  }

  datum_ptr &get_datum(index_t i) { return __data[i]; }
  std::vector<datum_ptr> &get_data() { return __data; }

  size_t vert_count() const { return __vert_begin.size(); }
  size_t face_count() const { return __face_begin.size(); }
  size_t corner_count() const { return __corners_next.size(); }

  index_t other(int id) const { return 2 * (id / 2) + (id + 1) % 2; };
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

  void set_vert_pair(index_t c0, index_t v0, index_t v1) {
    index_t c1 = other(c0);
    set_vert(c0, v0);
    set_vert(c1, v1);
  }

  void set_face(index_t cid, index_t f) {
    __corners_face[cid] = f;
    set_fbegin(f, cid);
  }

  void link(index_t c0, index_t c1) {
    set_next(c0, c1);
    set_prev(c1, c0);
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

  void remove_vertex(index_t i) {
    std::cout << "    " << __FUNCTION__ << " " << i << std::endl;
    __vert_begin[i] = -1;
  }

  void remove_face(index_t i) {
    std::cout << "    " << __FUNCTION__ << " " << i << std::endl;
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
                     std::function<void(index_t cid, manifold &m)> func) {
    int j0 = this->fbegin(i);
    int j_end = this->fend(i);
    bool it = true;
    while (it) {
      it = j0 != j_end;
      func(j0, *this);
      j0 = this->next(j0);
    }
  }

  void for_each_vertex(index_t i,
                       std::function<void(index_t cid, manifold &m)> func) {
    int j0 = this->vbegin(i);
    int j_end = this->vend(i);
    bool it = true;
    while (it) {
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

  void fupdate(index_t f) {
    for_each_face(f, [f](index_t cid, manifold &m) { m.set_face(cid, f); });
  }

  void vupdate(index_t v) {
    for_each_vertex(v, [v](index_t cid, manifold &m) { m.set_vert(cid, v); });
  }

  size_t fsize(index_t v) {
    size_t s = 0;
    for_each_face(v, [&s](index_t cid, manifold &m) { s++; });
    return s;
  }

  size_t vsize(index_t v) {
    size_t s = 0;
    for_each_vertex(v, [&s](index_t cid, manifold &m) { s++; });
    return s;
  }

  void fprint(index_t f) {
    std::cout << "f" << f << ": ";
    for_each_face(f, [](index_t cid, manifold &m) { std::cout << cid << " "; });
    std::cout << std::endl;
  }

  void fprintv(index_t f) {
    std::cout << "f valence:" << f << ": ";
    for_each_face(f, [this](index_t cid, manifold &m) {
      std::cout << this->vert(cid) << " ";
    });
    std::cout << std::endl;
  }

  void vprint(index_t v) {
    std::cout << "v" << v << ": ";
    for_each_vertex(v,
                    [](index_t cid, manifold &m) { std::cout << cid << " "; });
    std::cout << std::endl;
  }

  std::vector<index_t> __corners_next;
  std::vector<index_t> __corners_prev;
  std::vector<index_t> __corners_vert;
  std::vector<index_t> __corners_face;

  std::vector<index_t> __vert_begin;
  std::vector<index_t> __face_begin;
  std::vector<datum_ptr> __data;
};

struct datum {
public:
  typedef std::shared_ptr<datum> ptr;

  // static ptr create() { return std::make_shared<datum>(); }

  datum(prim_type type) : __type(type){};
  virtual ~datum(){};

  template <class... Types>
  void calc(index_t i, const manifold &M, Types... args) {
    std::vector<index_t> vector;
    add_args(vector, args...);
    assert(!empty(vector));
    do_calc(i, vector, M);
  }

  template <typename LastType>
  static void add_args(std::vector<index_t> &vector, LastType arg) {
    vector.push_back(arg);
  };

  template <typename FirstType, typename... OtherTypes>
  static void add_args(std::vector<index_t> &vector, FirstType const &firstArg,
                       OtherTypes... otherArgs) {
    vector.push_back(firstArg);
    add_args(vector, otherArgs...);
  };

  void alloc(size_t sz) { do_alloc(sz); }
  void resize(size_t sz) { do_resize(sz); }
  void map(index_t i, index_t it) { do_map(i, it); }

  virtual void do_alloc(const size_t &sz) = 0;
  virtual void do_resize(const size_t &sz) = 0;
  virtual void do_calc(const index_t &i, const std::vector<index_t> &vals,
                       const manifold &M) = 0;
  virtual void do_map(const index_t i, const index_t it) = 0;

  virtual void write(FILE *file) const = 0;

  virtual void compute_vert_split(size_t v0, size_t v1) {}

  virtual void read(FILE *file) = 0;
  virtual void clear() = 0;
  const prim_type &type() { return __type; };

  prim_type __type;
};

template <typename TYPE> struct datum_t : public datum {
public:
  typedef std::shared_ptr<datum_t<TYPE>> ptr;

  static ptr create(prim_type type, const std::vector<TYPE> &data) {
    return std::make_shared<datum_t<TYPE>>(type, data);
  }

  datum_t(prim_type type, const std::vector<TYPE> &data)
      : datum(type), __data(data){};
  virtual ~datum_t(){};

  std::vector<TYPE> &data() { return __data; }
  const std::vector<TYPE> &data() const { return __data; }

  virtual void do_calc(const index_t &i, const std::vector<index_t> &vals,
                       const manifold &M) {
    index_t ic0 = vals[0];
    index_t ic1 = M.other(ic0);
    index_t iv0 = M.vert(ic0);
    index_t iv1 = M.vert(ic1);

    TYPE v0 = __data[iv0];
    TYPE v1 = __data[iv1];

    __tmp[i] = 0.5 * (v0 + v1);
  }
  virtual void do_alloc(const size_t &sz) {
    __tmp = std::vector<TYPE>(sz, z::zero<TYPE>());
  }
  virtual void do_resize(const size_t &sz) { __data.resize(sz); }
  virtual void do_map(const index_t i, const index_t it) {
    __data[i] = __tmp[it];
  }

  unsigned long operator[](int i) const { return __data[i]; }
  unsigned long &operator[](int i) { return __data[i]; }

  virtual void write(FILE *file) const {
    size_t nData = __data.size();
    size_t e;
    e = fwrite((void *)&nData, sizeof(size_t), 1, file);
    e = fwrite((void *)__data.data(), sizeof(TYPE), nData, file);
  }

  virtual void read(FILE *file) {
    size_t nData;
    size_t e;
    e = fread((void *)&nData, sizeof(size_t), 1, file);
    __data.resize(nData);
    e = fread((void *)__data.data(), sizeof(TYPE), nData, file);
  }

  virtual void clear() { __data.clear(); }

  std::vector<TYPE> __data;
  std::vector<TYPE> __tmp;
};

} // namespace asawa
#endif