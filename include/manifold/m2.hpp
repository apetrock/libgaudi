/*
 *  m2Common.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 1/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __TWOMANIFOLDCOMMON__
#define __TWOMANIFOLDCOMMON__

#include <algorithm>
#include <assert.h>
#include <bitset>
#include <cstddef>
#include <iostream>
#include <list>
#include <math.h>
#include <stack>
#include <vector>

#include <zlib.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
//#include <GLUT/glut.h>
#else
#ifdef _WIN32
#include <windows.h>
#endif
//#include <GL/gl.h>
//#include <GL/glu.h>
//#include <GL/glew.h>
//#include <GL/glut.h>
#endif

#include <any>
#include <memory>

#include "geometry_types.hpp"
#include "trace.hpp"
#include "vec_addendum.h"
// this typedef list is ugly but useful!!!
#define M2_TYPEDEFS                                                            \
  typedef typename SPACE::real real;                                           \
  typedef typename SPACE::complex complex;                                     \
  typedef typename SPACE::coordinate_type coordinate_type;                     \
  typedef typename SPACE::swept_point_type swept_point_type;                   \
  typedef typename SPACE::swept_triangle_type swept_triangle_type;             \
  typedef typename SPACE::box_type box_type;                                   \
  typedef typename SPACE::quat quat;                                           \
  typedef typename SPACE::mat2 mat2;                                           \
  typedef typename SPACE::mat3 mat3;                                           \
  typedef typename SPACE::mat4 mat4;                                           \
  typedef typename SPACE::mat43 mat43;                                         \
  typedef typename SPACE::double_type T;                                       \
  typedef m2::edge_line<SPACE> line_type;                                      \
  typedef m2::edge_line_pair<SPACE> line_pair;                                 \
  typedef m2::face_triangle<SPACE> triangle_type;                              \
  typedef m2::face_triangle_pair<SPACE> triangle_pair;                         \
                                                                               \
  typedef m2::surf<SPACE> surf_type;                                           \
  typedef typename m2::surf<SPACE>::face face_type;                            \
  typedef typename m2::surf<SPACE>::edge edge_type;                            \
  typedef typename m2::surf<SPACE>::vertex vertex_type;                        \
  typedef typename m2::surf<SPACE>::face_vertex face_vertex_type;              \
                                                                               \
  typedef surf_type *surf_ptr;                                                 \
  typedef face_type *face_ptr;                                                 \
  typedef edge_type *edge_ptr;                                                 \
  typedef vertex_type *vertex_ptr;                                             \
  typedef face_vertex_type *face_vertex_ptr;                                   \
                                                                               \
  typedef surf_type &surf_ref;                                                 \
  typedef face_type &face_ref;                                                 \
  typedef edge_type &edge_ref;                                                 \
  typedef vertex_type &vertex_ref;                                             \
  typedef face_vertex_type &face_vertex_ref;                                   \
                                                                               \
  typedef vector<coordinate_type> coordinate_array;                            \
  typedef vector<face_ptr> face_array;                                         \
  typedef vector<edge_ptr> edge_array;                                         \
  typedef vector<vertex_ptr> vertex_array;                                     \
  typedef vector<face_vertex_ptr> face_vertex_array;                           \
                                                                               \
  typedef list<coordinate_type> coordinate_list;                               \
  typedef list<face_ptr> face_list;                                            \
  typedef list<edge_ptr> edge_list;                                            \
  typedef list<vertex_ptr> vertex_list;                                        \
  typedef list<face_vertex_ptr> face_vertex_list;                              \
                                                                               \
  typedef typename coordinate_array::iterator ca_iterator;                     \
  typedef typename vector<face_ptr>::iterator fa_iterator;                     \
  typedef typename vector<edge_ptr>::iterator ea_iterator;                     \
  typedef typename vector<vertex_ptr>::iterator va_iterator;                   \
  typedef typename vector<face_vertex_ptr>::iterator fva_iterator;             \
                                                                               \
  typedef typename list<face_ptr>::iterator fl_iterator;                       \
  typedef typename list<edge_ptr>::iterator el_iterator;                       \
  typedef typename list<vertex_ptr>::iterator vl_iterator;                     \
  typedef typename list<face_vertex_ptr>::iterator fvl_iterator;

using namespace std;

namespace m2 {

template <typename SPACE> class surf;
template <typename SPACE> struct edge_line;
template <typename SPACE> struct edge_line_pair;
template <typename SPACE> struct face_triangle;
template <typename SPACE> struct face_triangle_pair;

template <typename SPACE> struct flattened_surf;

template <typename SPACE> struct edge_line : SPACE::line_type {

public:
  M2_TYPEDEFS;

  int indices[2];
  int edgeId;
  edge_line() {
    indices[0] = -1;
    indices[1] = -1;
    edgeId = -1;
  };

  edge_line(coordinate_type p0, coordinate_type p1, int i0, int i1, int eid) {
    this->setP(p0, p1);
    indices[0] = i0;
    indices[1] = i1;
    edgeId = eid;
  };

  T dist(const edge_line &that) const {
    T d = edgeId == that.edgeId ? std::numeric_limits<T>::infinity()
                                : sqrt(this->avgSqdDist(that));

    return d;
  };
};

template <typename SPACE>
inline bool operator==(const edge_line<SPACE> &lhs,
                       const edge_line<SPACE> &rhs) {
  return (lhs.indices[0] == rhs.indices[0] &&
          lhs.indices[1] == rhs.indices[1] && lhs.edgeId == rhs.edgeId);
}

template <typename SPACE> struct edge_line_pair {
public:
  M2_TYPEDEFS;

  line_type A;
  line_type B;
  mutable T _norm = -1;
  mutable T _angle = -3;

  edge_line_pair(){};

  edge_line_pair(const line_type &a, const line_type &b) {
    A = a;
    B = b;
  };

  T dist() const {
    if (A.edgeId == B.edgeId)
      return std::numeric_limits<T>::infinity();

    if (_norm < 0)
      _norm = A.dist(B);

    return _norm;
  };
};

template <typename SPACE>
inline bool operator==(const edge_line_pair<SPACE> &lhs,
                       const edge_line_pair<SPACE> &rhs) {
  return (lhs.A == rhs.A && lhs.B == rhs.B) ||
         (lhs.A == rhs.B && lhs.B == rhs.A);
}

template <typename SPACE>
inline bool operator<(const edge_line_pair<SPACE> &lhs,
                      const edge_line_pair<SPACE> &rhs) {
  return lhs.dist() < rhs.dist();
  // return lhs.angle() < rhs.angle();
}

template <typename SPACE> struct face_triangle : SPACE::triangle_type {

public:
  M2_TYPEDEFS;

  int indices[3];
  int faceId;
  face_triangle() {
    indices[0] = -1;
    indices[1] = -1;
    indices[2] = -1;
    faceId = -1;
  };

  face_triangle(coordinate_type p0, coordinate_type p1, coordinate_type p2,
                int i0, int i1, int i2, int fid) {
    this->p[0] = p0;
    this->p[1] = p1;
    this->p[2] = p2;
    indices[0] = i0;
    indices[1] = i1;
    indices[2] = i2;
    faceId = fid;
  };

  T dist(const face_triangle &that) const {
    if (faceId == that.faceId)
      return std::numeric_limits<T>::infinity();
    return sqrt(this->avgSqdDist(that));
  };
};

template <typename SPACE>
inline bool operator==(const face_triangle<SPACE> &lhs,
                       const face_triangle<SPACE> &rhs) {
  return (lhs.indices[0] == rhs.indices[0] &&
          lhs.indices[1] == rhs.indices[1] &&
          lhs.indices[2] == rhs.indices[2] && lhs.faceId == rhs.faceId);
}

template <typename SPACE> struct face_triangle_pair {
public:
  M2_TYPEDEFS;

  triangle_type A;
  triangle_type B;
  mutable T _norm = -1;
  mutable T _angle = -3;

  face_triangle_pair(const triangle_type &a, const triangle_type &b) {
    A = a;
    B = b;
  };

  T dist() const {
    if (A.faceId == B.faceId)
      return std::numeric_limits<T>::infinity();

    if (_norm < 0)
      _norm = A.dist(B);

    return _norm;
  };

  T angle() const {
    if (A.faceId == B.faceId)
      return std::numeric_limits<T>::infinity();

    if (_angle < -2) {
      _angle = A.angle(B);
    }

    return _angle;
  };
};

template <typename SPACE>
inline bool operator==(const face_triangle_pair<SPACE> &lhs,
                       const face_triangle_pair<SPACE> &rhs) {
  return (lhs.A == rhs.A && lhs.B == rhs.B) ||
         (lhs.A == rhs.B && lhs.B == rhs.A);
}

template <typename SPACE>
inline bool operator<(const face_triangle_pair<SPACE> &lhs,
                      const face_triangle_pair<SPACE> &rhs) {
  return lhs.dist() < rhs.dist();
  // return lhs.angle() < rhs.angle();
}

struct colorRGB {
  double r, g, b, a;
  colorRGB() {
    r = 0.5;
    g = 0.5;
    b = 0.5;
    a = 1.0;
  }
  colorRGB(double rr, double gg, double bb, double aa) {
    r = rr;
    g = gg;
    b = bb;
    a = aa;
  }
};

/*
inline void gl_set_color(colorRGB in){
  glColor4f(in.r, in.g, in.b, in.a);
}
inline void gl_set_color_opaque(colorRGB in){
  glColor4f(in.r, in.g, in.b, 1.0);
}
*/
/*
template <typename T>
ostream &operator<<( ostream &out, const al::Vec<3,T>& in ) {
  out << " " << in[0] << " " <<  in[1] << " " <<  in[2] << " ";
  return out;
}
*/

template <typename T> inline T randd(T range_) {
  T output = (T)(rand() / (T(RAND_MAX) + 1.0)) * range_;
  return output;
}
/*
template <typename T>
inline Eigen::Matrix<T,4,1> cross(Eigen::Matrix<T,4,1> a,Eigen::Matrix<T,4,1>
b){ return Eigen::Matrix<T,4,1>( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y
- a.y*b.x, 1.0);
}
*/

template <typename SPACE> class surf {
public:
  template <class ITYPE> class data_node;
  class face_vertex;
  class face;
  class edge;
  class vertex;

  M2_TYPEDEFS

  template <class ITYPE> class data_node {
  public:
    data_node() { flags.reset(); }

    data_node(const data_node &other)
        : _ddata(other._ddata), _idx_size(other._idx_size), flags(other.flags) {
      this->set_dirty();
    }

    std::any &operator[](ITYPE i) { return _ddata[static_cast<int>(i)]; }

    const std::any &operator[](ITYPE i) const {
      return _ddata[static_cast<int>(i)];
    }

    int data_size() const { return _ddata_size; }
    int index_size() const { return _idx_size; }

    template <typename TYPE> TYPE get(const ITYPE &i) {
      int ii = static_cast<int>(i);

      if (!_ddata[ii].has_value()) {
        _ddata[ii] = z::zero<TYPE>();
      }
      return std::any_cast<TYPE>(_ddata[ii]);
    }

    template <typename TYPE> void set(const ITYPE &i, const TYPE &d) {
      _ddata[static_cast<int>(i)] = d;
      _ddata_size = sizeof(TYPE);
      //_ddata[static_cast<int>(i)] = datum_t<TYPE>(d);
    }

    void set_dirty(const ITYPE &i, bool val) {
      dirty[static_cast<int>(i)] = val;
    }

    void get_dirty(const ITYPE &i) { return dirty[static_cast<int>(i)]; }

    void set_dirty() {
      for (int i = 0; i < index_size(); i++) {
        dirty[i] = true;
      }
    }

    virtual void set_changed() {
      this->set_dirty();
      this->_topology_changed = true;
    }

    virtual void _update() = 0;

    virtual void update() {
      if (_topology_changed)
        this->_update();
      _topology_changed = false;
    }

    std::bitset<8> flags;
    int topologyChangeId = -1; // needs to be changed to stored value

  private:
    bool _topology_changed = false;
    std::any _ddata[static_cast<int>(ITYPE::MAXINDEX)];
    size_t _ddata_size = 0;
    bool dirty[static_cast<int>(ITYPE::MAXINDEX)];
    int _idx_size = static_cast<int>(ITYPE::MAXINDEX);
  };

  class edge : public data_node<typename SPACE::edge_index> {

  public:
    edge() {
      fv1 = NULL;
      fv2 = NULL;
      mSetPosition = -1;
      flag = 0;
      delete_flag = 0;
    };

    edge(face_vertex_ptr in1, face_vertex_ptr in2) {
      mSetPosition = -1;
      fv1 = in1;
      fv2 = in2;
      in1->edge() = this;
      in2->edge() = this;
      flag = 0;
      delete_flag = 0;
    }

    edge_ref operator=(edge_ref rhs) {
      edge_ptr out = new edge_type(*rhs.fv1, *rhs.fv2);
      return *out;
    }

    ~edge(){};

    int &position_in_set() { return mSetPosition; }
    int position_in_set() const { return mSetPosition; }
    int fv_position_in_set(face_vertex_ptr cv) const {
      return (cv == fv1) ? 2 * mSetPosition : 2 * mSetPosition + 1;
    }

    int group() const { return mGroup; }
    int &group() { return mGroup; }

    face_vertex_ptr &v1() { return fv1; }
    face_vertex_ptr v1() const { return fv1; }
    face_vertex_ptr &v2() { return fv2; }
    face_vertex_ptr v2() const { return fv2; }
    virtual void _update() {}

    face_vertex_ptr &other(const face_vertex_ptr &cv) {
      // std::cout << " cv " << cv << std::endl;
      // std::cout << " fv1 " << fv1 << std::endl;
      // std::cout << " fv2 " << fv2 << std::endl;

      if (cv == fv1) {
        return fv2;
      } else
        return fv1;
    }

    face_vertex_ptr other(const face_vertex_ptr &cv) const {
      if (cv == fv1) {
        return fv2;
      } else
        return fv1;
    }

    face_vertex_ptr this_fv(face_vertex_ptr cv) const {
      if (cv != fv1) {
        return fv2;
      } else
        return fv1;
    }

    int side(face_vertex_ptr cv) const {
      if (cv == fv1) {
        return 0;
      } else if (cv == fv2)
        return 1;
      else
        return -1;
    }

    virtual void set_changed() {
      data_node<typename SPACE::edge_index>::set_changed();
      if (this->v1())
        this->v1()->set_changed();
      if (this->v2())
        this->v2()->set_changed();
    }

    void set(face_vertex_ptr nfv1, face_vertex_ptr nfv2) {
      // lets be a little careful and make sure we maintain order.
      // its probably not a problem, but...

      if (nfv1 == fv2 || nfv2 == fv1) {
        fv1 = nfv2;
        fv2 = nfv1;
      } else {
        fv1 = nfv1;
        fv2 = nfv2;
      }

      fv1 = nfv1;
      fv2 = nfv2;

      this->set_changed();
      nfv1->edge() = this;
      nfv2->edge() = this;
    }

    void set_this(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex) {
      if (this_vertex == fv1)
        fv1 = new_vertex;
      else
        fv2 = new_vertex;
      new_vertex->edge() = this;

      this->set_changed();
    }

    void set_other(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex) {
      if (this_vertex == fv1)
        fv2 = new_vertex;
      else
        fv1 = new_vertex;
      new_vertex->edge() = this;

      this->set_changed();
    }

    void set_other(face_ptr cv, face_vertex_ptr ov) {
      if (cv->face() == fv1->face()) {
        return fv2 = ov;
      } else
        fv1 = ov;

      this->set_changed();
    }

    void set_empty(face_vertex_ptr c) {
      if (fv1) {
        fv2 = c;
      } else
        fv1 = c;
      this->set_changed();
    }

    void swap_corners() {
      face_vertex_ptr fv1t = this->fv1;
      face_vertex_ptr fv2t = this->fv2;
      this->fv1 = fv2t;
      this->fv2 = fv1t;
    }

    void update_vertex(face_vertex_ptr old_, face_vertex_ptr new_) {
      if (fv1 == old_) {
        fv1 = new_;
      } else {
        fv2 = new_;
      }
    }

    face_ptr coface(face_vertex_ptr this_vert) const {
      if (this_vert == fv1) {
        return fv2->face();
      } else {
        return fv1->face();
      }
    }

    void print() {
      std::cout << "e:" << this->position_in_set() << std::endl;
      std::cout << " -fvs:" << this->v1() << " " << this->v2() << std::endl;
      std::cout << " -fs:" << this->v1()->face()->position_in_set() << " "
                << this->v2()->face()->position_in_set() << std::endl;
      std::cout << " -vs:" << this->v1()->vertex()->position_in_set() << " "
                << this->v2()->vertex()->position_in_set() << std::endl;
    };

    void verify() {
      this->update();
      assert(fv1->edge() == this);
      assert(fv2->edge() == this);
    };

    void update_face_vertices() {
      fv1->edge() = this;
      fv2->edge() = this;
    }

    int vnum(face_vertex_ptr v) {
      if (v == fv1)
        return 1;
      else if (v == fv2)
        return 2;
      else
        return 0;
    }

    bool is_degenerate() {
      face_vertex_ptr fv0 = this->v1();
      face_vertex_ptr fv1 = this->v2();
      vertex_ptr v00 = fv0->vertex();
      vertex_ptr v10 = fv1->vertex();
      // test these suckers with face vertices, more robusta?
      if (v00 == v10) {
        return true;
      }
      /*
      if (fv0->prev()->vnext()->edge() == fv1->prev()->vnext()->edge())
        // means they are tetrahedra
        return true;
      */
      if (fv0->face()->is_degenerate() || fv1->face()->is_degenerate()) {
        return true;
      }

      return false;
    }

    face_vertex_ptr fv1;
    face_vertex_ptr fv2;

  protected:
    int mSetPosition;
    int mGroup;

  public:
    unsigned int idata;
    unsigned int flag;
    unsigned int delete_flag;
  };

  class face : public data_node<typename SPACE::face_index> {

  public:
    face() {
      flag = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 0.75;
      fHead = NULL;
      mArea = 0.0;
      mSetPosition = -1;
      mSize = 0;
      flag = 0;

      data = coordinate_type(0, 0, 0);
      data2 = coordinate_type(0, 0, 0);
      data3 = 1.0;
    }

    face(vertex_ref pnt) {
      // this makes a the face in a hypersphere from one point
      face_vertex_ptr nfv = new face_vertex_type();
      mSize = 1;
      mArea = 0.0;

      nfv->vertex() = &pnt;
      nfv->face() = this;

      nfv->set_next(nfv);
      fHead = nfv;
      this->update_all();

      flag = 0;

      color.r = 0.8 + randd(0.2);
      color.g = 0.0 + randd(0.1);
      color.b = 0.0 + randd(0.1);
      color.a = 1.0;
      mSetPosition = -1;

      data = coordinate_type(0, 0, 0);
      data2 = coordinate_type(0, 0, 0);
      data3 = 1.0;
    }

    face(face_vertex_ptr &fv0, face_vertex_ptr &fv1, face_vertex_ptr &fv2) {
      // this makes a the face in a hypersphere from one point
      face_vertex_ptr nfv = new face_vertex_type();
      mSize = 1;
      mArea = 0.0;
      fv0->set_next(fv1);
      fv1->set_next(fv2);
      fv2->set_next(fv0);

      this->set_front(fv0);

      data = coordinate_type(0, 0, 0);
      data2 = coordinate_type(0, 0, 0);
      data3 = 1.0;
    }

    ~face() {
      // std::cout << "erasing face: " << ID() << std::endl;
    }

    face_ref operator=(face_ref rhs) {
      face_ptr out = new face_type();
      out->mNormal = rhs.mNormal;
      out->mCenter = rhs.mCenter;
      out->mVertices = rhs.mVertices;
      return *out;
    }

    bool is_degenerate() { return this->size() < 3; }

    bool is_null() { return this->size() == 1 && !fHead->edge(); }

    int &position_in_set() { return mSetPosition; }
    int position_in_set() const { return mSetPosition; }

    // TODO: delete setHead... actually make all fronts() set_fronts()
    void set_front(face_vertex_ptr head) {
      this->set_changed();
      head->face() = this;
      fHead = head;
    }

    face_vertex_ptr get_front() { return fHead; }

    face_vertex_ptr fbegin() { return get_front(); }
    face_vertex_ptr fend() { return get_front()->prev(); }
    face_vertex_ptr fendr() { return get_front()->next(); }

    void for_each(
        std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {

      if (!this->fbegin())
        return;
      if (!this->fend())
        return;

      face_vertex_ptr fvb = this->fbegin();
      face_vertex_ptr fve = this->fend();
      bool iterating = true;
      while (iterating) {
        iterating = fvb != fve;
        func(fvb);
        fvb = fvb->next();
      }
    }

    void for_each_reverse(
        std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {

      if (!this->fbegin())
        return;
      if (!this->fendr())
        return;

      face_vertex_ptr fvb = this->fbegin();
      face_vertex_ptr fve = this->fendr();
      bool iterating = true;
      while (iterating) {
        iterating = fvb != fve;
        func(fvb);
        fvb = fvb->prev();
      }
    }

    face_vertex_ptr get_corner_on_vertex(vertex_ptr v) {
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();

      face_vertex_ptr ito = fHead;
      bool at_head = false;
      while (!at_head) {
        if (itb == ite) {
          at_head = true;
        }
        if (itb->vertex() == v)
          ito = itb;
        itb = itb->next();
      }
      return ito;
    }

    face_vertex_array face_vertex_trace() {
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      face_vertex_array array_out;

      bool at_head = false;
      while (!at_head) {
        if (itb == ite) {
          at_head = true;
        }
        array_out.push_back(itb);
        itb = itb->next();
      }
      return array_out;
    }

    face_vertex_list vertex_trace() {
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      face_vertex_list array_out;

      bool at_head = false;
      while (!at_head) {
        if (itb == ite) {
          at_head = true;
        }
        array_out.push_back(itb);
        itb = itb->next();
      }
      return array_out;
    }

    void update_all() {
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();

      int i = 0;
      bool at_head = false;

      while (!at_head) {
        at_head = itb == ite;
        itb->face() = this;
        itb->position_in_face() = i;
        itb = itb->next();
        i++;
      }
      this->mSize = i;
    }

    int size() {
      this->update();
      return mSize;
    }

    virtual void _update() {
      // std::cout << " updating f: " << this->position_in_set() << std::endl;
      this->update_all();
    }

    bool has_vertex(vertex_ptr v) {
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();

      bool iterating = true;
      while (iterating) {
        iterating = itb != ite;

        if (itb->vertex() == v)
          return true;
        itb = itb->next();
      }
      return false;
    }

    void verify() {

      this->update();
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      bool iterating = true;
      // std::cout << " v: " << this->position_in_set() << " - ";
      while (iterating) {
        iterating = itb != ite;
        /*
        if (itb->face() != this) {
          trace::print_trace();
        }
        */
        // std::cout << itb->vID << " ";
        // itb->vertex()->verify();
        if (this->mSize > 1)
          assert(itb != itb->next());
        itb = itb->next();
      }
      // std::cout << endl;
    }

    void print() {

      this->update();
      cout << "face " << mSetPosition
           << ": number of vertices: " << this->size() << endl;
      cout << "vertices: " << endl;

      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      size_t new_id = 0;
      bool iterating = true;
      while (iterating) {
        iterating = itb != ite;
        cout << itb->vertex()->position_in_set() << ", ";
        itb = itb->next();
      }
      cout << endl;
    }

    void flag_vertices() {

      this->update();
      face_vertex_ptr fvA = this->fbegin();
      face_vertex_ptr fvAe = this->fend();
      bool has = false;
      bool itf = true;
      while (itf) {
        itf = fvA != fvAe;
        fvA->vertex()->flag = 1;
        fvA = fvA->next();
      }
    }

    void flag_edges(int flag) {

      this->update();
      face_vertex_ptr fvA = this->fbegin();
      face_vertex_ptr fvAe = this->fend();
      bool has = false;
      bool itf = true;
      while (itf) {
        itf = fvA != fvAe;
        fvA->edge()->flags[flag] = 1;
        fvA = fvA->next();
      }
    }

    void print_vert_ids(bool reverse = false) {
      face_vertex_ptr fvA = this->fbegin();
      face_vertex_ptr fvAe = reverse ? this->fendr() : this->fend();
      bool itf = true;
      while (itf) {
        itf = fvA != fvAe;
        std::cout << "   " << fvA->vertex()->position_in_set() << " ";
        fvA = reverse ? fvA->prev() : fvA->next();
      }
      std::cout << std::endl;
    }

    void print_shared_edge() {
      face_vertex_ptr fvA = this->fbegin();
      face_vertex_ptr fvAe = this->fend();
      bool itf = true;
      while (itf) {
        itf = fvA != fvAe;
        bool eA = fvA->edge()->v1()->vertex()->position_in_set() ==
                  fvA->edge()->v2()->vertex()->position_in_set();
        std::cout << "   " << eA;
        fvA = fvA->next();
      }
      std::cout << std::endl;
    }

    int group() const { return mGroup; }
    int &group() { return mGroup; }

  protected:
    // face_vertex_array	mVertices;
    int mGroup;

    face_vertex_ptr fHead;
    coordinate_type mCenter;
    coordinate_type mNormal;
    T mArea;

    int mSize;

    int mSetPosition;

    bool calc_normal;

  public:
    colorRGB color;
    colorRGB ncolor;
    unsigned int flag;
    coordinate_type data;
    coordinate_type data2;
    T data3;
  };

  class face_vertex : public data_node<typename SPACE::face_vertex_index> {

  public:
    face_vertex() : mVertex(NULL), mEdge(NULL), mFace(NULL) {
      vID = 0;
      fID = 0;
      flag = 0;
      data = 0;
    }

    face_vertex(face_vertex_ref rhs) {
      this->mEdge = rhs.mEdge;
      //			mEdge->set_this(&rhs,this);
      this->mFace = rhs.mFace;
      this->mVertex = rhs.mVertex;
      this->nxt_face = rhs.nxt_face;
      this->prv_face = rhs.prv_face;
      this->data = rhs.data;
      fID = 0;

      flag = 0;
    }

    ~face_vertex(){
        // mVertex->remove_face_vertex(mVertexPosition);
    };

    int position_in_set() { return mEdge->fv_position_in_set(this); }

    bool operator==(face_vertex_ref rhs) {
      if (mEdge == rhs.mEdge && mFace == rhs.mFace && mVertex == rhs.mVertex) {
        return true;
      } else
        return false;
    }

    face_vertex_ref operator=(face_vertex_ref rhs) {
      face_vertex_ptr out = new face_vertex_type(rhs);
      out->mEdge->update_vertex(&rhs, this);
      return *out;
    }

    virtual void _update() {}

    virtual void set_changed() {
      data_node<typename SPACE::face_vertex_index>::set_changed();
      if (this->face())
        this->face()->set_changed();
      if (this->vertex())
        this->vertex()->set_changed();
    }

    face_vertex_ptr coedge() { return mEdge->other(this); }

    size_t position_in_flat_set() {
      return 2 * mEdge->position_in_set() + mEdge->side(this);
    }

    face_ref get_face() { return *mFace; }

    face_vertex_ptr add_next() {
      face_vertex_ptr out = new face_vertex(*this);
      face_vertex_ptr nxt = this->nxt_face;
      out->set_next(nxt);
      this->set_next(out);

      out->face() = this->face();
      out->face()->set_front(out);
      this->vertex()->add_face_vertex(out);
      this->vertex()->set_front(out);

      this->set_changed();
      return out;
    }

    face_vertex_ptr add_prev() {
      face_vertex_ptr out = new face_vertex(*this);
      face_vertex_ptr prv = this->prv_face;
      prv->set_next(out);
      out->set_next(this);

      out->face() = this->face();
      out->face()->set_front(out);
      this->vertex()->add_face_vertex(out);
      this->vertex()->set_front(out);

      this->set_changed();
      return out;
    }

    void set_next(face_vertex_ptr fv) {
      fv->prev() = this;
      this->next() = fv;
      this->set_changed();
    }

    void set_prev(face_vertex_ptr fv) {
      fv->next() = this;
      this->prev() = fv;
      this->set_changed();
    }

    face_vertex_ptr &next() { return nxt_face; }
    face_vertex_ptr next() const { return nxt_face; }
    face_vertex_ptr &prev() { return prv_face; }
    face_vertex_ptr prev() const { return prv_face; }

    face_vertex_ptr vnext() {
      if (mEdge == NULL) {
        return this;
      } else {
        face_vertex_ptr out = mEdge->other(this);
        return out->next();
      }
    }

    face_vertex_ptr vprev() {
      if (this->prev()->edge() == NULL) {
        return this;
      } else {
        face_vertex_ptr out = this->prev();
        if (out->mEdge)
          return out->mEdge->other(out);
        else
          return this;
      }
    }

    edge_ptr &edge() { return mEdge; }
    edge_ptr edge() const { return mEdge; }
    face_ptr &face() { return mFace; }
    face_ptr face() const { return mFace; }
    face_ptr coface() const { return mEdge->other(this)->face(); }
    vertex_ptr &vertex() { return mVertex; }
    vertex_ptr vertex() const { return mVertex; }

    int face_ID() const { return this->mFacePosition; }
    int &face_ID() { return this->mFacePosition; }

    int &position_in_face() { return mFacePosition; }
    int position_in_face() const { return mFacePosition; }
    void set_edge(edge_ref input) { mEdge = &input; };
    void set_face(face_ref input) { mFace = &input; };
    void set_vertex(vertex_ref input) { mVertex = &input; };

    int &group() { return mGroup; }
    int group() const { return mGroup; }

    int vID = -1;

  protected:
    int mGroup;
    int fID = -1;
    face_vertex_ptr nxt_face;
    face_vertex_ptr prv_face;

    edge_ptr mEdge;
    face_ptr mFace;
    vertex_ptr mVertex;

    int mFacePosition;

  public:
    unsigned int flag;
    T data;
  };

  class vertex : public data_node<typename SPACE::vertex_index> {

  public:
    int graphColor;
    bool smooth;
    int isDegenerate;
    colorRGB color;

    vertex() {
      int graphColor = -1;

      flag = 0;
      mSize = 0;
      color.r = 0.4 + randd(0.1);
      color.g = 0.4 + randd(0.1);
      color.b = 0.4 + randd(0.1);
      color.a = 1.0;
      pinned = false;
      mFront = NULL;
      mSetPosition = -1;
      winding = 0;
      data2 = 0;
      isDegenerate = -1;
    }

    void init() {
      int graphColor = -1;
      face_ptr new_face = new face_type(*this);
      face_vertex_ptr new_fv = new_face->fbegin();
      new_fv->face() = new_face;
      this->add_face_vertex(new_fv);
      pinned = false;
      data2 = 0;

      this->set_changed();
    }

    int &position_in_set() { return mSetPosition; }

    int &group() { return mGroup; }
    int group() const { return mGroup; }

    int &ID() { return mSetPosition; }
    int ID() const { return mSetPosition; }

    int size() {
      this->update();
      return mSize;
    }

    face_vertex_ptr make_face_vertex() {
      face_vertex_ptr fv = new face_vertex_type();
      mFront = fv;
      fv->vertex() = this;
      this->set_changed();
      return fv;
    }

    void add_face_vertex(face_vertex_ptr new_fv) {
      new_fv->vertex() = this;
      this->set_changed();
    }

    void remove_face_vertex(face_vertex_ptr fv) {
      fv->vertex() = NULL;
      this->set_changed();
    }

    void update_all() {

      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      mSize = 0;
      while (iterating) {
        // std::cout << mSize << std::endl;
        mSize++;
        // std::cout << itb << std::endl;
        iterating = itb != ite;
        // std::cout << itb->vID << " " << itb->next()->vID << std::endl;
        itb->vertex() = this;

        itb->vID = itb->vertex()->position_in_set();

        itb = itb->vnext();
      }
    }

    virtual void _update() {
      // std::cout << " updating v: " << this->position_in_set() << std::endl;
      this->update_all();
    }

    void set_front(face_vertex_ptr fv) {
      this->set_changed();
      fv->vertex() = this;
      mFront = fv;
    }

    face_vertex_ptr get_front() { return mFront; }

    face_vertex_ptr fbegin() { return this->get_front(); }
    face_vertex_ptr fend() { return this->get_front()->vprev(); }

    face_vertex_ptr fendr() { return this->get_front()->vnext(); }

    bool is_degenerate() {
      bool degenerate = this->size() < 4;
      // if (degenerate)
      //   std::cout << "degenerate! " << this->size() << std::endl;
      return degenerate;
    }

    void for_each(
        std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
      using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

      if (!this->fbegin())
        return;
      if (!this->fend())
        return;

      face_vertex_ptr fvb = this->fbegin();
      face_vertex_ptr fve = this->fend();
      bool iterating = true;
      while (iterating) {
        iterating = fvb != fve;
        func(fvb);
        fvb = fvb->vnext();
      }
    }

    void for_each_reverse(
        std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
      using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

      if (!this->fbegin())
        return;

      if (!this->fendr())
        return;

      face_vertex_ptr fvb = this->fbegin();
      face_vertex_ptr fve = this->fendr();
      bool iterating = true;
      while (iterating) {
        iterating = fvb != fve;
        func(fvb);
        fvb = fvb->vprev();
      }
    }

    bool has_vertex(vertex_ptr v) {
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      bool iterating = true;
      while (iterating) {
        iterating = itb != ite;
        if (itb->next()->vertex() == v)
          return true;
        itb = itb->vnext();
      }
      return false;
    }

    bool shares_edge_with(vertex_ptr vi) {
      if (mSize == 0)
        return false;
      vertex_ptr v = this;
      face_vertex_ptr fvb = v->fbegin();
      face_vertex_ptr fve = v->fend();
      bool iterating = true;
      bool share = false;
      while (iterating) {
        iterating = fvb != fve;
        if (fvb->next()->vertex() == vi)
          share = true;
        fvb = fvb->vnext();
      }
      return share;
    }

    std::vector<edge_ptr> get_shared_edges(vertex_ptr vB) {
      int shared_vertices = 0;
      std::vector<edge_ptr> edges;
      face_vertex_ptr fvb = this->fbegin();
      face_vertex_ptr fve = this->fend();
      bool itA = true;
      int count = 0;
      // std::cout << " edges: " << std::endl;
      while (itA) {
        itA = fvb != fve;
        vertex_ptr vA = fvb->coedge()->vertex();
        /*
        std::cout << " va: " << fvb << " " <<fvb->coedge() << std::endl;
        std::cout << " va: " << vA->position_in_set() << std::endl;
        std::cout << " vb: " << vB->position_in_set() << std::endl;
        */
        if (vA == vB)
          edges.push_back(fvb->edge());
        fvb = fvb->vnext();
      }
      // std::cout << " done: " << edges.size() <<  std::endl;
      return edges;
    }

    edge_ptr get_shared_edge(vertex_ptr vi) {
      if (mSize == 0)
        return false;
      vertex_ptr v = this;
      face_vertex_ptr fvb = v->fbegin();
      face_vertex_ptr fve = v->fend();
      bool iterating = true;
      bool share = false;
      while (iterating) {
        iterating = fvb != fve;
        if (fvb->coedge()->vertex() == vi)
          return fvb->edge();
        fvb = fvb->vnext();
      }
      return NULL;
    }

    void verify() {
      this->update();
      if (mSize == 0)
        return;
      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      int i = 0;

      while (iterating) {
        iterating = itb != ite;
        if (itb->vertex() != this) {
          // trace::print_trace();

          std::cout << i++ << " asserted: " << position_in_set() << " v- "
                    << itb->vertex()->position_in_set() << " "
                    << " e- " << itb->edge()->position_in_set() << std::endl;
          std::cout << "          e eq: "
                    << bool(itb->edge()->v1() == itb->edge()->v2())
                    << std::endl;
          assert(itb != itb->vnext());
        }
        if (i++ > 20)
          break;
        itb = itb->vnext();
      }
    }

    void print() {

      this->update();
      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      std::cout << " - vertex: " << mSetPosition << ", size: " << this->mSize
                << std::endl;
      int i = 0;
      while (iterating) {
        if (i > 20)
          break;
        int ei = -1;
        if (itb->edge())
          ei = itb->edge()->position_in_set();

        iterating = itb != ite;
        // std::cout << itb->next()->vertex()->position_in_set() << ", ";
        // std::cout << "   " << i << ": v- " << itb->vID << ": vn- "
        //          << itb->next()->vID << std::endl;
        std::cout << "   " << i << ": v- " << itb->vertex()->position_in_set()
                  << " "
                  << " - e- " << ei << " - f- "
                  << itb->face()->position_in_set() << ": vn- "
                  << itb->next()->vertex()->position_in_set() << ": vp- "
                  << itb->prev()->vertex()->position_in_set() << std::endl;
        itb = itb->vnext();
        i++;
      }
      std::cout << std::endl;
    }

    void print_adj_sz() {
      if (mSize == 0)
        return;
      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      std::cout << " - vertex: " << mSetPosition << ", size: " << this->mSize
                << std::endl;
      int i = 0;
      while (iterating) {
        if (i > 20)
          break;
        iterating = itb != ite;
        std::cout << "    szn: " << itb->next()->vertex()->position_in_set()
                  << " " << itb->next()->vertex()->size() << std::endl;
        itb = itb->vnext();
        i++;
      }
      std::cout << std::endl;
    }

  protected:
    face_vertex_ptr mFront = NULL;

    int mSetPosition;
    int mSize;
    int mGroup;

  public:
    int pinned;
    unsigned int flag;
    coordinate_type data;
    T data2;
    T winding;
  };

  surf() {
    maxGraphColor = 0;
    manual_clean_up = false;
  }

  surf(const surf_ref rhs) {
    maxGraphColor = 0;
    manual_clean_up = false;
    // mFaces.clear();
    // mVertices.clear();
    // mEdges.clear();cl
    mFaces.resize(rhs.mFaces.size());
    mVertices.resize(rhs.mVertices.size());
    mEdges.resize(rhs.mEdges.size());

    for (int i = 0; i < rhs.mVertices.size(); i++) {
      coordinate_type nc(rhs.mVertices[i]->coordinate());
      mVertices[i] = new vertex_type(nc);
      mVertices[i]->position_in_set() = i;
    }

    for (int i = 0; i < rhs.mEdges.size(); i++) {
      mEdges[i] = new edge_type();
      mEdges[i]->position_in_set() = i;
    }

    for (int i = 0; i < rhs.mFaces.size(); i++) {
      if (rhs.mFaces[i]) {
        face_ptr nf = new face_type();
        nf->position_in_set() = i;

        face_ptr of = rhs.mFaces[i];
        nf->size() = of->size();
        face_vertex_ptr itb = of->fbegin();
        face_vertex_ptr ite = of->fend();
        face_vertex_ptr fv0, fv1;

        fv0 = new face_vertex_type();
        vertex_ptr nv0 = mVertices[ite->vertex()->position_in_set()];
        fv0->vertex() = nv0;
        nv0->add_face_vertex(fv0);

        edge_ptr ne0 = mEdges[ite->edge()->position_in_set()];
        fv0->edge() = ne0;
        fv0->face() = nf;
        nf->set_front(fv0);
        if (!ne0->v1())
          ne0->v1() = fv0;
        else
          ne0->v2() = fv0;

        bool iterating = true;
        while (iterating) {
          iterating = itb != ite->prev();
          fv1 = new face_vertex_type();
          vertex_ptr nv1 = mVertices[itb->vertex()->position_in_set()];
          fv1->vertex() = nv1;
          nv1->add_face_vertex(fv1);
          edge_ptr ne1 = mEdges[itb->edge()->position_in_set()];
          fv1->edge() = ne1;
          fv1->face() = nf;

          if (!ne1->v1())
            ne1->v1() = fv1;
          else
            ne1->v2() = fv1;

          fv0->next() = fv1;
          fv1->prev() = fv0;
          fv0 = fv1;
          itb = itb->next();
        }
        fv1 = nf->fbegin();
        fv0->next() = fv1;
        fv1->prev() = fv0;
        mFaces[i] = nf;
      }
    }
    mFaceRecycle = rhs.mFaceRecycle;
    mVertexRecycle = rhs.mVertexRecycle;
    mEdgeRecycle = rhs.mEdgeRecycle;
    bool here = true;
  }

  surf_ref operator=(const surf_ref rhs) {
    std::cout << "deep copy equals" << std::endl;
    // deep copy
    maxGraphColor = 0;
    manual_clean_up = false;
    mFaces.clear();
    mVertices.clear();
    mEdges.clear();
    mFaces.resize(rhs.mFaces.size());
    mVertices.resize(rhs.mVertices.size());
    mEdges.resize(rhs.mEdges.size());

    for (int i = 0; i < rhs.mVertices.size(); i++) {
      coordinate_type nc(rhs.mVertices[i]->coordinate());
      mVertices[i] = new vertex_type(nc);
      mVertices[i]->position_in_set() = i;
    }

    for (int i = 0; i < rhs.mEdges.size(); i++) {
      // std::cout << rhs.mEdges[i]->v1()->face()->position_in_set() << " "
      // 	  << rhs.mEdges[i]->v2()->face()->position_in_set() <<
      // std::endl;
      mEdges[i] = new edge_type();
      mEdges[i]->position_in_set() = i;
    }

    for (int i = 0; i < rhs.mFaces.size(); i++) {
      if (rhs.mFaces[i]) {
        face_ptr nf = new face_type();
        nf->position_in_set() = i;

        face_ptr of = rhs.mFaces[i];
        nf->size() = of->size();
        face_vertex_ptr itb = of->fbegin();
        face_vertex_ptr ite = of->fend();
        vector<face_vertex_ptr> tmpFaceArray;
        tmpFaceArray.reserve(3);
        bool iterating = true;
        int fs = 0;
        while (iterating) {
          iterating = itb != ite;
          face_vertex_ptr fv1 = new face_vertex_type();
          edge_ptr oe = itb->edge();
          edge_ptr ne = mEdges[oe->position_in_set()];
          vertex_ptr nv = mVertices[itb->vertex()->position_in_set()];

          // std::cout << i << " " << oe->position_in_set() << std::endl;

          if (oe->v1() == itb)
            ne->v1() = fv1;
          else
            ne->v2() = fv1;

          // if(!ne1->v1()) ne1->v1() = fv1; else ne1->v2() = fv1;
          fv1->vertex() = nv;
          nv->add_face_vertex(fv1);

          fv1->edge() = ne;
          fv1->face() = nf;

          tmpFaceArray.push_back(fv1);
          fs++;
          itb = itb->next();
        }
        for (int j = 0; j < tmpFaceArray.size(); j++) {
          face_vertex_ptr fvi = tmpFaceArray[j];
          face_vertex_ptr fvn = tmpFaceArray[(j + 1) % fs];
          fvi->next() = fvn;
          fvn->prev() = fvi;
        }
        nf->set_front(tmpFaceArray[0]);
        mFaces[i] = nf;
      }
    }
    mFaceRecycle = rhs.mFaceRecycle;
    mVertexRecycle = rhs.mVertexRecycle;
    mEdgeRecycle = rhs.mEdgeRecycle;
    bool here = true;
    this->update_all();
    return *this;
  }

  ~surf() {
    for (int i = 0; i < mFaces.size(); i++) {
      if (mHasFace[i]) {
        delete mFaces[i];
      }
    }

    for (int i = 0; i < mEdges.size(); i++) {
      if (mHasEdge[i]) {
        face_vertex_ptr fv0 = mEdges[i]->v1();
        face_vertex_ptr fv1 = mEdges[i]->v2();
        delete mEdges[i];
        delete fv0;
        delete fv1;
      }
    }

    for (int i = 0; i < mVertices.size(); i++) {
      if (mHasVertex[i]) {
        delete mVertices[i];
      }
    }
  }

  //		face<T>& get_face(size_t ind){
  //			return *mFaces[ind];
  //		}

  face_ptr &face(size_t ind) { return mFaces[ind]; }
  edge_ptr &edge(size_t ind) { return mEdges[ind]; }
  vertex_ptr &vertex(size_t ind) { return mVertices[ind]; }

  face_array &get_faces() { return mFaces; }
  edge_array &get_edges() { return mEdges; }
  vertex_array &get_vertices() { return mVertices; }
  size_t nverts() { return mVertices.size(); }
  size_t nedges() { return mEdges.size(); }
  size_t nfaces() { return mFaces.size(); }

  void merge(surf_ref other) {
    this->pack();
    other.pack();
    for (int i = 0; i < other.mFaces.size(); i++) {
      if (other.mFaces[i]) {
        this->push_face(other.mFaces[i]);
      }
    }
    for (int i = 0; i < other.mEdges.size(); i++) {
      if (other.mEdges[i]) {
        this->push_edge(other.mEdges[i]);
      }
    }
    for (int i = 0; i < other.mVertices.size(); i++) {
      if (other.mVertices[i]) {
        this->push_vertex(other.mVertices[i]);
      }
    }
    other.mFaces.clear();
    other.mEdges.clear();
    other.mVertices.clear();
  }

  vertex_ptr insert_vertex(bool make_face = false) {
    vertex_ptr v = new vertex_type();
    v->init();

    this->push_vertex(v);
    if (make_face) {
      face_vertex_ptr fv = v->make_face_vertex();
      face_ptr f = new face_type();
      fv->set_next(fv);
      f->set_front(fv);
      v->set_front(fv);
      this->push_face(f);
    }
    return v;
  }

  bool has_face(size_t ind) { return mHasFace[ind]; }
  bool has_edge(size_t ind) { return mHasEdge[ind]; }
  bool has_vertex(size_t ind) { return mHasVertex[ind]; }

  void push_vertex(vertex_ptr in) {
    in->position_in_set() = mVertices.size();
    mVertices.push_back(in);
    // std::cout << " ins v: " << in->position_in_set() << " " << in <<
    // std::endl;
    mHasVertex.push_back(true);
  }

  void push_edge(edge_ptr in) {
    in->position_in_set() = mEdges.size();
    mEdges.push_back(in);
    mHasEdge.push_back(true);
  }

  void push_face(face_ptr in) {
    in->position_in_set() = mFaces.size();
    mFaces.push_back(in);

    // std::cout << " ins f: " << in->position_in_set() << " " << in <<
    // std::endl;
    mHasFace.push_back(true);
  }

  void brute_force_verify(vertex_ptr v) {

    for (int i = 0; i < mEdges.size(); i++) {
      if (this->has_edge(i)) {
        edge_ptr e = mEdges[i];
        if (e->v1()->vertex() == v) {
          std::cout << " oopsy edge v1: " << std::endl;
          e->print();
        }
        if (e->v2()->vertex() == v) {
          std::cout << " oopsy edge v2: " << std::endl;
          e->print();
        }
      }
    }

    for (int i = 0; i < mFaces.size(); i++) {
      if (this->has_face(i)) {
        face_ptr f = mFaces[i];
        f->for_each([&v, &f](face_vertex_ptr fv) {
          if (fv->vertex() == v) {
            std::cout << " oopsy face: " << std::endl;
            f->print();
          }
        });
      }
    }
  }

  void remove_vertex(int i) {
    if (!mVertices[i])
      return;
    vertex_ptr v = mVertices[i];
    // std::cout << " rm v: " << i << " " << v << std::endl;
    // brute_force_verify(v);
    mVertices[i] = 0;
    mHasVertex[i] = false;
    delete v;
  }

  void remove_edge(int i) {
    edge_ptr e = mEdges[i];
    // std::cout << " rm e: " << i << " " << e << std::endl;
    mEdges[i] = 0;
    mHasEdge[i] = false;
    delete mEdges[i];
  }

  void remove_face(int i) {
    face_ptr f = mFaces[i];
    // std::cout << " rm f: " << i << " " << f << std::endl;
    mFaces[i] = 0;
    mHasFace[i] = false;
    delete f;
  }

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
    //}
    // mFaceRecycle.clear();

    // if (mVertexRecycle.size() > 0) {
    vertex_array tVertices;
    j = 0;
    for (int i = 0; i < mVertices.size(); i++) {
      if (has_vertex(i)) {
        // mVertices[i]->pack();
        tVertices.push_back(mVertices[i]);
        tVertices.back()->position_in_set() = j;
        j++;
      }
    }
    mHasVertex = std::vector<bool>(tVertices.size(), true);
    swap(mVertices, tVertices);
    //}
    // mVertexRecycle.clear();

    // if (mEdgeRecycle.size() > 0) {
    edge_array tEdges;
    j = 0;
    for (int i = 0; i < mEdges.size(); i++) {
      if (has_edge(i)) {
        tEdges.push_back(mEdges[i]);
        tEdges.back()->position_in_set() = j;
        j++;
      }
    }
    mHasEdge = std::vector<bool>(tEdges.size(), true);
    swap(mEdges, tEdges);
    //}
    // mEdgeRecycle.clear();
    // this->pack_vertices();
  }

#if 1
  void groupElements() {

    bool coloring = true;
    maxGraphColor = 0;
    for (int i = 0; i < mVertices.size(); i++) {
      if (mVertices[i])
        mVertices[i]->group() = -1;
    }
    for (int i = 0; i < mEdges.size(); i++) {
      if (mEdges[i]) {
        mEdges[i]->group() = -1;
        mEdges[i]->v1()->group() = -1;
        mEdges[i]->v2()->group() = -1;
      }
    }
    for (int i = 0; i < mFaces.size(); i++) {
      if (mFaces[i])
        mFaces[i]->group() = -1;
    }

    int currentGroup = 0;
    T r = (T)rand() / (T)RAND_MAX;
    T g = (T)rand() / (T)RAND_MAX;
    T b = (T)rand() / (T)RAND_MAX;
    mGroupColors.push_back(coordinate_type(r, g, b));

    for (int i = 0; i < mVertices.size(); i++) {
      vertex_ptr vi = mVertices[i];
      if (vi->group() == -1) {
        currentGroup++;
        T r = (T)rand() / (T)RAND_MAX;
        T g = (T)rand() / (T)RAND_MAX;
        T b = (T)rand() / (T)RAND_MAX;
        if (currentGroup > mGroupColors.size() - 1)
          mGroupColors.push_back(coordinate_type(r, g, b));
      } else
        continue;
      std::stack<int> stack;
      stack.push(i);
      while (stack.size() > 0) {
        int pId = stack.top();
        stack.pop();
        vertex_ptr v = mVertices[pId];
        v->group() = currentGroup;
        if (v && v->size() > 0) {
          face_vertex_ptr fvb = v->fbegin();
          face_vertex_ptr fve = v->fend();
          bool iterating = true;
          while (iterating) {
            iterating = fvb != fve;
            fvb->group() = currentGroup;
            fvb->edge()->group() = currentGroup;
            fvb->face()->group() = currentGroup;
            if (fvb->next()->vertex()->group() == -1) {
              stack.push(fvb->next()->vertex()->position_in_set());
            }
            // fvb->face()->color.r = mGroupColors[currentGroup][0];
            // fvb->face()->color.g = mGroupColors[currentGroup][1];
            // fvb->face()->color.b = mGroupColors[currentGroup][2];
            fvb = fvb->vnext();
          }
        }
      }
    }
  }
#endif

  void colorVertices() {
    // greedy coloring
    bool coloring = true;
    maxGraphColor = 0;
    for (int i = 0; i < mVertices.size(); i++) {
      if (mVertices[i])
        mVertices[i]->graphColor = -1;
    }

    vertex_array permVerts = mVertices;
    for (int i = 0; i < permVerts.size(); i++) {
      int card = rand() % permVerts.size();
      vertex_ptr vt = permVerts[i];
      permVerts[i] = permVerts[card];
      permVerts[card] = vt;
    }

    for (int i = 0; i < permVerts.size(); i++) {
      vertex_ptr v = permVerts[i];
      if (v && v->size() > 0) {
        int maxNeighborColor = 0;
        int minNeighborColor = 0;
        face_vertex_ptr fvb = v->fbegin();
        face_vertex_ptr fve = v->fend();
        bool iterating = true;
        while (iterating) {
          iterating = fvb != fve;
          int neighborColor = fvb->next()->vertex()->graphColor;
          maxNeighborColor = neighborColor > maxNeighborColor
                                 ? neighborColor
                                 : maxNeighborColor;
          if (neighborColor > 0)
            minNeighborColor = neighborColor < minNeighborColor
                                   ? neighborColor
                                   : minNeighborColor;
          fvb = fvb->vnext();
        }
        // std::cout << minNeighborColor << std::endl;
        // std::cout << maxNeighborColor << std::endl;
        if (minNeighborColor - 1 > -1) {
          v->graphColor = minNeighborColor - 1;
        } else {
          v->graphColor = maxNeighborColor + 1;
        }
        maxGraphColor = maxGraphColor > maxNeighborColor + 1
                            ? maxGraphColor
                            : maxNeighborColor + 1;
      }
    }
  }

  void setFaceFlags(int k) {
    for (int i = 0; i < mFaces.size(); i++) {
      if (mFaces[i])
        mFaces[i]->flag = k;
    }
  }

  void verify() {

    for (int i = 0; i < mVertices.size(); i++) {
      if (this->has_vertex(i)) {
        mVertices[i]->verify();
      }
    }

    for (int i = 0; i < mEdges.size(); i++) {
      if (this->has_edge(i)) {
        mEdges[i]->verify();
      }
    }

    for (int i = 0; i < mFaces.size(); i++) {
      if (this->has_face(i)) {
        mFaces[i]->verify();
      }
    }
  }

  void reset_flags() {
    for (int i = 0; i < mFaces.size(); i++) {
      if (mFaces[i]) {
        mFaces[i]->flag = 0;
      }
    }

    for (int i = 0; i < mEdges.size(); i++) {
      if (mEdges[i]) {
        mEdges[i]->flag = 0;
        mEdges[i]->v1()->flag = 0;
        mEdges[i]->v2()->flag = 0;

        mEdges[i]->flags[0] = 0;
        mEdges[i]->flags[1] = 0;
      }
    }

    for (int i = 0; i < mVertices.size(); i++) {
      if (mVertices[i]) {
        mVertices[i]->flag = 0;
        mVertices[i]->topologyChangeId = -1;
      }
    }
  }

  void color_dead_pointers() {
    for (int i = 0; i < mFaces.size(); i++) {
      if (mFaces[i]) {
        if (!mFaces[i]->fend() || !mFaces[i]->fbegin()) {
          mFaces[i]->color.r = 1.0;
          mFaces[i]->color.g = 0.0;
          mFaces[i]->color.b = 0.0;
          std::cout << "bad bad face" << std::endl;
        }
      }
    }

    for (int i = 0; i < mEdges.size(); i++) {
      if (mEdges[i]) {
        if (!mEdges[i]->v2()) {
          mEdges[i]->v1()->face()->color.r = 0.0;
          mEdges[i]->v1()->face()->color.g = 0.5;
          mEdges[i]->v1()->face()->color.b = 0.5;
          std::cout << "bad edge v2 pntr" << std::endl;
        }
        if (!mEdges[i]->v1()) {
          mEdges[i]->v2()->face()->color.r = 0.0;
          mEdges[i]->v2()->face()->color.g = 0.5;
          mEdges[i]->v2()->face()->color.b = 0.0;
          std::cout << "bad edge v1 pntr" << std::endl;
        }
      }
    }
    for (int i = 0; i < mVertices.size(); i++) {
      if (mVertices[i]) {
        vector<face_vertex_ptr> &fva = mVertices[i]->get_face_vertices();
        for (int j = 0; j < fva.size(); j++) {
          face_vertex_ptr &fv = fva[j];
#if 1
          if (!fv) {
            // fv->face()->color.r = 0.20;
            // fv->face()->color.g = 0.10;
            // fv->face()->color.b = 0.5;
            break;
          }
          face_vertex_ptr fve = fv->vprev();
          if (!fve) {
            // fv->face()->color.r = 0.30;
            // fv->face()->color.g = 0.40;
            // fv->face()->color.b = 0.1;
            break;
          }
          int i = 0;
          int maxIt = 100;
          while (fv != fve && i < maxIt) {
            if (!fv) {
              break;
            }
            if (!fv->prev()) {
              fv->face()->color.r = 1.0;
              fv->face()->color.g = 0.5;
              fv->face()->color.b = 0.0;
              std::cout << "bad prev pntr" << std::endl;
              break;
            }
            if (!fv->next()) {
              fv->face()->color.r = 1.0;
              fv->face()->color.g = 0.0;
              fv->face()->color.b = 0.5;
              std::cout << "bad next pntr" << std::endl;
              break;
            }
            if (!fv->vnext()) {
              fv->face()->color.r = 0.30;
              fv->face()->color.g = 0.10;
              fv->face()->color.b = 0.50;
              std::cout << "bad vnext pntr" << std::endl;
              break;
            }
            if (!fv->vprev()) {
              fv->face()->color.r = 0.6;
              fv->face()->color.g = 0.5;
              fv->face()->color.b = 0.2;
              std::cout << "bad vprev pntr" << std::endl;
              break;
            }
            if (i > maxIt / 2) {
              fv->face()->color.r = 0.1;
              fv->face()->color.g = 0.4;
              fv->face()->color.b = 0.5;
              fv->vertex()->color.r = 0.6;
              fv->vertex()->color.g = 0.4;
              fv->vertex()->color.b = 0.1;
              std::cout << "vertex exceeded maximum iterations" << std::endl;
              break;
            }
            fv = fv->vnext();
            i++;
          }

#endif
        }
      }
    }
  }

  void update_all() {
    fa_iterator fit_b = mFaces.begin();
    fa_iterator fit_e = mFaces.end();
    for (int i = 0; i < mFaces.size(); i++) {
      if (this->has_face(i)) {
        mFaces[i]->update_all();
      }
    }
    for (int i = 0; i < mVertices.size(); i++) {
      if (this->has_vertex(i)) {
        mVertices[i]->update_all();
      }
    }
  }

  void set_edge_delete_func(std::function<void(edge_ptr)> func) {
    this->mEdgeDeleteFunc = func;
  }

  void print() {
    cout << "----begin dump----" << endl;
    cout << "number of faces: " << mFaces.size() << endl;
    cout << "number of edges: " << mEdges.size() << endl;
    cout << "number of vertices: " << mVertices.size() << endl;
  }

  void print_stack() {
    cout << "----begin dump----" << endl;
    cout << "number of faces: " << mFaces.size() << endl;
    cout << "number of edges: " << mEdges.size() << endl;
    cout << "number of vertices: " << mVertices.size() << endl;
    size_t i = 0;
    // fa_iterator itb = mFaces.begin();
    // fa_iterator ite = mFaces.end();
    // while (itb != ite) {
    // 	if(*itb){
    // 	  cout << "Stack Number: "<< i << endl;
    // 	  i++;
    // 	  face_ptr cur_face = *itb;
    // 	  cout << cur_face << endl;
    // 	  cur_face->print();
    // 	}
    // 	++itb;
    // }

    va_iterator itb = mVertices.begin();
    va_iterator ite = mVertices.end();
    while (itb != ite) {
      if (*itb) {
        cout << "Vertex Number: " << i << endl;
        i++;
        cout << (*itb)->coordinate() << endl;
      }
      ++itb;
    }
    cout << "----end dump----" << endl;
  }

  void print_edge() {
    cout << "----begin dump----" << endl;
    cout << "number of faces: " << mFaces.size() << endl;
    cout << "number of edges: " << mEdges.size() << endl;
    cout << "number of vertices: " << mVertices.size() << endl;

    ea_iterator itb = mEdges.begin();
    ea_iterator ite = mEdges.end();
    while (itb != ite) {
      // cout <<"pointer addres: "<< *itb << endl;
      face_ptr cur_face = (*itb)->v1()->face();
      cout << cur_face << endl;
      cur_face->print();
      ++itb;
    }
    cout << "----end dump----" << endl;
  }

protected:
  bool manual_clean_up;
  int numGroups;

  vector<int> groupSizes;

  face_array mFaces;
  vector<bool> mHasFace;
  vector<int> mFaceRecycle;

  edge_array mEdges;
  vector<bool> mHasEdge;
  vector<int> mEdgeRecycle;

  vertex_array mVertices;
  vector<bool> mHasVertex;
  vector<int> mVertexRecycle;

  coordinate_array mGroupColors;

  std::function<void(edge_ptr)> mEdgeDeleteFunc;

public:
  int maxGraphColor;
}; // namespace m2

template <typename SPACE>
void for_each_vertex(
    typename surf<SPACE>::vertex_ptr v,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
  v->for_each(func);
}

template <typename SPACE>
void for_each_vertex_reverse(
    typename surf<SPACE>::vertex_ptr v,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
  v->for_each_reverse(func);
}

template <typename SPACE>
void for_each_face(
    typename surf<SPACE>::face_ptr f,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
  f->for_each(func);
}

template <typename SPACE>
void for_each_face_reverse(
    typename surf<SPACE>::face_ptr f,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
  f->for_each_reverse(func);
}

template <typename SPACE>
void for_each_face_except_0(
    typename surf<SPACE>::face_ptr f,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {

  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  if (!f->fbegin())
    return;
  if (!f->fend())
    return;

  face_vertex_ptr fvb = f->fbegin()->next();
  face_vertex_ptr fve = f->fend();
  bool iterating = true;
  while (iterating) {
    iterating = fvb != fve;
    func(fvb);
    fvb = fvb->next();
  }
}

template <typename SPACE>
void set_vertex_flag(typename surf<SPACE>::face_ptr f, int flag, int val) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  m2::for_each_face<SPACE>(
      f, [&](face_vertex_ptr fv) { fv->vertex()->flags[flag] = val; });
}

template <typename SPACE>
void set_vertex_flag(typename surf<SPACE>::vertex_ptr v, int flag, int val) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  m2::for_each_vertex<SPACE>(
      v, [&](face_vertex_ptr fv) { fv->vertex()->flags[flag] = val; });
}

/////////////////////////////////////////////////////////////////////////////
// flat
/////////////////////////////////////////////////////////////////////////////

template <typename SPACE, typename ITYPE> struct flattened_data_vector {
public:
  flattened_data_vector(){};
  virtual ~flattened_data_vector(){};

  virtual void push_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) = 0;

  virtual void apply_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) = 0;

  virtual void write(FILE *file) const = 0;
  virtual void read(FILE *file) = 0;
  virtual void clear() = 0;
};

template <typename SPACE, typename ITYPE, typename TYPE>
struct flattened_data_vector_t : public flattened_data_vector<SPACE, ITYPE> {
public:
  flattened_data_vector_t(){};
  virtual ~flattened_data_vector_t(){};

  virtual void push_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) {
    for (auto n : nodes) {
      TYPE data = n->template get<TYPE>(i);
      _data.push_back(data);
    }
  };

  virtual void apply_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) {
    int idx = 0;
    for (auto n : nodes) {
      TYPE data = _data[idx++];
      n->template set<TYPE>(i, data);
    }
  };

  virtual void write(FILE *file) const {
    size_t nData = _data.size();
    size_t e;
    e = fwrite((void *)&nData, sizeof(size_t), 1, file);
    e = fwrite((void *)_data.data(), sizeof(TYPE), nData, file);
  }

  virtual void read(FILE *file) {
    size_t nData;
    size_t e;
    e = fread((void *)&nData, sizeof(size_t), 1, file);
    _data.resize(nData);
    e = fread((void *)_data.data(), sizeof(TYPE), nData, file);
  }

  virtual void clear() { _data.clear(); }

  vector<TYPE> _data;
};

template <typename SPACE, typename ITYPE> struct flattened_data {
public:
  int _size = static_cast<int>(ITYPE::MAXINDEX);
  flattened_data() {}

  ~flattened_data() {
    for (auto d : _data) {
      delete d;
    }
  }

  virtual void init_data() {
    using flat_data = flattened_data_vector<SPACE, ITYPE>;
    for (int i = 0; i < static_cast<int>(ITYPE::MAXINDEX); i++) {
      ITYPE ii = static_cast<ITYPE>(i);

      storage_type type = SPACE::get_type(ii);

      flat_data *data;
      switch (type) {
      case storage_type::REAL:
        data =
            new flattened_data_vector_t<SPACE, ITYPE, typename SPACE::real>();
        break;
      case storage_type::VEC3:
        data = new flattened_data_vector_t<SPACE, ITYPE,
                                           typename SPACE::coordinate_type>();
        break;
      default:;
        // do nothing
      }
      _data.push_back(data);
    }
  }

  virtual void
  flatten(const std::vector<typename surf<SPACE>::template data_node<ITYPE> *>
              &nodes) {
    using flat_data = flattened_data_vector<SPACE, ITYPE>;
    _data = std::vector<flat_data *>();
    this->init_data();
    for (int i = 0; i < static_cast<int>(ITYPE::MAXINDEX); i++) {
      ITYPE ii = static_cast<ITYPE>(i);
      _data[i]->push_data(nodes, ii);
    }
  };

  void flatten(const std::vector<typename surf<SPACE>::face_ptr> &faces) {
    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::face_index> *>
        face_nodes;
    for (auto f : faces)
      face_nodes.push_back(f);
    std::cout << " flatten face" << std::endl;
    this->flatten(face_nodes);
  };

  void flatten(const std::vector<typename surf<SPACE>::edge_ptr> &edges) {
    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::edge_index> *>
        edge_nodes;
    for (auto e : edges)
      edge_nodes.push_back(e);
    this->flatten(edge_nodes);
  };

  void flatten(const std::vector<typename surf<SPACE>::vertex_ptr> &vertices) {

    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::vertex_index> *>
        vertex_nodes;
    for (auto v : vertices)
      vertex_nodes.push_back(v);
    this->flatten(vertex_nodes);
  };

  virtual void
  inflate(const std::vector<typename surf<SPACE>::template data_node<ITYPE> *>
              &nodes) {

    for (int i = 0; i < static_cast<int>(ITYPE::MAXINDEX); i++) {
      ITYPE ii = static_cast<ITYPE>(i);
      flattened_data_vector<SPACE, ITYPE> *data = _data[i];
      data->apply_data(nodes, ii);
    }
  };

  void inflate(const std::vector<typename surf<SPACE>::face_ptr> &faces) {
    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::face_index> *>
        face_nodes;
    for (auto f : faces)
      face_nodes.push_back(f);

    std::cout << "face nodes: " << face_nodes.size() << std::endl;

    this->inflate(face_nodes);
  };

  void inflate(const std::vector<typename surf<SPACE>::edge_ptr> &edges) {
    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::edge_index> *>
        edge_nodes;

    for (auto e : edges)
      edge_nodes.push_back(e);

    this->inflate(edge_nodes);
  };

  void inflate(const std::vector<typename surf<SPACE>::vertex_ptr> &vertices) {

    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::vertex_index> *>
        vertex_nodes;
    for (auto v : vertices)
      vertex_nodes.push_back(v);

    this->inflate(vertex_nodes);
  };

  virtual void write(FILE *file) const {
    for (auto datum : _data)
      datum->write(file);
  }

  virtual void read(FILE *file) {
    this->init_data();
    for (auto datum : _data)
      datum->read(file);
  }

  virtual void clear() {
    for (auto datum : _data)
      datum->clear();
  }

  std::vector<flattened_data_vector<SPACE, ITYPE> *> _data;
};

template <typename SPACE> struct flattened_surf {
  // class used for serialization of surfaces and associated data
public:
  M2_TYPEDEFS;
  flattened_surf() {}

  flattened_surf(const surf_ptr surf) { this->from_surf(surf); }

  ~flattened_surf() {}

  void from_surf(const surf_ptr surf) {

    nFaces = surf->get_faces().size();
    nVertices = surf->get_vertices().size();
    nCorners = 2 * surf->get_edges().size();
    corner_verts = std::vector<size_t>(2 * surf->get_edges().size(), 0);
    corner_faces = std::vector<size_t>(2 * surf->get_edges().size(), 0);
    corner_next = std::vector<size_t>(2 * surf->get_edges().size(), 0);

    std::vector<edge_ptr> edges = surf->get_edges();

    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::face_vertex_index> *>
        corners;

    auto push_fv = [&](face_vertex_ptr fv) {
      size_t c = fv->position_in_flat_set();
      vertex_ptr v = fv->vertex();
      face_ptr f = fv->face();
      face_vertex_ptr fvn = fv->next();
      size_t cn = fvn->position_in_flat_set();

      corner_next[c] = cn;
      corner_faces[c] = f->position_in_set();
      corner_verts[c] = v->position_in_set();
      corners.push_back(fv);
    };

    for (auto e : edges) {
      push_fv(e->v1());
      push_fv(e->v2());
    }

    corner_data.flatten(corners);
    edge_data.flatten(surf->get_edges());
    face_data.flatten(surf->get_faces());
    vertex_data.flatten(surf->get_vertices());
  }

  surf_ptr to_surf() {
    surf_ptr surf = new surf_type();

    vertex_array verts;
    for (int i = 0; i < nVertices; i++) {
      vertex_ptr v = new vertex_type();
      verts.push_back(v);
      surf->push_vertex(v);
    }

    face_array faces;
    for (int i = 0; i < nFaces; i++) {
      face_ptr f = new face_type();
      faces.push_back(f);
      surf->push_face(f);
    }

    std::vector<typename m2::surf<SPACE>::template data_node<
        typename SPACE::face_vertex_index> *>
        corners;

    std::vector<face_vertex_ptr> fvs;
    for (int i = 0; i < corner_verts.size(); i += 2) {
      face_vertex_ptr fv0 = new face_vertex_type();
      face_vertex_ptr fv1 = new face_vertex_type();

      edge_ptr e = new edge_type();
      e->set(fv0, fv1);
      surf->push_edge(e);
      verts[corner_verts[i]]->set_front(fv0);
      verts[corner_verts[i + 1]]->set_front(fv1);
      faces[corner_faces[i]]->set_front(fv0);
      faces[corner_faces[i + 1]]->set_front(fv1);
      fvs.push_back(fv0);
      fvs.push_back(fv1);
      corners.push_back(fv0);
      corners.push_back(fv1);
    }

    int i = 0;

    for (auto fv : fvs) {
      fv->set_next(fvs[corner_next[i++]]);
    }

    corner_data.inflate(corners);
    edge_data.inflate(surf->get_edges());
    face_data.inflate(surf->get_faces());
    vertex_data.inflate(surf->get_vertices());

    surf->update_all();
    return surf;
  }

  ///////////////////////////////////////////////////////////////////////
  // write out a field to a file stream
  ///////////////////////////////////////////////////////////////////////
  void clear() {
    nFaces = 0;
    nVertices = 0;
    corner_verts.clear();
    corner_faces.clear();
    corner_next.clear();

    corner_data.clear();
    edge_data.clear();
    face_data.clear();
    vertex_data.clear();
  }

  ///////////////////////////////////////////////////////////////////////
  // write out a field to a file stream
  ///////////////////////////////////////////////////////////////////////
  void write(FILE *file) const {
    size_t e;
    e = fwrite((void *)&nFaces, sizeof(size_t), 1, file);
    e = fwrite((void *)&nVertices, sizeof(size_t), 1, file);
    e = fwrite((void *)&nCorners, sizeof(size_t), 1, file);

    // always write out as a double

    e = fwrite((void *)corner_verts.data(), sizeof(size_t), nCorners, file);
    e = fwrite((void *)corner_faces.data(), sizeof(size_t), nCorners, file);
    e = fwrite((void *)corner_next.data(), sizeof(size_t), nCorners, file);
    corner_data.write(file);
    edge_data.write(file);
    face_data.write(file);
    vertex_data.write(file);
  }

  ///////////////////////////////////////////////////////////////////////
  // read in a field from a file stream
  ///////////////////////////////////////////////////////////////////////
  void read(FILE *file) {
    // read dimensions
    size_t e;
    e = fread((void *)&nFaces, sizeof(size_t), 1, file);
    e = fread((void *)&nVertices, sizeof(size_t), 1, file);
    e = fread((void *)&nCorners, sizeof(size_t), 1, file);
    corner_verts.resize(nCorners);
    corner_faces.resize(nCorners);
    corner_next.resize(nCorners);

    e = fread((void *)corner_verts.data(), sizeof(size_t), nCorners, file);
    e = fread((void *)corner_faces.data(), sizeof(size_t), nCorners, file);
    e = fread((void *)corner_next.data(), sizeof(size_t), nCorners, file);
    corner_data.read(file);
    edge_data.read(file);
    face_data.read(file);
    vertex_data.read(file);
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  void write(string filename) const {
    FILE *file;
    file = fopen(filename.c_str(), "wb");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " FIELD_3D write failed! " << endl;
      cout << " Could not open file " << filename.c_str() << endl;
      exit(0);
    }

    cout << " Writing file " << filename.c_str() << " ... ";
    flush(cout);

    // write to the stream
    write(file);

    // close the stream
    fclose(file);

    cout << " done. " << endl;
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  void read(string filename) {
    int size = filename.size();

    /*
    if (filename[size - 1] == 'z' && filename[size - 2] == 'g') {
      readGz(filename);
      return;
    }
    */

    FILE *file;
    file = fopen(filename.c_str(), "rb");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " FIELD_3D read failed! " << endl;
      cout << " Could not open file " << filename.c_str() << endl;
      exit(0);
    }
    // read from the stream
    read(file);

    // close the file
    fclose(file);
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  size_t nCorners = 0;
  size_t nFaces = 0;
  size_t nVertices = 0;
  std::vector<size_t> corner_verts;
  std::vector<size_t> corner_faces;
  std::vector<size_t> corner_next;

  flattened_data<SPACE, typename SPACE::face_vertex_index> corner_data;
  flattened_data<SPACE, typename SPACE::edge_index> edge_data;
  flattened_data<SPACE, typename SPACE::face_index> face_data;
  flattened_data<SPACE, typename SPACE::vertex_index> vertex_data;
};

} // namespace m2
//#undef TYPEDEF_LIST
#endif
