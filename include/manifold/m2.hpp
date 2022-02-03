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
#include <iostream>
#include <list>
#include <math.h>
#include <stack>
#include <vector>
#include <bitset>

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
#include "manifold_singleton.h"
#include "vec_addendum.h"

// this typedef list is ugly but useful!!!
#define M2_TYPEDEFS                                                            \
  typedef typename SPACE::real real;                                           \
  typedef typename SPACE::complex complex;                                     \
  typedef typename SPACE::coordinate_type coordinate_type;                     \
  typedef typename SPACE::line_type line_type;                                 \
  typedef typename SPACE::swept_point_type swept_point_type;                   \
  typedef typename SPACE::swept_triangle_type swept_triangle_type;             \
  typedef typename SPACE::box_type box_type;                                   \
  typedef typename SPACE::quat quat;                                           \
  typedef typename SPACE::mat3 mat3;                                           \
  typedef typename SPACE::mat4 mat4;                                           \
  typedef typename SPACE::double_type T;                                       \
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

template <typename SPACE> struct face_triangle;

template <typename SPACE> struct face_triangle_pair;

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
        : _ddata(other._ddata), _size(other._size), flags(other.flags) {
    }

    std::any &operator[](ITYPE i) { return _ddata[static_cast<int>(i)]; }

    const std::any &operator[](ITYPE i) const {
      return _ddata[static_cast<int>(i)];
    }

    int size() const { return _size; }

    template <typename TYPE> TYPE get(const ITYPE &i) {
      return std::any_cast<TYPE>(_ddata[static_cast<int>(i)]);
    }

    template <typename TYPE> void set(const ITYPE &i, const TYPE &d) {
      _ddata[static_cast<int>(i)] = d;
      //_ddata[static_cast<int>(i)] = datum_t<TYPE>(d);
    }

    void set_dirty(const ITYPE &i, bool val) {
      dirty[static_cast<int>(i)] = val;
    }

    void get_dirty(const ITYPE &i) { return dirty[static_cast<int>(i)]; }

    void set_dirty() {
      for (int i = 0; i < size(); i++) {
        dirty[i] = true;
      }
    }

    std::bitset<8> flags;
    int topologyChange = -1; // needs to be changed to stored value
    
  private:
    std::any _ddata[static_cast<int>(ITYPE::MAXINDEX)];
    bool dirty[static_cast<int>(ITYPE::MAXINDEX)];
    int _size = static_cast<int>(ITYPE::MAXINDEX);
  };

  class edge : public data_node<typename SPACE::edge_index> {

  public:
    edge() {
      fv1 = NULL;
      fv2 = NULL;
      mSetPosition = -1;
      m2::ID &manager = m2::ID::get_instance();
      mID = manager.new_edge_id();
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
      m2::ID &manager = m2::ID::get_instance();
      out->mID = manager.new_edge_id();
      return *out;
    }

    ~edge(){};

    int &position_in_set() { return mSetPosition; }
    int position_in_set() const { return mSetPosition; }

    int group() const { return mGroup; }
    int &group() { return mGroup; }

    size_t ID() const { return this->mID; }

    face_vertex_ptr &v1() { return fv1; }
    face_vertex_ptr v1() const { return fv1; }
    face_vertex_ptr &v2() { return fv2; }
    face_vertex_ptr v2() const { return fv2; }

    face_vertex_ptr& other(const face_vertex_ptr &cv) {
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
      nfv1->edge() = this;
      nfv2->edge() = this;
    }

    void set_this(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex) {
      if (this_vertex == fv1)
        fv1 = new_vertex;
      else
        fv2 = new_vertex;
      new_vertex->edge() = this;
    }

    void set_other(face_vertex_ptr this_vertex, face_vertex_ptr new_vertex) {
      if (this_vertex == fv1)
        fv2 = new_vertex;
      else
        fv1 = new_vertex;
      new_vertex->edge() = this;
    }

    void set_other(face_ptr cv, face_vertex_ptr ov) {
      if (cv->face() == fv1->face()) {
        return fv2 = ov;
      } else
        fv1 = ov;
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

    void verify() {
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

    face_vertex_ptr fv1;
    face_vertex_ptr fv2;

  protected:
    size_t mID;
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
      m2::ID &manager = m2::ID::get_instance();
      mID = manager.new_face_id();
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
      m2::ID &manager = m2::ID::get_instance();
      mID = manager.new_face_id();

      face_vertex_ptr nfv = new face_vertex_type();
      mSize = 1;
      mArea = 0.0;

      nfv->vertex() = &pnt;
      nfv->face() = this;

      nfv->prev() = nfv;
      nfv->next() = nfv;

      fHead = nfv;
      this->renumber_vertex_IDs();

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
      m2::ID &manager = m2::ID::get_instance();
      mID = manager.new_face_id();

      face_vertex_ptr nfv = new face_vertex_type();
      mSize = 1;
      mArea = 0.0;
      fv0->next() = fv1;
      fv1->next() = fv2;
      fv2->next() = fv0;

      fv2->prev() = fv1;
      fv1->prev() = fv0;
      fv0->prev() = fv2;

      fv0->face() = this;
      fv1->face() = this;
      fv2->face() = this;

      this->fbegin() = fv0;

      data = coordinate_type(0, 0, 0);
      data2 = coordinate_type(0, 0, 0);
      data3 = 1.0;
    }

    ~face() {
      // std::cout << "erasing face: " << ID() << std::endl;
    }

    int &size() { return mSize; }
    int size() const { return mSize; }

    int &position_in_set() { return mSetPosition; }
    int position_in_set() const { return mSetPosition; }

    // TODO: delete setHead... actually make all fronts() set_fronts()
    void setHead(face_vertex_ptr head) { fHead = head; }
    face_vertex_ptr &front() { return fHead; }

    face_vertex_ptr &fbegin() { return fHead; }

    face_vertex_ptr fbegin() const { return fHead; }

    face_vertex_ptr fend() { return fHead->prev(); }
    face_vertex_ptr fendr() { return fHead->next(); }

    face_ref operator=(const face_ref rhs) {
      face_ptr out = new face_type();
      m2::ID &manager = m2::ID::get_instance();
      out->mID = manager.new_face_id();
      out->mNormal = rhs.mNormal;
      out->mCenter = rhs.mCenter;
      out->mVertices = rhs.mVertices;

      return *out;
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
      this->update_vertex_faces();
      this->renumber_vertex_IDs();
    }

    void renumber_vertex_IDs() {
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      size_t new_id = 0;
      bool at_head = false;
      bool iterating = true;
      int i = 0;
      while (iterating && i < 200) {
        iterating = itb != ite;
        itb->position_in_face() = new_id;
        itb = itb->next();
        new_id++;
        i++;
      }
      this->mSize = new_id;
    }

    void update_vertex_faces() {
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      int i = 0;
      bool at_head = false;
      while (!at_head) {
        if (itb->vertex() == NULL)
          std::cout << "NULL VERTEX: " << std::endl;
        if (itb->next() == NULL)
          std::cout << "NULL NEXT: " << std::endl;
        if (itb->prev() == NULL)
          std::cout << "NULL PREV: " << std::endl;
        at_head = itb == ite;
        itb->face() = this;
        itb = itb->next();
        i++;
      }
    }

    void update_vertex_normals() {
      face_vertex_ptr itb = fHead;
      face_vertex_ptr ite = fHead->prev();
      int i = 0;
      bool at_head = false;
      while (!at_head && i < 200) {
        if (itb == ite) {
          at_head = true;
        }
        itb->vertex()->update_normal();
        itb = itb->next();
        i++;
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
        itb = itb->next();
      }
      return false;
    }

    int count_shared_vertices(face_ptr fB) {
      int shared_vertices = 0;

      face_vertex_ptr fvA = this->fbegin();
      face_vertex_ptr fvAe = this->fend();
      bool itA = true;
      while (itA) {
        itA = fvA != fvAe;
        vertex_ptr vA = fvA->vertex();
        face_vertex_ptr fvB = fB->fbegin();
        face_vertex_ptr fvBe = fB->fend();
        bool itB = true;
        while (itB) {
          itB = fvB != fvBe;
          vertex_ptr vB = fvB->vertex();
          shared_vertices += int(vA == vB);
          fvB = fvB->next();
        }
        fvA = fvA->next();
      }
      return shared_vertices;
    }

    int count_adjacent_shared_vertices(face_ptr fB) {
      int shared_vertices = 0;

      face_vertex_ptr fvA = this->fbegin();
      face_vertex_ptr fvAe = this->fend();

      bool itA = true;
      while (itA) {
        itA = fvA != fvAe;
        vertex_ptr vA = fvA->vertex();

        face_vertex_ptr fvvA = fvA->vnext();
        face_vertex_ptr fvvAe = fvA->vprev();
        bool itv = true;
        while (itv) {
          itv = fvvA != fvvAe;
          vertex_ptr vvA = fvvA->next()->vertex();

          face_vertex_ptr fvB = fB->fbegin();
          face_vertex_ptr fvBe = fB->fend();
          bool itB = true;
          while (itB) {
            itB = fvB != fvBe;
            vertex_ptr vB = fvB->vertex();
            shared_vertices += int(vvA == vB);
            fvB = fvB->next();
          }

          fvvA = fvvA->vnext();
        }
        fvA = fvA->next();
      }
      return shared_vertices;
    }

    void verify() {
      face_vertex_ptr itb = fbegin();
      face_vertex_ptr ite = fend();
      bool iterating = true;
      while (iterating) {
        iterating = itb != ite;

        assert(itb->prev()->next() == itb);
        assert(itb->next()->prev() == itb);
        assert(itb->face() == this);
        itb = itb->next();
      }
    }

    void print() {
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

    size_t ID() const { return this->mID; }

    int group() const { return mGroup; }
    int &group() { return mGroup; }

  protected:
    // face_vertex_array	mVertices;
    int mGroup;

    face_vertex_ptr fHead;
    coordinate_type mCenter;
    coordinate_type mNormal;
    T mArea;
    int mID;
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
    face_vertex() {
      mVertex = NULL;
      mEdge = NULL;
      mFace = NULL;

      m2::ID &manager = m2::ID::get_instance();
      mID = manager.new_face_vertex_id();
      fID = 0;
      flag = 0;
      data = 0;
    }

    face_vertex(const face_vertex_ref rhs) {
      this->mEdge = rhs.mEdge;
      //			mEdge->set_this(&rhs,this);
      this->mFace = rhs.mFace;
      this->mVertex = rhs.mVertex;
      this->nxt_face = rhs.nxt_face;
      this->prv_face = rhs.prv_face;
      this->data = rhs.data;
      fID = 0;

      m2::ID &manager = m2::ID::get_instance();
      mID = manager.new_face_vertex_id();

      flag = 0;
    }

    ~face_vertex(){
        // mVertex->remove_face_vertex(mVertexPosition);
    };

    bool operator==(const face_vertex_ref rhs) {
      if (mID == rhs.mID) {
        return true;
      } else
        return false;
    }

    face_vertex_ref operator=(const face_vertex_ref rhs) {
      face_vertex_ptr out = new face_vertex_type(rhs);
      out->mEdge->update_vertex(&rhs, this);
      return *out;
    }

    face_vertex_ptr coedge() { return mEdge->other(this); }

    face_ref get_face() { return *mFace; }

    face_vertex_ptr add_next() {
      face_vertex_ptr out = new face_vertex(*this);
      face_vertex_ptr nxt = this->nxt_face;
      out->next() = nxt;
      nxt->prev() = out;
      out->prev() = this;
      this->next() = out;

      mFace->size() += 1;

      out->face() = this->face();
      out->face()->fbegin() = out;
      this->vertex()->add_face_vertex(out);
      this->vertex()->front() = out;

      return out;
    }

    face_vertex_ptr add_prev() {
      face_vertex_ptr out = new face_vertex(*this);
      face_vertex_ptr prv = this->prv_face;
      out->prev() = prv;
      prv->next() = out;
      out->next() = this;
      this->prev() = out;

      // prv->edge()->set_this(this,out);
      mFace->size() += 1;

      out->face() = this->face();
      out->face()->fbegin() = out;
      this->vertex()->add_face_vertex(out);
      this->vertex()->front() = out;
      return out;
    }

    face_vertex_ptr &next() { return nxt_face; }
    face_vertex_ptr next() const { return nxt_face; }
    face_vertex_ptr &prev() { return prv_face; }
    face_vertex_ptr prev() const { return prv_face; }

    face_vertex_ptr vnext() {
      if (mEdge == NULL) {
        return NULL;
      } else {
        face_vertex_ptr out = mEdge->other(this);
        return out->next();
      }
    }

    face_vertex_ptr vprev() {
      if (this->prev() == NULL || this->prev()->mEdge == NULL) {
        return NULL;
      } else {
        face_vertex_ptr out = this->prev();
        if (out->mEdge)
          return out->mEdge->other(out);
        else
          return NULL;
      }
    }

    void draw(T off) {
      this->draw_vertex(off);
      this->draw_tail(off);
    }

    edge_ptr &edge() { return mEdge; }
    edge_ptr edge() const { return mEdge; }
    face_ptr &face() { return mFace; }
    face_ptr face() const { return mFace; }
    face_ptr coface() const { return mEdge->other(this)->face(); }
    vertex_ptr &vertex() { return mVertex; }
    vertex_ptr vertex() const { return mVertex; }

    int vertex_ID() const { return mVertex->ID(); }
    int ID() const { return this->mID; }
    int face_ID() const { return this->mFacePosition; }
    int &face_ID() { return this->mFacePosition; }

    int &position_in_face() { return mFacePosition; }
    int position_in_face() const { return mFacePosition; }
    void set_edge(edge_ref input) { mEdge = &input; };
    void set_face(face_ref input) { mFace = &input; };
    void set_vertex(vertex_ref input) { mVertex = &input; };

    int &group() { return mGroup; }
    int group() const { return mGroup; }

  protected:
    int mGroup;
    int mID;
    int fID;

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

      m2::ID &manager = m2::ID::get_instance();
      mID = manager.new_vertex_id();

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
    }

    int &position_in_set() { return mSetPosition; }

    int &group() { return mGroup; }
    int group() const { return mGroup; }

    int &ID() { return mSetPosition; }
    int ID() const { return mSetPosition; }

    int calc_size() {
      face_vertex_ptr fvb = this->fbegin();
      face_vertex_ptr fve = this->fend();
      bool iterating = true;
      int sz = 0;
      while (iterating) {
        iterating = fvb != fve;
        fvb = fvb->vnext();
        sz++;
      }
      return sz;
    }

    int size() { return mSize; }

    face_vertex_ptr make_face_vertex() {
      face_vertex_ptr fv = new face_vertex_type();
      mFront = fv;
      fv->vertex() = this;
      mSize++;
      return fv;
    }

    void add_face_vertex(face_vertex_ptr new_fv) {
      mFront = new_fv;
      mSize++;
      // mFaceVertices.push_back(new_fv);
      // fvl_iterator end = mFaceVertices.end();
      // end--;
      // new_fv->position_in_vertex() = end;
      new_fv->vertex() = this;
    }

    void remove_face_vertex(face_vertex_ptr fv) {
      // mFront = fv->vnext();
      mSize--;
      // mFaceVertices.erase(it);

      // if (mFaceVertices[in]) {
      // 	face_vertex_ptr fv = mFaceVertices[in];
      // 	mFaceVertices[in] = NULL;
      // 	mRecycle.push_back(in);
      // }
    }

    void update_all() {

      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;

      while (iterating) {
        mSize++;
        iterating = itb != ite;
        itb->vertex() = this;
        itb = itb->vnext();
      }
    }
    // void pack(){
    //   if (mRecycle.size() > 0) {
    // 	vector<face_vertex_ptr> tFaceVertices;
    // 	int j = 0;
    // 	for (int i = 0; i < mFaceVertices.size(); i++) {
    // 	  if(mFaceVertices[i]){
    // 	    tFaceVertices.push_back(mFaceVertices[i]);
    // 	    tFaceVertices.back()->position_in_vertex() = j;
    // 	    j++;
    // 	  }
    // 	}
    // 	swap(mFaceVertices,tFaceVertices);
    //   mRecycle.clear();
    //   }
    // }

    face_vertex_ptr &front() { return mFront; }
    face_vertex_ptr fbegin() { return this->front(); }
    face_vertex_ptr fend() { return this->front()->vprev(); }

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
      if (mSize == 0)
        return;
      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      while (iterating) {
        iterating = itb != ite;
        if (itb->vertex() != this)
          std::cout << " asserted: " << position_in_set() << std::endl;
        assert(itb->vertex() == this);
        itb = itb->vnext();
      }
    }

    void print() {
      if (mSize == 0)
        return;
      face_vertex_ptr itb = this->fbegin();
      face_vertex_ptr ite = this->fend();
      bool iterating = true;
      std::cout << " - vertex: " << mSetPosition << ", size: " << this->mSize
                << std::endl;
      std::cout << " edges: ";
      while (iterating) {
        iterating = itb != ite;
        // std::cout << itb->next()->vertex()->position_in_set() << ", ";
        std::cout << itb->edge()->position_in_set() << ", ";
        itb = itb->vnext();
      }
      std::cout << std::endl;
    }

  protected:
    face_vertex_ptr mFront = NULL;

    size_t mID;
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
    // mEdges.clear();
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
        nf->fbegin() = fv0;
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
        nf->fbegin() = tmpFaceArray[0];
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

  ~surf() {}

  //		face<T>& get_face(size_t ind){
  //			return *mFaces[ind];
  //		}

  face_ptr &face(size_t ind) { return mFaces[ind]; }
  edge_ptr &edge(size_t ind) { return mEdges[ind]; }
  vertex_ptr &vertex(size_t ind) { return mVertices[ind]; }

  bool has_face(size_t ind) {
    if (mFaces[ind])
      return true;
    else
      return false;
  }

  bool has_edge(size_t ind) {
    if (mEdges[ind])
      return true;
    else
      return false;
  }

  bool has_vertex(size_t ind) {
    if (mVertices[ind])
      return true;
    else
      return false;
  }

  face_array &get_faces() { return mFaces; }
  edge_array &get_edges() { return mEdges; }
  vertex_array &get_vertices() { return mVertices; }

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

  vertex_ptr insert_vertex() {
    vertex_ptr new_vert = new vertex_type();
    new_vert->init();
    face_ptr new_face = new_vert->front()->face();
    this->push_face(new_face);
    this->push_vertex(new_vert);
    return new_vert;
  }

  void push_vertex(vertex_ptr in) {
    if (mVertexRecycle.size() > 0 && !manual_clean_up) {
      int i = mVertexRecycle.back();
      in->position_in_set() = i;
      mVertexRecycle.pop_back();
      mVertices[i] = in;
    } else {
      mVertices.push_back(in);
      in->position_in_set() = mVertices.size() - 1;
    }
  }

  void push_edge(edge_ptr in) {
    if (mEdgeRecycle.size() > 0 && !manual_clean_up) {
      int i = mEdgeRecycle.back();
      in->position_in_set() = i;
      mEdgeRecycle.pop_back();
      mEdges[i] = in;
    } else {
      in->position_in_set() = mEdges.size();
      mEdges.push_back(in);
    }
  }

  void push_face(face_ptr in) {
    if (mFaceRecycle.size() > 0 && manual_clean_up) {
      int i = mFaceRecycle.back();
      in->position_in_set() = i;
      mFaceRecycle.pop_back();
      mFaces[i] = in;
    } else {
      mFaces.push_back(in);
      in->position_in_set() = mFaces.size() - 1;
    }
  }

  void remove_vertex(int i) {
    vertex_ptr v = mVertices[i];
    // std::cout << "       deleting vert: " << i << std::endl;
    // std::cout << " removing vert: " << i << " " << v  << " " << v->size() <<
    // std::endl;
    mVertexRecycle.push_back(i);
    if (!manual_clean_up) {
      // if(v->size() > 1) throw("fuck you fucktard");
      mVertices[i] = NULL;
      delete v;
    }
  }

  void remove_edge(int i) {
    edge_ptr e = mEdges[i];
    // std::cout << "       deleting edge: " << i << std::endl;
    // if(i == 175548) throw('gdb here we come!');
    if (this->mEdgeDeleteFunc)
      this->mEdgeDeleteFunc(e);
    mEdgeRecycle.push_back(i);

    if (!manual_clean_up) {
      delete mEdges[i];
      mEdges[i] = NULL;
    }
  }

  void remove_face(int i) {
    // std::cout << "       deleting face: " << i << std::endl;
    face_ptr f = mFaces[i];
    mFaceRecycle.push_back(i);
    if (!manual_clean_up) {
      mFaces[i] = NULL;
      delete f;
    };
  }

  void toggle_clean_up() { manual_clean_up ^= true; }

  void clean_up() {
    // cleanup globally deletes pointers after an operation that needs them to
    // exist.
    for (int i = 0; i < mFaceRecycle.size(); i++) {
      int ii = mFaceRecycle[i];
      delete mFaces[ii];
      mFaces[ii] = NULL;
    }
    for (int i = 0; i < mEdgeRecycle.size(); i++) {
      int ii = mEdgeRecycle[i];
      delete mEdges[ii];
      mEdges[ii] = NULL;
    }
    for (int i = 0; i < mVertexRecycle.size(); i++) {
      int ii = mVertexRecycle[i];
      delete mVertices[ii];
      mVertices[ii] = NULL;
    }
    this->toggle_clean_up();
  }

  void pack() {
    // TODO: safer pack, iterating from mRecycle[i] to mRecycle[i+1]
    if (mFaceRecycle.size() > 0) {
      face_array tFaces;
      int j = 0;
      for (int i = 0; i < mFaces.size(); i++) {
        if (mFaces[i]) {
          tFaces.push_back(mFaces[i]);
          tFaces.back()->position_in_set() = j;
          j++;
        }
      }
      swap(mFaces, tFaces);
    }
    mFaceRecycle.clear();

    if (mVertexRecycle.size() > 0) {
      vertex_array tVertices;
      int j = 0;
      for (int i = 0; i < mVertices.size(); i++) {
        if (mVertices[i]) {
          // mVertices[i]->pack();
          tVertices.push_back(mVertices[i]);
          tVertices.back()->position_in_set() = j;
          j++;
        }
      }
      swap(mVertices, tVertices);
    }
    mVertexRecycle.clear();

    if (mEdgeRecycle.size() > 0) {
      edge_array tEdges;
      int j = 0;
      for (int i = 0; i < mEdges.size(); i++) {
        if (mEdges[i]) {
          tEdges.push_back(mEdges[i]);
          tEdges.back()->position_in_set() = j;
          j++;
        }
      }
      swap(mEdges, tEdges);
    }
    mEdgeRecycle.clear();
    this->pack_vertices();
  }

  void pack_vertices() {
    for (int i = 0; i < this->mVertices.size(); i++) {
      // mVertices[i]->pack();
    }
  }

  void validateFaceVertices() {
    for (int i = 0; i < mVertices.size(); i++) {
      vertex_ptr v = mVertices[i];
      if (v) {
        face_vertex_ptr fvb = v->fbegin();
        face_vertex_ptr fve = v->fend();
        bool iterating = true;
        while (iterating) {
          if (fvb->vertex() != v) {
            std::cout << "face vertex: " << fvb
                      << " has improperly assigned vertex!" << std::endl;
          }
          iterating = fvb != fve;
          fvb = fvb->vnext();
        }
      }
    }

    for (int i = 0; i < mFaces.size(); i++) {
      face_ptr f = mFaces[i];
      if (f) {
        face_vertex_ptr fvb = f->fbegin();
        face_vertex_ptr fve = f->fend();
        bool iterating = true;
        while (iterating) {
          if (fvb->face() != f) {
            std::cout << "face vertex: " << fvb
                      << " has improperly assigned face!" << std::endl;
          }
          iterating = fvb != fve;
          fvb = fvb->next();
        }
      }
    }
    for (int i = 0; i < mEdges.size(); i++) {
      edge_ptr e = mEdges[i];
      if (e) {
        if (e->v1()->edge() != e) {
          std::cout << "face vertex: " << e->v1()->edge()
                    << " has improperly assigned edge!" << std::endl;
        }
        if (e->v2()->edge() != e) {
          std::cout << "face vertex: " << e->v2()->edge()
                    << " has improperly assigned edge!" << std::endl;
        }
      }
    }
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
    for (int i = 0; i < mFaces.size(); i++) {
      if (mFaces[i]) {
        mFaces[i]->verify();
      }
    }

    for (int i = 0; i < mEdges.size(); i++) {
      if (mEdges[i]) {
        mEdges[i]->verify();
      }
    }

    for (int i = 0; i < mVertices.size(); i++) {
      if (mVertices[i]) {
        mVertices[i]->verify();
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
        mVertices[i]->topologyChange = -1;
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
      if (mFaces[i]) {
        mFaces[i]->update_all();
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
  vector<int> mFaceRecycle;
  edge_array mEdges;
  vector<int> mEdgeRecycle;
  vertex_array mVertices;
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
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  if (!v->fbegin())
    return;
  if (!v->fend())
    return;

  face_vertex_ptr fvb = v->fbegin();
  face_vertex_ptr fve = v->fend();
  bool iterating = true;
  while (iterating) {
    iterating = fvb != fve;
    func(fvb);
    fvb = fvb->vnext();
  }
}

template <typename SPACE>
void for_each_vertex_reverse(
    typename surf<SPACE>::vertex_ptr v,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  if (!v->fbegin())
    return;

  if (!v->fendr())
    return;

  face_vertex_ptr fvb = v->fbegin();
  face_vertex_ptr fve = v->fendr();
  bool iterating = true;
  while (iterating) {
    iterating = fvb != fve;
    func(fvb);
    fvb = fvb->vprev();
  }
}

template <typename SPACE>
void for_each_face(
    typename surf<SPACE>::face_ptr f,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {

  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  if (!f->fbegin())
    return;
  if (!f->fend())
    return;

  face_vertex_ptr fvb = f->fbegin();
  face_vertex_ptr fve = f->fend();
  bool iterating = true;

  while (iterating) {

    iterating = fvb != fve;
    func(fvb);
    fvb = fvb->next();

  }

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
void for_each_face_reverse(
    typename surf<SPACE>::face_ptr f,
    std::function<void(typename surf<SPACE>::face_vertex_ptr fv)> func) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  if (!f->fbegin())
    return;

  if (!f->fendr())
    return;

  face_vertex_ptr fvb = f->fbegin();
  face_vertex_ptr fve = f->fendr();
  bool iterating = true;
  while (iterating) {
    iterating = fvb != fve;
    func(fvb);
    fvb = fvb->prev();
  }
}

} // namespace m2
//#undef TYPEDEF_LIST
#endif
