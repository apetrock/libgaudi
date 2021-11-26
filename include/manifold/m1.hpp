/*
 *  graph_forwards.hpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 9/26/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __M1FORWARDS__
#define __M1FORWARDS__

#include <cassert>
#include <cmath>
#include <iostream>

#define M1_TYPEDEFS                                                            \
  typedef Eigen::Matrix<T, 4, 1> coordinate_type;                              \
  typedef Eigen::Quaternion<T> quat;                                           \
  typedef m1::edge<T> edge_type;                                               \
  typedef m1::vertex<T> vertex_type;                                           \
  typedef m1::control<T> surf_type;                                            \
                                                                               \
  typedef m1::edge<T> *edge_ptr;                                               \
  typedef m1::vertex<T> *vertex_ptr;                                           \
  typedef m1::control<T> *surf_ptr;                                            \
                                                                               \
  typedef m1::edge<T> &edge_ref;                                               \
  typedef m1::vertex<T> &vertex_ref;                                           \
  typedef m1::control<T> &surf_ref;                                            \
                                                                               \
  typedef vector<coordinate_type> coordinate_array;                            \
  typedef vector<edge_ptr> edge_array;                                         \
  typedef vector<vertex_ptr> vertex_array;                                     \
  typedef typename edge_array::iterator ea_iterator;                           \
  typedef typename vertex_array::iterator va_iterator;                         \
                                                                               \
  typedef list<edge_ptr> edge_list;                                            \
  typedef list<vertex_ptr> vertex_list;

namespace m1 {
template <typename T> class vertex;

template <typename T> class edge;

template <typename T> class control;

template <typename T> class control {
  /*TODO: extend the SPACE typedef methods to this class*/
public:
  M1_TYPEDEFS
  control() {}
  /*
  control(m2::surf<T> &in) {
    in.pack();
    in.pack_vertices();
    std::vector<m2::surf<T>::vertex*> V = in.get_vertices();
    std::vector<m2::surf<T>::edge*> E = in.get_edges();
    for (long i = 0; i < V.size(); i++) {
      m2::surf<T>::vertex* vc = V[i];
      std::vector<m2::surf<T>::face_vertex*> vfv = vc->get_face_vertices();
      long sz = vfv.size();
      vertex_ptr vn = new vertex_type();
      coordinate_type c;
      c = V[i]->coordinate();
      vn->coordinate() = c;
      this->push_vertex(vn);
    }

    for (long i = 0; i < E.size(); i++) {
      edge_ptr en = new edge_type();
      long v0 = E[i]->v0()->vertex()->position_in_set();
      long v1 = E[i]->v1()->vertex()->position_in_set();
      en->v0() = mVertices[v0];
      en->v1() = mVertices[v1];
      mVertices[v0]->push_edge(en);
      mVertices[v1]->push_edge(en);
      this->push_edge(en);
    }
  }
  */
  void reset_flags() {
    for (long i = 0; i < mVertices.size(); i++) {
      if (mVertices[i])
        mVertices[i]->visited = 0;
    }
  }

  size_t size_edge() { return mEdges.size(); }
  size_t size_vertex() { return mVertices.size(); }

  bool is_sorted() const { return is_sorted_; }
  bool &is_sorted() { return is_sorted_; }

  vertex_array &get_vertices() { return mVertices; }
  edge_array &get_edges() { return mEdges; }

  vertex_array get_vertices() const { return mVertices; }
  edge_array get_edges() const { return mEdges; }

  vertex_ptr N(size_t index_) { return mVertices[index_]; }
  edge_ptr E(size_t index_) { return mEdges[index_]; }
  vertex_ptr vertex(size_t index_) { return mVertices[index_]; }
  edge_ptr edge(size_t index_) { return mEdges[index_]; }

  void update_edges() {
    mEdges.clear();
    for (size_t i = 0; i < mVertices.size(); i++) {
      this->insert_edges(*mVertices[i]);
    }
  }

  void remove_vertex(size_t vert_num) {
    vertex_ptr v = mVertices[vert_num];
    edge_array &vedges = v->get_edges();
    for (int i = 0; i < vedges.size(); i++) {
    }
  }

  void remove_edge(size_t edge_num) {
    std::cout << "removing " << edge_num << std::endl;

    edge_ptr e = mEdges[edge_num];
    if (e) {
      vertex_ptr v0 = e->v0();
      vertex_ptr v1 = e->v1();

      if (v0 == v1) {
        if (e->v0_socket() == e->v1_socket())
          v0->remove_edge(e->v0_socket());
        else {
          v0->remove_edge(e->v0_socket());
          v1->remove_edge(e->v1_socket());
        }
      } else {
        v0->remove_edge(e->this_socket(v0));
        v1->remove_edge(e->this_socket(v1));
      }
      this->delete_edge(edge_num);
    }
  }

  void collapse_edge(size_t edge_num) {

    edge_ptr e = mEdges[edge_num];
    if (e) {
      vertex_ptr v0 = e->v0();
      vertex_ptr v1 = e->v1();
      // this->remove_edge(e->position_in_set());

      coordinate_type p0 = v0->coordinate();
      coordinate_type p1 = v1->coordinate();
      vertex_ptr nv = new vertex_type;
      v0->coordinate() = 0.5 * p0 + 0.5 * p1;

      edge_array collectedEdges;
      vertex_array collectedVertices;

      edge_array &ea1 = v1->get_edges();
      bool inset = false;
      for (int i = 0; i < ea1.size(); i++) {
        edge_ptr ei = ea1[i];
        if (ei) {
          if (ei == e)
            inset = true;
          collectedEdges.push_back(ei);
          collectedVertices.push_back(ei->other(v1));
        }
      }
      assert(inset);

      for (int i = 0; i < collectedVertices.size(); i++) {
        vertex_ptr vi = collectedVertices[i];
        if (mVertices[vi->position_in_set()])
          if (vi != v0)
            this->connect(vi->position_in_set(), v0->position_in_set());
      }

      for (int i = 0; i < collectedEdges.size(); i++) {
        edge_ptr ei = collectedEdges[i];
        std::cout << ei->position_in_set() << " ";
        this->remove_edge(ei->position_in_set());
      }
      std::cout << endl;
      this->delete_vertex(v1->position_in_set());
      // this->update_placement();
    }
  }

  void merge_vertices(size_t i0, size_t i1) {
    vertex_ptr v0 = mVertices[i0];
    vertex_ptr v1 = mVertices[i1];
    if (v0 && v1) {
      // this->remove_edge(e->position_in_set());

      coordinate_type p0 = v0->coordinate();
      coordinate_type p1 = v1->coordinate();
      vertex_ptr nv = new vertex_type;
      v0->coordinate() = 0.5 * p0 + 0.5 * p1;

      edge_array collectedEdges;
      vertex_array collectedVertices;

      edge_array &ea1 = v1->get_edges();
      for (int i = 0; i < ea1.size(); i++) {
        edge_ptr ei = ea1[i];
        if (ei) {
          collectedEdges.push_back(ei);
          collectedVertices.push_back(ei->other(v1));
        }
      }

      for (int i = 0; i < collectedVertices.size(); i++) {
        vertex_ptr vi = collectedVertices[i];
        if (mVertices[vi->position_in_set()])
          if (vi != v0)
            this->connect(vi->position_in_set(), v0->position_in_set());
      }

      for (int i = 0; i < collectedEdges.size(); i++) {
        edge_ptr ei = collectedEdges[i];
        this->remove_edge(ei->position_in_set());
      }
      this->delete_vertex(v1->position_in_set());
      // this->update_placement();
    }
  }

  void delete_vertex(long i) {
    vertex_ptr v = mVertices[i];
    if (v) {
      mVertices[i] = NULL;
      mVertexRecycle.push_back(i);
    }
    // delete v;
  }

  void delete_edge(long i) {
    edge_ptr e = mEdges[i];
    std::cout << " deleting: " << i << std::endl;
    assert(e->position_in_set() == i);
    if (e) {
      mEdges[i] = NULL;
      mEdgeRecycle.push_back(i);
    }
    delete e;
  }

  void push_edge(edge_ptr edge_in) {
    mEdges.push_back(edge_in);
    edge_in->position_in_set() = mEdges.size() - 1;
    assert(mEdges[edge_in->position_in_set()] == edge_in);
  }

  void push_vertex(vertex_ptr vert_in) {
    mVertices.push_back(vert_in);
    vert_in->position_in_set() = mVertices.size() - 1;
    assert(mVertices[vert_in->position_in_set()] == vert_in);
  }

  void pack() {
    if (mVertexRecycle.size() > 0) {
      vertex_array tVertices;
      long j = 0;
      for (long i = 0; i < mVertices.size(); i++) {
        if (mVertices[i]) {
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
      long j = 0;
      for (long i = 0; i < mEdges.size(); i++) {
        if (mEdges[i]) {
          tEdges.push_back(mEdges[i]);
          tEdges.back()->position_in_set() = j;
          j++;
        }
      }
      swap(mEdges, tEdges);
    }
    mEdgeRecycle.clear();
  }

  virtual edge_ref connect(long i1, long i2) {
    vertex_ptr v0 = mVertices[i1];
    vertex_ptr v1 = mVertices[i2];
    edge_ptr ne = new edge_type;
    ne->v0() = v0;
    ne->v1() = v1;
    v0->push_edge(ne);
    v1->push_edge(ne);
    assert(ne == v0->get_edges()[ne->v0_socket()]);
    assert(ne == v1->get_edges()[ne->v1_socket()]);
    // ne->update_orientation();
    this->push_edge(ne);
    return *ne;
  }

  virtual edge_ref non_redundant_connect(long i1, long i2) {
    vertex_ptr v0 = mVertices[i1];
    vertex_ptr v1 = mVertices[i2];

    bool hasV0 = false;
    edge_array ea1 = v1->get_edges();
    for (int i = 0; i < ea1.size(); i++) {
      if (ea1[i])
        if (ea1[i]->other(v1) == v0)
          hasV0 = true;
    }

    bool hasV1 = false;
    edge_array ea0 = v0->get_edges();
    for (int i = 0; i < ea0.size(); i++) {
      if (ea0[i])
        if (ea0[i]->other(v0) == v1)
          hasV1 = true;
    }

    if (!hasV1 && !hasV0) {
      return connect(i1, i2);
    };
  }

  void update_placement() {

    typename edge_array::iterator e_itb = mEdges.begin();
    typename edge_array::iterator e_ite = mEdges.end();

    int i = 0;
    while (e_itb != e_ite) {
      if ((*e_itb) != NULL) {
        (*e_itb)->position_in_set() = i;
      }
      i++;
      e_itb++;
    }

    typename vertex_array::iterator n_itb = mVertices.begin();
    typename vertex_array::iterator n_ite = mVertices.end();

    i = 0;
    while (n_itb != n_ite) {
      if (*n_itb != 0) {
        (*n_itb)->position_in_set() = i;
      }
      i++;
      n_itb++;
    }
  }

  void clear_color() {

    typename edge_array::iterator e_itb = mEdges.begin();
    typename edge_array::iterator e_ite = mEdges.end();

    while (e_itb != e_ite) {
      if ((*e_itb) != NULL) {
        (*e_itb)->color() = 0;
      }
      e_itb++;
    }

    typename vertex_array::iterator n_itb = mVertices.begin();
    typename vertex_array::iterator n_ite = mVertices.end();

    while (n_itb != n_ite) {
      (*n_itb)->color() = 0;
      n_itb++;
    }
  }

  void remove_duplicate_edges() {
    for (long i = 0; i < mVertices.size(); i++) {

      vertex_ptr vi = mVertices[i];
      if (vi) {
        vector<edge_ptr> &eai = vi->get_edges();
        vector<edge_ptr> no_dups;
        for (long j = 0; j < vi->size_edge(); j++) {
          bool is_dup = false;
          for (long k = j + 1; k < vi->size_edge(); k++) {
            edge_ptr ej = vi->E(j);
            edge_ptr ek = vi->E(k);
            if (ej && ek) {
              if (ej->v0() == ek->v0() && ej->v1() == ek->v1()) {
                if (ej->v0() != ej->v1())
                  is_dup = true;
              }
              if (ej->v1() == ek->v0() && ej->v0() == ek->v1()) {
                if (ej->v0() != ej->v1())
                  is_dup = true;
              }
              if (is_dup) {
                this->remove_edge(ek->position_in_set());
              }
            }
          }
        }
      }
    }
  };

  void draw() {
    glDisable(GL_DEPTH_TEST); // Enables Depth Testing
    glEnable(GL_BLEND);
    glLineWidth(0.5);
    long chk = mEdges.size();
    for (long i = 0; i < mEdges.size(); i++) {
      if (mEdges[i] != NULL) {
        edge_ptr e = mEdges[i];
        coordinate_type &c1 = e->v0()->coordinate();
        coordinate_type &c2 = e->v1()->coordinate();
        glBegin(GL_LINES);
        T r = 0.5, g = 0.4, b = 0.6;

        glColor4f(r, g, b, 0.5);

        glVertex3f(c1[0], c1[1], c1[2]);
        glVertex3f(c2[0], c2[1], c2[2]);
        glEnd();
      }
    }

    for (long i = 0; i < mVertices.size(); i++) {
      vertex_ptr v0 = mVertices[i];
      if (mVertices[i] != NULL) {
        coordinate_type &c1 = mVertices[i]->coordinate();
        coordinate_type cov = v0->mCovarianceVal;
        T frac = cov[1] / (cov[0] + cov[1] + cov[2]);
        glPointSize(3.0);
        glBegin(GL_POINTS);
        //
        bool covColoring = false;
        if (covColoring) {
          if (v0->num_edges() == 2)
            glColor3f(0.9, 0.0, 0.0);
          else
            glColor3f(0.0, 0.9, 0.0);
        } else
          glColor3f(v0->drawColor[0], v0->drawColor[1], v0->drawColor[2]);
        glVertex3f(c1[0], c1[1], c1[2]);
        glEnd();

        coordinate_type &c2 = v0->mAvgLoc;
#if 0
	  glBegin(GL_POINTS);			

	  if(frac > 1e-2)
	    glColor3f(0.9, 0.0, 0.0);
	  else
	    glColor3f(0.0, 0.9, 0.0);

	  glVertex3f(c2[0],c2[1],c2[2]);
	  glEnd();
#endif
#if 1
        T denom = v0->mCovarianceVal[0] + v0->mCovarianceVal[1] +
                  v0->mCovarianceVal[2];
        for (int g = 0; g < 3; g++) {
          coordinate_type dp0 = c1;
          coordinate_type dp1 = c1 + 0.005 * v0->mCovarianceVal[g] / denom *
                                         v0->mCovarianceVec[g];

          glBegin(GL_LINES);
          glColor4f(0.0f, 0.8f, 0.5f, 0.5f);
          glVertex3f(dp0[0], dp0[1], dp0[2]);
          glVertex3f(dp1[0], dp1[1], dp1[2]);
          glEnd();
        }
#endif
      }
    }
    // glEnable(GL_DEPTH_TEST);					// Enables Depth
    // Testing glDisable(GL_BLEND);
  }

private:
  edge_array mEdges;
  vector<int> mEdgeRecycle;
  vertex_array mVertices;
  vector<int> mVertexRecycle;
  bool is_sorted_;

  void clear_edges() { mEdges.clear(); }
};

template <typename T> class edge {
public:
  M1_TYPEDEFS
  bool tagForRemoval;
  edge() { tagForRemoval = false; }

  edge(vertex_ref o_vertex, vertex_ref d_vertex, size_t o_sock, size_t d_sock) {

    if (&o_vertex == NULL || &d_vertex == NULL) {
      throw("what the fuck is going on???");
    }
    tagForRemoval = false;
    flag = 0;
    mV0 = &o_vertex;
    mV1 = &d_vertex;
    placement_in_v0 = o_sock;
    placement_in_v1 = d_sock;
    mColor = 0;
  }

  ~edge() { cout << "destroying the edge: " << this << endl; }

  virtual vertex_ptr v0() const { return mV0; }
  virtual vertex_ptr v1() const { return mV1; }
  virtual vertex_ptr &v0() { return mV0; }
  virtual vertex_ptr &v1() { return mV1; }

  bool incoming(const vertex_ptr cur_vertex) {
    bool out = false;
    if (mV1 != NULL && mV1 == cur_vertex) {
      out = true;
    }
    return out;
  }

  bool outgoing(const vertex_ptr cur_vertex) {
    bool out = false;
    if (mV0 != NULL && mV0 == cur_vertex) {
      out = true;
    }
    return out;
  }

  vertex_ptr get_other(const vertex_ptr cur_vertex) {
    vertex_ptr out = mV0;
    if (mV0 != NULL && mV0 == cur_vertex) {
      out = mV1;
    }
    return out;
  }

  vertex_ptr &other(const vertex_ptr cur_vertex) {
    if (mV0 != NULL && mV0 == cur_vertex) {
      return mV1;
    } else {
      return mV0;
    }
  }

  vertex_ptr other(vertex_ptr cur_vertex) const {
    if (mV0 != NULL && mV0 == cur_vertex) {
      return mV1;
    } else {
      return mV0;
    }
  }

  vertex_ptr &this_v(vertex_ptr cur_vertex) {
    if (mV0 != NULL && mV0 == cur_vertex) {
      return mV0;
    } else {
      return mV1;
    }
  }

  size_t other_socket(vertex_ptr cur_vertex) const {
    size_t out = placement_in_v0;
    if (mV0 != NULL && mV0 == cur_vertex) {
      out = placement_in_v1;
    }
    return out;
  }

  size_t this_socket(vertex_ptr cur_vertex) const {
    size_t out = placement_in_v1;
    if (mV0 != NULL && mV0 == cur_vertex) {
      out = placement_in_v0;
    }
    return out;
  }

  size_t stack_place() const { return placement_in_stack; }

  size_t &stack_place() { return placement_in_stack; }

  virtual int color() const { return mColor; }
  virtual int &color() { return mColor; }

  bool check_null() {
    bool is_null = false;
    if (mV0 == NULL || mV1 == NULL) {
      is_null = true;
    }
    return is_null;
  }

  // virtual al::Quat<T>& orientation()       {return mOrient;}
  // virtual al::Quat<T>  orientation() const {return mOrient;}

  // virtual void update_orientation(){
  //   al::Vec<3,T> de = mV0->coordinate() - mV1->coordinate();
  //   al::Quat<T>  q(de);
  //   q.normalize();
  //   mOrient = q;
  // }
  virtual size_t v0_socket() const { return placement_in_v0; }
  virtual size_t &v0_socket() { return placement_in_v0; }
  virtual size_t v1_socket() const { return placement_in_v1; }
  virtual size_t &v1_socket() { return placement_in_v1; }

  virtual size_t position_in_set() const { return placement_in_stack; }
  virtual size_t &position_in_set() { return placement_in_stack; }

  T k;
  long flag;

private:
  int mColor;

  quat mOrient;
  vertex_ptr mV0;
  vertex_ptr mV1;

  size_t placement_in_stack;
  size_t mID;

  size_t placement_in_v0;
  size_t placement_in_v1;
};

template <typename T> class vertex_iterator {
  M1_TYPEDEFS
public:
  vertex_iterator(vertex_ptr in) {
    cur_vertex = in;
    vert_index = 0;
  }
  vertex_ptr &operator*() {
    return cur_vertex->get_edge(vert_index)->other(cur_vertex);
  }
  vertex_ptr operator*() const {
    return cur_vertex->get_edge(vert_index)->other(cur_vertex);
  }
  vertex_iterator &operator++() {
    vert_index++;
    return *this;
  }
  vertex_iterator operator++(int) {
    vert_index++;
    return *this;
  }
  vertex_iterator &operator--() {
    vert_index--;
    return *this;
  }
  vertex_iterator operator--(int) {
    vert_index--;
    return *this;
  }
  bool operator==(vertex_iterator &other) {
    bool out = this->vert_index == other.vert_index;
    return out;
  }
  bool operator!=(vertex_iterator &other) {
    bool out = this->vert_index != other.vert_index;
    return out;
  }
  bool operator<=(vertex_iterator &other) {
    bool out = this->vert_index <= other.vert_index;
    return out;
  }
  bool operator<(vertex_iterator &other) {
    bool out = this->vert_index < other.vert_index;
    return out;
  }
  bool operator>=(vertex_iterator &other) {
    bool out = this->vert_index >= other.vert_index;
    return out;
  }
  bool operator>(vertex_iterator &other) {
    bool out = this->vert_index > other.vert_index;
    return out;
  }
  edge_ptr edge() const { return cur_vertex->get_edge(vert_index); }
  void set_end() { vert_index = cur_vertex->get_edges().size(); }
  void set_begin() { vert_index = 0; }

private:
  vertex_ptr cur_vertex;
  long vert_index;
};

template <typename T> class vertex {
public:
  M1_TYPEDEFS

  typedef std::vector<edge_ptr> edge_container;
  typedef typename std::vector<edge_ptr>::iterator edge_iterator;
  typedef std::vector<bool> bool_container;

  vertex() {
    item_in_vertex = new T[1];
    this->allocate(1);
    hasIn = 0;
    hasOut = 0;
    // id_singleton id;
  }

  vertex(T in_) {
    item_in_vertex = new T[1];
    *item_in_vertex = in_;
    this->allocate(0);
  }

  ~vertex() { mEdges.clear(); }

  void allocate(size_t cap_) {
    visited = false;
    mNumEdges = 0;
    cur_socket = 0;
    mColor = 0;
  }

  void reserve(size_t cap_) { mEdges.reserve(cap_); };

  void pack() {
    if (deallocated.size() > 0) {
      edge_array tEdges;
      long j = 0;
      for (long i = 0; i < mEdges.size(); i++) {
        if (mEdges[i]) {
          edge_ptr ei = mEdges[i];
          tEdges.push_back(ei);
          if (ei->v0() == this) {
            tEdges.back()->v0_socket() = j;
            j++;
          } else {
            tEdges.back()->v1_socket() = j;
            j++;
          }
        }
      }
      swap(mEdges, tEdges);
    }
    deallocated.clear();
  }

  vertex_iterator<T> &end() {
    vertex_iterator<T> *out = new vertex_iterator<T>(this);
    out->set_end();
    return *out;
  }
  vertex_iterator<T> &begin() {
    vertex_iterator<T> *out = new vertex_iterator<T>(this);
    out->set_begin();
    return *out;
  }

  virtual void do_stuff() {}

  virtual int color() const { return mColor; }
  virtual int &color() { return mColor; }
  T &val() { return *item_in_vertex; }
  T val() const { return *item_in_vertex; }

  virtual void set_current_socket(size_t sock_in) { cur_socket = sock_in; }

  virtual bool is_incoming(size_t in_) {
    bool out_ = !is_outgoing(in_);
    return out_;
  }

  virtual bool is_outgoing(size_t in_) {
    bool out_ = false;
    edge_iterator itb = mEdges.begin();
    int i = 0;
    while (i != in_) {
      i++;
      itb++;
    }
    if (*itb != NULL) {
      out_ = (*itb)->outgoing(this);
    }
    return out_;
  }

  virtual bool has_incoming() {
    bool checking_ = true;
    bool return_ = false;

    if (mEdges.size() > 0) {
      checking_ == true;
      edge_iterator itb = mEdges.begin();
      edge_iterator ite = mEdges.end();
      while (itb != ite) {
        edge_ptr check_vertex = (*itb);
        if (check_vertex != NULL) {
          return_ = check_vertex->incoming(this);
        }
        itb++;

        if (return_)
          itb = ite;
      }
    }
    return return_;
  }

  virtual bool has_outgoing() {
    bool checking_ = true;
    bool return_ = false;

    if (mEdges.size() > 0) {
      checking_ == true;
      edge_iterator itb = mEdges.begin();
      edge_iterator ite = mEdges.end();
      while (itb != ite) {
        edge_ptr check_vertex = (*itb);
        if (check_vertex != NULL) {
          return_ = check_vertex->outgoing(this);
        }
        itb++;

        if (return_)
          itb = ite;
      }
    }
    return return_;
  }

  /////////////////////////////////////
  // adding edges
  /////////////////////////////////////

  virtual void push_edge(edge_ptr new_edge) {
    if (deallocated.size() == 0) {
      mEdges.push_back(new_edge);
      mBEdges.push_back(true);
      mNumEdges++;
    } else {
      size_t recycle = deallocated.back();
      deallocated.pop_back();
      mEdges[recycle] = new_edge;
      mBEdges[recycle] = true;
    }
    this->update_edge_sockets();
    int thisSock = new_edge->this_socket(this);
    assert(new_edge == mEdges[thisSock]);
  }

  virtual void set_edge(edge_ptr new_edge, long i1) {
    mEdges[i1] = new_edge;
    this->update_edge_sockets();
  }

  edge_ptr get_edge(int ei) {
    if (mEdges[ei])
      return mEdges[ei];
    else {
      int eCount = 0;
      int retrieve = 0;
      for (int i = 0; i < mEdges.size(); i++) {
        if (mEdges[i])
          eCount++;
        if (eCount == ei)
          retrieve = i;
      }
      return mEdges[retrieve];
    }
  }
  virtual void push_singular_edge(edge_ptr new_edge) {
    bool is_singular = true;
    for (long i = 0; i < mEdges.size(); i++) {
      if (new_edge->v0() == mEdges[i]->v0() &&
          new_edge->v1() == mEdges[i]->v1()) {
        is_singular = false;
      }
      if (new_edge->v1() == mEdges[i]->v0() &&
          new_edge->v0() == mEdges[i]->v1()) {
        is_singular = false;
      }
    }
    if (is_singular) {
      this->push_edge(new_edge);
    }
  }

  /////////////////////////////////////
  // setting incoming and outgoing vertexs
  /////////////////////////////////////

  virtual void update_edge_sockets() {
    edge_iterator itb = mEdges.begin();
    edge_iterator ite = mEdges.end();
    for (long i = 0; i < mEdges.size(); i++) {
      edge_ptr e = mEdges[i];
      if (e) {
        if (e->incoming(this) == true) {
          e->v1_socket() = i;
          if (e->outgoing(this)) {
            e->v0_socket() = i;
          }
        } else {
          e->v0_socket() = i;
        }
      }
    }
  }

  bool has_connection(edge_ptr check) {
    bool out_ = false;
    edge_iterator itb = mEdges.begin();
    edge_iterator ite = mEdges.end();
    while (itb != ite) {
      edge_ptr check_vertex = (*itb);
      if (check->v0() == (*itb)->v0() && check->v1() == (*itb)->v1()) {
        out_ == true;
      }
      itb++;
      return true;
    }
  }

  bool has_edge(int i) const { return mBEdges[i]; }

  virtual void remove_edge(size_t edge_num) {
    if (this->num_edges() > 0) {
      edge_ptr e_ = this->get_edge(edge_num);

      size_t this_sock = e_->this_socket(this);
      std::cout << " deallocating: " << mEdges[this_sock]->position_in_set()
                << std::endl;
      mEdges[this_sock] = NULL;
      mBEdges[this_sock] = false;
      deallocated.push_back(this_sock);

      std::cout << mNumEdges << std::endl;
      mNumEdges--;
      this->update_edge_sockets();
    }
  }

  /////////////////////////////////////
  // getting connected vertexs
  /////////////////////////////////////

  virtual bool is_connected(size_t socket) const {
    size_t out = true;
    if (mEdges.size() < socket)
      out = false;
    return out;
  }

  // virtual size_t num_edges(){return  mEdges.size() - deallocated.size();};
  // virtual size_t size_edge(){return  mEdges.size() - deallocated.size();};
  virtual size_t num_edges() { return mNumEdges; };
  virtual size_t size_edge() { return mNumEdges; };
  // virtual edge_ptr& get_edge_vector(size_t index){return *mEdges;}
  virtual edge_ptr E(size_t index) { return get_edge(index); }
  virtual vertex_ref N(size_t index) { return get_connected(index); }

  virtual edge_ptr get_edge(size_t index) {
    if (this->num_edges() > 0) {
      edge_ptr e_ = mEdges[index];
      return e_;
    }
  }

  virtual edge_array &get_edges() { return mEdges; }

  virtual edge_ref get_last() {
    if (this->num_edges() > 0) {
      edge_iterator ite = mEdges.end();
      --ite;
      return **ite;
    }
  }

  virtual vertex_ref get_connected(size_t sock_num) {

    vertex_ptr adj_vertex;
    edge_ptr cur_edge = this->get_edge(sock_num);

    if (cur_edge->incoming(this)) {
      adj_vertex = cur_edge->v0();
    } else
      adj_vertex = cur_edge->v1();

    return *adj_vertex;
  }

  virtual size_t get_connected_socket(size_t sock_num) {
    if (is_connected(sock_num) == false) {
      throw("size exceeds number of edges");
    }
    //		if ( mEdges[sock_num] ==NULL) {
    //			throw("NULL SOCKET!!!");
    //		}

    edge_ptr cur_edge = this->get_edge(sock_num);
    size_t out;
    if (cur_edge->incoming(this)) {
      out = cur_edge->v0_socket();
    } else
      out = cur_edge->v1_socket();
    return out;
  }

  size_t &position_in_set() { return placement_in_stack; }
  size_t position_in_set() const { return placement_in_stack; }

  coordinate_type &coordinate() { return mCoordinate; }
  coordinate_type coordinate() const { return mCoordinate; }

  T flag;
  coordinate_type mCovarianceVal;
  coordinate_type mCovarianceVec[3];
  coordinate_type mAvgLoc;
  int clusterId;
  int isBiFurcation;
  coordinate_type drawColor;
  bool visited;
  int hasIn;
  int hasOut;

private:
  T *item_in_vertex;

  int mColor;

  size_t mID;
  size_t mNumEdges;
  size_t cur_socket;

  coordinate_type mCoordinate;
  edge_container mEdges;
  vector<bool> mBEdges;
  std::vector<size_t> deallocated;
  size_t placement_in_stack;
};

}; // namespace m1
#endif
