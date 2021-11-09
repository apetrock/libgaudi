#ifndef __M2MOVING__
#define __M2MOVING__

#include <stack>

//#ifdef _OPENMP
//# include <omp.h>
//#endif

#include "bins.hpp"
#include "conj_grad.hpp"
#include "debugger.h"
#include "geometry_types.hpp"
#include "m2Includes.h"
#include "modify.hpp"
#include "remesh.hpp"

#include "tree_code.hpp"

#include "TIMER.h"

#include "surface_calculator.hpp"
#include "surface_filters.hpp"
#include <cmath>

namespace m2 {

template <typename SPACE> class set_operations {
  M2_TYPEDEFS;

public:
  static bool flip_edges(control_ptr obj) {
    // TIMER function//TIMER(__FUNCTION__);
    edge_array &edges = obj->get_edges();
    m2::construct<SPACE> cons;
    bool flipped = false;
    edge_array permEdges;
    for (int i = 0; i < edges.size(); i++) {
      if (!obj->has_edge(i))
        continue;
      edge_ptr e = edges[i];
      if (e->v1()->face()->size() == 3 && e->v2()->face()->size() == 3)
        permEdges.push_back(edges[i]);
    }

    for (int i = 0; i < permEdges.size(); i++) {
      int card = rand() % permEdges.size();
      edge_ptr et = permEdges[i];
      permEdges[i] = permEdges[card];
      permEdges[card] = et;
    }

    for (int i = 0; i < permEdges.size(); i++) {
      edge_ptr e = permEdges[i];
      bool pinned = e->v1()->vertex()->pinned == true &&
                    e->v2()->vertex()->pinned == true;
      bool notPinned = e->v1()->vertex()->pinned != true &&
                       e->v2()->vertex()->pinned != true;
      if (pinned)
        continue;
      if (notPinned)
        continue;

      // if(e->v1()->vertex()->size() < 4) continue;
      // if(e->v2()->vertex()->size() < 4) continue;

      face_vertex_ptr v0 = e->v1();
      face_vertex_ptr v1 = v0->prev();
      face_vertex_ptr v2 = e->v2();
      face_vertex_ptr v3 = v2->prev();

      coordinate_type c0 = v0->coordinate();
      coordinate_type c1 = v1->coordinate();
      coordinate_type c2 = v2->coordinate();
      coordinate_type c3 = v3->coordinate();

      face_vertex_ptr fvEdge = NULL;
      T m01 = 1.0 / (c0 - c1).norm();
      T m12 = 1.0 / (c1 - c2).norm();
      T m23 = 1.0 / (c2 - c3).norm();
      T m30 = 1.0 / (c3 - c0).norm();

      T cos0 = (c1 - c0).dot(c3 - c0) * m01 * m30;
      T cos1 = (c0 - c1).dot(c2 - c1) * m01 * m12;
      T cos2 = (c1 - c2).dot(c3 - c2) * m12 * m23;
      T cos3 = (c0 - c3).dot(c2 - c3) * m30 * m23;
      // half angle cos^2(2a) = 0.5*(1+cos(a))
      T cSame = acos(cos1) + acos(cos3); // corresponds to flipped edge
      T cFlip = acos(cos0) + acos(cos2); // corresponds to flipped edge
      // T cSame =  (M_PI - (acos(cos1) + acos(cos3))); //corresponds to flipped
      // edge

      T eFlip = cFlip * cFlip;
      T eSame = cSame * cSame;
#if 0
	T sin0 = va::cross(c1-c0,c3-c0).norm();
	T sin1 = va::cross(c0-c1,c2-c1).norm();
	T sin2 = va::cross(c1-c2,c3-c2).norm();
	T sin3 = va::cross(c0-c3,c2-c3).norm();
	bool div0 = (sin0 < 1e-12 || sin1 < 1e-12 || sin2 < 1e-12 || sin3 < 1e-12);
	if(!div0){
	  //curvature penalty
	  T cot0 = cos0/sin0;
	  T cot1 = cos1/sin1;
	  T cot2 = cos2/sin2;
	  T cot3 = cos3/sin3; 
	  //if we flip we change this much	  

	  T eCurveFlip = (cot1 + cot3)*(cot1 + cot3) + (cot0 + cot2)*(cot0 + cot2);
	  T eCurveSame = 0.0;

	  T C = 0.00000001;
	  eFlip += C*eCurveFlip;
	  eSame += C*eCurveSame;
	}
#endif

      // if(cSame > M_PI){
      if (cFlip < cSame) {
        m2::construct<SPACE> cons;
        cons.flip_edge(obj, e);
        // e->v1()->data *= -1;
        // e->v2()->data *= -1;
        flipped = true;
      }
    }
    return flipped;
  }

}; // set_operations

template <typename SPACE> class join_mesh {
  M2_TYPEDEFS;

public:
  //---------------------------------------------------------------------------
  // local structures
  //---------------------------------------------------------------------------
  struct contact_manifold {
  public:
    int p1, p2;
    int p1type, p2type; // means that we can querie the collision type
    control_ptr mMesh;
    T length() {
      if (p1type == p2type) {
        vector<vertex_ptr> &verts = this->mMesh->get_vertices();
        T d0 = norm2(verts[p1]->coordinate() - verts[p2]->coordinate());
        return d0;
      }
    }
  };

  struct contact_manifold_sort {
  public:
    bool operator()(contact_manifold c0, contact_manifold c1) {
      // point 0
      // edge  1
      // tri   2
      T d0, d1;
      if (c0.p1type == c0.p2type) {
        vector<vertex_ptr> &verts = this->mMesh->get_vertices();
        d0 = norm2(verts[c0.p1]->coordinate() - verts[c0.p2]->coordinate());
      }
      if (c1.p1type == c1.p2type) {
        vector<vertex_ptr> &verts = this->mMesh->get_vertices();
        d1 = norm2(verts[c1.p1]->coordinate() - verts[c1.p2]->coordinate());
      }

      return (d0 < d1);
    }

    control_ptr mMesh;
  } mContactSorter;

  // in order to run unit test
  join_mesh(){};

  join_mesh(control_ptr mesh, T tolerance) {
    tol = tolerance;
    mMesh = mesh;
  };

  //---------------------------------------------------------------------------
  // constructor and init
  //---------------------------------------------------------------------------

  void init() {
    // initializes the trees necessary to make this all work
    vector<face_ptr> &faces = mMesh->get_faces();
    for (int i = 0; i < faces.size(); i++) {
      triangle_type tempFace =
          makeTriangle(faces[i]->fbegin(), faces[i]->fbegin()->next(),
                       faces[i]->fbegin()->prev());
      oldFaces.push_back(tempFace);
    }
    exteriorTree.build(oldFaces, 10);
  }

  inline void dimPolynomial(coordinate_type a, coordinate_type va,
                            coordinate_type b, coordinate_type vb,
                            coordinate_type c, coordinate_type vc, int i, int j,
                            int k, T &o, T &p, T &q, T &r) {
    o = va[i] * vb[j] * vc[k];
    p = a[i] * vb[j] * vc[k] + va[i] * b[j] * vc[k] + va[i] * vb[j] * c[k];
    q = va[i] * b[j] * c[k] + a[i] * vb[j] * c[k] + a[i] * b[j] * vc[k];
    r = a[i] * b[j] * c[k];
    // std::cout << o << " " << p << " " << q << " " << r << std::endl;
  }

  void buildPolynomial(coordinate_type a, coordinate_type va, coordinate_type b,
                       coordinate_type vb, coordinate_type c,
                       coordinate_type vc, T &p, T &q, T &r) {
    T o0, o1, o2, o3, o4, o5;
    T p0, p1, p2, p3, p4, p5;
    T q0, q1, q2, q3, q4, q5;
    T r0, r1, r2, r3, r4, r5;
    dimPolynomial(a, va, b, vb, c, vc, 1, 2, 0, o0, p0, q0, r0);
    dimPolynomial(a, va, b, vb, c, vc, 2, 1, 0, o1, p1, q1, r1);
    dimPolynomial(a, va, b, vb, c, vc, 2, 0, 1, o2, p2, q2, r2);
    dimPolynomial(a, va, b, vb, c, vc, 0, 2, 1, o3, p3, q3, r3);
    dimPolynomial(a, va, b, vb, c, vc, 0, 1, 2, o4, p4, q4, r4);
    dimPolynomial(a, va, b, vb, c, vc, 1, 0, 2, o5, p5, q5, r5);
    T o = o0 - o1 + o2 - o3 + o4 - o5;
    p = p0 - p1 + p2 - p3 + p4 - p5;
    q = q0 - q1 + q2 - q3 + q4 - q5;
    r = r0 - r1 + r2 - r3 + r4 - r5;
    // std::cout << o << " " << p << " " << q << " " << r << std::endl;
    // o = o < 1e-12 ? 1e-12 : o;
    o += 1e-12;
    p /= o;
    q /= o;
    r /= o;
    // std::cout << o << " " << p << " " << q << " " << r << std::endl;
  }

  T fRand() {
    T f = (T)rand() / (T)RAND_MAX;
    return f;
  }

  T qRand() {
    T f = (T)rand() / (T)RAND_MAX;
    return 1.0 - 2.0 * f;
  }

  void testCollision() {

    // we'll build a random triangle, and a random point from the interior of
    // the triangle

    coordinate_type T0(qRand(), qRand(), qRand());
    coordinate_type T1(qRand(), qRand(), qRand());
    coordinate_type T2(qRand(), qRand(), qRand());
    T q0 = 0.33 * fRand();
    T q1 = 0.33 * fRand();
    T q2 = 1.0 - q0 - q1;

    // coordinate_type T0(-1.0, 0.5, 0.2);
    // coordinate_type T1( 0.5, 1.0, 0.2);
    // coordinate_type T2( 1.0,-1.0,-0.2);

    // T q0 = 0.33;
    // T q1 = 0.33;
    // T q2 = 1.0 - q0 - q1;

    coordinate_type x = q0 * T0 + q1 * T1 + q2 * T2;
    coordinate_type tri[3] = {T0, T1, T2};

    m2::distance_calculator<SPACE> calc;
    T dist = calc.distanceFromTriangle(tri, x);
    std::cout << " dist from tri: " << dist << std::endl;

    coordinate_type v0(qRand(), qRand(), qRand());
    coordinate_type v1(qRand(), qRand(), qRand());
    coordinate_type v2(qRand(), qRand(), qRand());
    coordinate_type vx(qRand(), qRand(), qRand());

    // coordinate_type v0(-0.1,-0.22,-0.1);
    // coordinate_type v1(-0.2,-0.1,-0.11);
    // coordinate_type v2(-0.1,-0.2, 0.3);
    // coordinate_type vx( 0.1, 0.1, 0.2);

    // choose a random time, displace the triangle and point a random vector
    // away
    T t = fRand();
    coordinate_type Tp0 = T0 + v0 * t;
    coordinate_type Tp1 = T1 + v1 * t;
    coordinate_type Tp2 = T2 + v2 * t;
    coordinate_type xp = x + vx * t;

    // now run the continuous collision detection backwards and see if we get a
    // negative time as one of the roots

    T p, q, r, roots[3];
    buildPolynomial(Tp1 - Tp0, -v1 + v0, Tp2 - Tp0, -v2 + v0, xp - Tp0,
                    -vx + v0, p, q, r);
    int rootcount =
        magnet::math::cubicSolve(p, q, r, roots[0], roots[1], roots[2]);
    if (rootcount == 0)
      std::cout << " no roots  " << std::endl;
    if (rootcount == 1)
      std::cout << " poly   : " << p << " " << q << " " << r << " t: " << t
                << " roots: " << roots[0] << std::endl;
    if (rootcount == 2)
      std::cout << " poly   : " << p << " " << q << " " << r << " t: " << t
                << " roots: " << roots[0] << " " << roots[1] << std::endl;
    if (rootcount == 3)
      std::cout << " poly   : " << p << " " << q << " " << r << " t: " << t
                << " roots: " << roots[0] << " " << roots[1] << " " << roots[2]
                << std::endl;
  }

  void firstImpulse(vector<coordinate_type> &vels, T eps) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<contact_manifold> collisions;

    vector<vertex_ptr> &vertices = mMesh->get_vertices();
    vector<coordinate_type> points;
    points.resize(vertices.size());
    for (int i = 0; i < vertices.size(); i++) {
      points[i] = vertices[i]->coordinate();
    }
    vector<face_ptr> &faces = mMesh->get_faces();
    vector<int> triIndex;
    triIndex.resize(3 * faces.size());
    vector<triangle_type> tris;
    tris.resize(faces.size());

    for (int i = 0; i < faces.size(); i++) {
      tris[i][0] = faces[i]->fbegin()->coordinate();         // 0
      tris[i][1] = faces[i]->fbegin()->next()->coordinate(); // 1
      tris[i][2] =
          faces[i]->fbegin()->prev()->coordinate(); // 2 {next->next == prev}

      triIndex[3 * i + 0] = faces[i]->fbegin()->vertex()->position_in_set();
      triIndex[3 * i + 1] =
          faces[i]->fbegin()->next()->vertex()->position_in_set();
      triIndex[3 * i + 2] =
          faces[i]->fbegin()->prev()->vertex()->position_in_set();
    }

    typedef aabb_tree<SPACE, triangle_type> triangle_tree;
    typedef typename triangle_tree::node_type edge_node_type;
    triangle_tree tree(tris);
    for (int i = 0; i < points.size(); i++) {
      vector<int> collectedTris;
      vector<int> filteredTris;
      coordinate_type ci = points[i];
      swept_point_type li(ci, vels[i]);
      li.dt = dt;
      coordinate_type veli = vels[i];
      int indexMin = 0;
      coordinate_type cmin;
      T d = this->minLength;
      getAllNearest<SPACE, triangle_type, coordinate_type>(
          ci, tree, tris, collectedTris, 2.0 * d);
      coordinate_type avgNormal(0.0, 0.0, 0.0);
      coordinate_type avgVelocity(0.0, 0.0, 0.0);
      T accumWeight = 0;

      for (int j = 0; j < collectedTris.size(); j++) {
        int pj = tree.permutation[j];
        int j0 = triIndex[3 * collectedTris[j] + 0];
        int j1 = triIndex[3 * collectedTris[j] + 1];
        int j2 = triIndex[3 * collectedTris[j] + 2];
        if (j0 == i)
          continue;
        if (j1 == i)
          continue;
        if (j2 == i)
          continue;
        vertex_ptr v4 = vertices[i];
        vertex_ptr vj[3] = {vertices[j0], vertices[j1], vertices[j2]};
        coordinate_type velj[3] = {vels[j0], vels[j1], vels[j2]};

        for (int k = 0; k < 3; k++) {
          int jk = triIndex[3 * collectedTris[j] + k];
          if (v4->shares_edge_with(vj[k]))
            continue;
          if (jk == i)
            continue;

          coordinate_type cj = vertices[jk]->coordinate();

          T dist = norm2(ci - cj);

          T w = 1.0 / powf(dist * dist + d * d, 1.0);
          avgNormal += w * vj[k]->normal();
          avgVelocity += w * velj[k];
          accumWeight += w;
        }
      }

      if (accumWeight > 1e-16) {
        avgNormal /= accumWeight;
        avgVelocity /= accumWeight;

        T J = (eps - d) / dt - va::dot(avgNormal, avgVelocity - veli);
        std::cout << accumWeight << " " << avgVelocity << " " << dt << " " << J
                  << std::endl;
        vels[i] += J * avgNormal;
      }
    }
  }

  vector<contact_manifold> getCollisions(const vector<coordinate_type> &vels) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &vertices = mMesh->get_vertices();
    vector<coordinate_type> points;
    points.resize(vertices.size());
    vector<contact_manifold> collisions;
    for (int i = 0; i < vertices.size(); i++) {
      points[i] = vertices[i]->coordinate();
    }

    vector<face_ptr> &faces = mMesh->get_faces();
    vector<int> triIndex;
    triIndex.resize(3 * faces.size());
    vector<swept_triangle_type> tris;
    tris.resize(faces.size());

    for (int i = 0; i < faces.size(); i++) {
      tris[i][0] = faces[i]->fbegin()->coordinate();         // 0
      tris[i][1] = faces[i]->fbegin()->next()->coordinate(); // 1
      tris[i][2] =
          faces[i]->fbegin()->prev()->coordinate(); // 2 {next->next == prev}

      tris[i].v[0] = vels[faces[i]->fbegin()->vertex()->position_in_set()]; // 0
      tris[i].v[1] =
          vels[faces[i]->fbegin()->next()->vertex()->position_in_set()]; // 1
      tris[i].v[2] = vels[faces[i]
                              ->fbegin()
                              ->prev()
                              ->vertex()
                              ->position_in_set()]; // 2
                                                    // {next->next
                                                    // ==
                                                    // prev}
      tris[i].dt = this->dt;

      triIndex[3 * i + 0] = faces[i]->fbegin()->vertex()->position_in_set();
      triIndex[3 * i + 1] =
          faces[i]->fbegin()->next()->vertex()->position_in_set();
      triIndex[3 * i + 2] =
          faces[i]->fbegin()->prev()->vertex()->position_in_set();
    }

    typedef aabb_tree<SPACE, swept_triangle_type> triangle_tree;
    typedef typename triangle_tree::node_type edge_node_type;
    triangle_tree tree(tris);
    // tid = omp_get_thread_num();

    for (int i = 0; i < points.size(); i++) {
      vector<int> collectedTris;
      vector<int> filteredTris;
      T minDist = 99999;
      coordinate_type ci = points[i];
      swept_point_type li(ci, vels[i]);
      li.dt = dt;
      coordinate_type veli = velocities[i];
      int indexMin = 0;
      coordinate_type cmin;
      getAllNearest<SPACE, swept_triangle_type, swept_point_type>(
          li, tree, tris, collectedTris, 0.025 * this->minLength);

      for (int j = 0; j < collectedTris.size(); j++) {
        int pj = tree.permutation[j];
        if (triIndex[3 * collectedTris[j] + 0] == i)
          continue;
        if (triIndex[3 * collectedTris[j] + 1] == i)
          continue;
        if (triIndex[3 * collectedTris[j] + 2] == i)
          continue;
        for (int k = 0; k < 3; k++) {
          int index = triIndex[3 * collectedTris[j] + k];
          coordinate_type cj = vertices[index]->coordinate();
          T dist = norm2(ci - cj);
          if (dist < minDist) {
            minDist = dist;
            indexMin = index;
            cmin = cj;
          }
        }
      }

      if (!vertices[i])
        continue;
      if (!vertices[indexMin])
        continue;
      if (vertices[indexMin]->flag > 0 || vertices[i]->flag > 0)
        continue;
      vertices[indexMin]->flag += 1;
      vertices[i]->flag += 1;
      cmin = vertices[indexMin]->coordinate();
      vertices[i]->update_normal();
      coordinate_type N = vertices[i]->normal();
      coordinate_type dc = ci - cmin;
      dc.normalize();
      if (minDist < 999) {
        // if(minDist < 999 && abs(dot(dc,N)) > 0.0){
        contact_manifold cm;
        cm.mMesh = mMesh;
        cm.p1 = i;
        cm.p2 = indexMin;
        cm.p1type = 0;
        cm.p2type = 0;
        collisions.push_back(cm);
      }
    }

    return collisions;
  }

  vector<contact_manifold> getCollisions() {
    // TIMER function//TIMER(__FUNCTION__);
    vector<contact_manifold> collisions;
    std::vector<contact_manifold *> local;
    local.resize(8);
    std::vector<int> numCollisions;
    numCollisions.resize(8);
#if 0
#pragma omp parallel
      {
	int np = omp_get_num_threads();
	int tid  = omp_get_thread_num();
#endif
    std::vector<contact_manifold> localCollisions;

    vector<vertex_ptr> &vertices = mMesh->get_vertices();
    vector<coordinate_type> points;
    points.resize(vertices.size());

    //#pragma omp parallal for
    for (int i = 0; i < vertices.size(); i++) {
      points[i] = vertices[i]->coordinate();
    }
    vector<face_ptr> &faces = mMesh->get_faces();
    vector<int> triIndex;
    triIndex.resize(3 * faces.size());
    vector<triangle_type> tris;
    tris.resize(faces.size());

    //#pragma omp parallal for
    for (int i = 0; i < faces.size(); i++) {
      tris[i][0] = faces[i]->fbegin()->coordinate();         // 0
      tris[i][1] = faces[i]->fbegin()->next()->coordinate(); // 1
      tris[i][2] =
          faces[i]->fbegin()->prev()->coordinate(); // 2 {next->next == prev}

      triIndex[3 * i + 0] = faces[i]->fbegin()->vertex()->position_in_set();
      triIndex[3 * i + 1] =
          faces[i]->fbegin()->next()->vertex()->position_in_set();
      triIndex[3 * i + 2] =
          faces[i]->fbegin()->prev()->vertex()->position_in_set();
    }

    typedef aabb_tree<SPACE, triangle_type> triangle_tree;
    typedef typename triangle_tree::node_type edge_node_type;
    triangle_tree tree(tris);

    //#pragma omp parallal for
    for (int i = 0; i < points.size(); i++) {
#if 0
	  int tid = omp_get_thread_num();
#endif
      vector<int> collectedTris;
      T minDist = 99999;
      coordinate_type ci = points[i];
      int indexMin = 0;
      coordinate_type cmin;
      getAllNearestTriPoint<SPACE>(ci, tree, tris, collectedTris, tol);
      for (int j = 0; j < collectedTris.size(); j++) {
        int pj = tree.permutation[j];
        if (triIndex[3 * collectedTris[j] + 0] == i)
          continue;
        if (triIndex[3 * collectedTris[j] + 1] == i)
          continue;
        if (triIndex[3 * collectedTris[j] + 2] == i)
          continue;
        for (int k = 0; k < 3; k++) {
          int index = triIndex[3 * collectedTris[j] + k];
          coordinate_type cj = vertices[index]->coordinate();
          T dist = norm2(ci - cj);
          if (dist < minDist) {
            minDist = dist;
            indexMin = index;
            cmin = cj;
          }
        }
      }

      if (!vertices[i])
        continue;
      if (!vertices[indexMin])
        continue;
      if (vertices[indexMin]->flag > 0 || vertices[i]->flag > 0)
        continue;
      vertices[indexMin]->flag += 1;
      vertices[i]->flag += 1;
      cmin = vertices[indexMin]->coordinate();
      vertices[i]->update_normal();
      coordinate_type N = vertices[i]->normal();
      coordinate_type dc = ci - cmin;
      dc.normalize();

      if (minDist < 999 && abs(va::dot(dc, N)) > 0.0) {
        contact_manifold cm;
        cm.mMesh = mMesh;
        cm.p1 = i;
        cm.p2 = indexMin;
        cm.p1type = 0;
        cm.p2type = 0;
        // localCollisions.push_back(cm);
        collisions.push_back(cm);
      }
    }
#if 0
	local[omp_get_thread_num()] = new contact_manifold[localCollisions.size()];
	numCollisions[omp_get_thread_num()] = localCollisions.size();
	for(int k = 0; k < localCollisions.size(); k++){
	  local[omp_get_thread_num()][k] = localCollisions[k];
	}
#endif
#if 0	
      }
#endif
#if 0	    
      for(int p = 0; p < 8; p++){
	for(int i = 0; i < numCollisions[p]; i++){
	  collisions.push_back(local[p][i]);
	}
	delete local[p];
      }
#endif
    return collisions;
  }

  bool checkEdgeInSet(face_vertex_ptr fv0, face_vertex_ptr fv1,
                      face_vertex_ptr fva[4]) {
    vertex_ptr v0 = fv0->vertex();
    vertex_ptr v1 = fv1->vertex();
    for (int i = 0; i < 4; i++) {
      vertex_ptr v0i = fva[i]->vertex();
      vertex_ptr v1i = fva[i]->next()->vertex();
      if (v0 == v0i && v1 == v1i)
        return true;
      if (v1 == v0i && v0 == v1i)
        return true;
    }
    return false;
  }

  void pipeEdge(edge_ptr e0, edge_ptr e1, T d00, T d01, T d10, T d11,
                dynamic_octree<SPACE, triangle_type> &faceTree,
                vector<triangle_type> &corners) {
    // if(e0->v1()->vertex() == e1->v1()->vertex()) return;
    // if(e0->v1()->vertex() == e1->v2()->vertex()) return;
    // if(e0->v2()->vertex() == e1->v1()->vertex()) return;
    // if(e0->v2()->vertex() == e1->v2()->vertex()) return;
    coordinate_type c00 = e0->v1()->vertex()->coordinate();
    coordinate_type c01 = e0->v2()->vertex()->coordinate();
    coordinate_type c10 = e1->v1()->vertex()->coordinate();
    coordinate_type c11 = e1->v2()->vertex()->coordinate();
    coordinate_type cAvg = 0.25 * (c00 + c01 + c10 + c11);
    m2::Debugger &debug = m2::Debugger::get_instance();
    // std::cout << " current join: " << debug.labelId << std::endl;
    debug.add_id(cAvg[0], cAvg[1], cAvg[2]);

    T e00 = norm(c00 - c10);
    T e11 = norm(c01 - c11);
    T e01 = norm(c00 - c11);
    T e10 = norm(c01 - c10);
    face_vertex_ptr fv00 = e0->v2()->next();
    face_vertex_ptr fv10 = e1->v2()->next();
    face_vertex_ptr fv01 = e0->v1()->next();
    face_vertex_ptr fv11 = e1->v1()->next();
    face_vertex_ptr fv00p = e0->v1()->prev();
    face_vertex_ptr fv10p = e1->v1()->prev();
    face_vertex_ptr fv01p = e0->v2()->prev();
    face_vertex_ptr fv11p = e1->v2()->prev();

    face_vertex_ptr fva[8];
    fva[0] = fv00;
    fva[1] = fv10;
    fva[2] = fv01;
    fva[3] = fv11;
    fva[4] = fv00p;
    fva[5] = fv10p;
    fva[6] = fv01p;
    fva[7] = fv11p;
    for (int i = 0; i < 8; i++) {
      if (e0->v1() == fva[i])
        return;
      if (e0->v2() == fva[i])
        return;
      if (e1->v1() == fva[i])
        return;
      if (e1->v2() == fva[i])
        return;
    }
    // std::cout << e0->v1() << " " << e0->v2() << " " << e1->v1() << " " <<
    // e1->v2() << std::endl; std::cout << fv00 << " " << fv11 << " " << fv01 <<
    // " " << fv10 << std::endl; std::cout << fv00p << " " << fv11p << " " <<
    // fv01p << " " << fv10p << std::endl;

    construct<SPACE> cons;
    face_vertex_ptr fva0[4];
    face_vertex_ptr fva1[4];

    fva0[0] = e0->v1()->next();
    fva0[1] = e0->v1()->prev();
    fva0[2] = e0->v2()->next();
    fva0[3] = e0->v2()->prev();
    fva1[0] = e1->v1()->next();
    fva1[1] = e1->v1()->prev();
    fva1[2] = e1->v2()->next();
    fva1[3] = e1->v2()->prev();

    int jb = 0;
    T minE = 999; // beginning i and j indices
    for (int i = 0; i < 4; i++) {
      fva0[i]->vertex()->flag = 1;
      fva1[i]->vertex()->flag = 1;
    }
#if 1
    for (int i = 0; i < 4; i++) {
      T E = 0;
      for (int k = 0; k < 4; k++) {
        int km = (k - 1 + 4) % 4;
        int j = (i - k + 4) % 4;
        int jm = (i - k + 1 + 4) % 4;
        coordinate_type cik = fva0[k]->coordinate();
        coordinate_type cimk = fva0[km]->coordinate();
        coordinate_type cjk = fva1[j]->coordinate();

        bool checkComp0 = checkEdgeInSet(fva0[k], fva1[j], fva0);
        bool checkComp1 = checkEdgeInSet(fva0[k], fva1[j], fva1);
        if (!checkEdgeInSet(fva0[k], fva1[j], fva0) &&
            !checkEdgeInSet(fva0[k], fva1[j], fva1)) {
          coordinate_type dc = cik - cjk;
          E += va::dot(dc, dc);
        }
#if 1
        if (!checkEdgeInSet(fva0[km], fva1[j], fva0) &&
            !checkEdgeInSet(fva0[km], fva1[j], fva1)) {
          coordinate_type dcm = cimk - cjk;
          E += va::dot(dcm, dcm);
        }
#endif
      }
      if (E < minE) {
        minE = E;
        jb = i;
      }
    }

#endif

    vector<int> ppairs;
    for (int i = 0; i < 4; i++) {
      int ii = (i) % 4;
      int jj = (jb - i + 4) % 4;
      int im = (i - 1 + 4) % 4;
      int jm = (jb - i + 1 + 4) % 4;
      if (!checkEdgeInSet(fva0[ii], fva1[jj], fva0) &&
          !checkEdgeInSet(fva0[ii], fva1[jj], fva1)) {
        ppairs.push_back(ii);
        ppairs.push_back(jj);
      }
#if 1

      if (!checkEdgeInSet(fva0[im], fva1[jj], fva0) &&
          !checkEdgeInSet(fva0[im], fva1[jj], fva1)) {
        ppairs.push_back(im);
        ppairs.push_back(jj);
      }
#endif
    }
    std::cout << std::endl;
    bool hit = false;
    for (int i = 0; i < ppairs.size() / 2; i++) {
      face_vertex_ptr fv0 = fva0[ppairs[2 * i + 0]];
      face_vertex_ptr fv1 = fva1[ppairs[2 * i + 1]];

      coordinate_type ep0 = fv0->coordinate();
      coordinate_type ep1 = fv1->coordinate();
      if (m2::intersectLineTest<SPACE, triangle_type>(ep0, ep1, faceTree,
                                                      corners))
        hit = true;
      if (hit)
        return;
    }
    // the runs that made it through the gauntlet
    debug.DebugLines0.push_back(e0->v1()->coordinate());
    debug.DebugLines0.push_back(e0->v2()->coordinate());
    debug.DebugLines0.push_back(e1->v1()->coordinate());
    debug.DebugLines0.push_back(e1->v2()->coordinate());

    debug.DebugTriangles0.push_back(fva0[0]->coordinate());
    debug.DebugTriangles0.push_back(fva0[1]->coordinate());
    debug.DebugTriangles0.push_back(fva0[2]->coordinate());
    debug.DebugTriangles0.push_back(fva0[2]->coordinate());
    debug.DebugTriangles0.push_back(fva0[3]->coordinate());
    debug.DebugTriangles0.push_back(fva0[0]->coordinate());

    debug.DebugTriangles0.push_back(fva1[0]->coordinate());
    debug.DebugTriangles0.push_back(fva1[1]->coordinate());
    debug.DebugTriangles0.push_back(fva1[2]->coordinate());
    debug.DebugTriangles0.push_back(fva1[2]->coordinate());
    debug.DebugTriangles0.push_back(fva1[3]->coordinate());
    debug.DebugTriangles0.push_back(fva1[0]->coordinate());

    face_ptr f0 = cons.delete_edge(mMesh, e0);
    face_ptr f1 = cons.delete_edge(mMesh, e1);
    vector<edge_ptr> collectedEdges;
    for (int i = 0; i < ppairs.size() / 2; i++) {
      face_vertex_ptr fv0 = fva0[ppairs[2 * i + 0]];
      face_vertex_ptr fv1 = fva1[ppairs[2 * i + 1]];
      m2::Debugger &debug = m2::Debugger::get_instance();
      debug.DebugLines.push_back(fv0->coordinate());
      debug.DebugLines.push_back(fv1->coordinate());
      edge_ptr ei = cons.insert_edge(mMesh, fv0, fv1);
      collectedEdges.push_back(ei);
    }

    for (int i = 0; i < collectedEdges.size(); i++) {
      edge_ptr ei = collectedEdges[i];
      face_vertex_ptr fv0 = ei->v1();
      face_vertex_ptr fv1 = ei->v2();
      // m2::Debugger& debug = m2::Debugger::get_instance();
      // debug.DebugTriangles.push_back(fv0->coordinate());
      // debug.DebugTriangles.push_back(fv0->next()->coordinate());
      // debug.DebugTriangles.push_back(fv0->next()->next()->coordinate());
      // debug.DebugTriangles.push_back(fv1->coordinate());
      // debug.DebugTriangles.push_back(fv1->next()->coordinate());
      // debug.DebugTriangles.push_back(fv1->next()->next()->coordinate());
      corners.push_back(makeTriangle(fv0, fv0->next(), fv0->next()->next()));
      corners.push_back(makeTriangle(fv1, fv1->next(), fv1->next()->next()));
      faceTree.insert(corners, corners.size() - 1, 8);
      faceTree.insert(corners, corners.size() - 2, 8);
      // faceTree.insert(makeTriangle(fv1, fv1->next(),fv1->next()->next()));
      // if(fv0->face()->size() > 3){
      //   face_vertex_ptr fvn = fv0->next()->next();
      //   debug.DebugTriangles.push_back(fvn->coordinate());
      //   debug.DebugTriangles.push_back(fvn->next()->coordinate());
      //   debug.DebugTriangles.push_back(fvn->next()->next()->coordinate());
      //   insertedTriangles.push_back(triangle(fvn,
      //   fvn->next(),fvn->next()->next()));
      // }
    }
  }

  bool join() {
    bool topology_change = false;
    vector<contact_manifold> collisions;
    if (this->velocities.size() > 0)
      collisions = this->getCollisions(this->velocities);
    else
      collisions = this->getCollisions();
    vector<vertex_ptr> verts = mMesh->get_vertices();

    for (int i = 0; i < collisions.size(); i++) {
      coordinate_type c1 = verts[collisions[i].p1]->coordinate();
      coordinate_type c2 = verts[collisions[i].p2]->coordinate();
      Debugger &debug = Debugger::get_instance();
      debug.DebugPoints.push_back(c1);
      debug.DebugPoints.push_back(c2);
      debug.DebugLines.push_back(c1);
      debug.DebugLines.push_back(c2);
    }
    if (collisions.size() > 0)
      topology_change = true;

    std::cout << " - joining: " << collisions.size() << " pairs" << std::endl;
    mContactSorter.mMesh = mMesh;
    std::sort(collisions.begin(), collisions.end(), mContactSorter);
#if 1
    if (collisions.size() == 0)
      return false;
    int sharedEdges = 0;
    int sepEdges = 0;
    for (int i = 0; i < collisions.size(); i++) {
      vertex_ptr v1 = verts[collisions[i].p1];
      vertex_ptr v2 = verts[collisions[i].p2];
      construct<SPACE> cons;
      // std::cout << v1->position_in_set() << " " << v1->size() << " | "
      //	  << v2->position_in_set() << " " << v2->size() << " " <<
      // std::endl;
      if (v1 == v2)
        continue;
      if (v1->size() < 3)
        continue;
      if (v2->size() < 3)
        continue;

      if (v1->shares_edge_with(v2)) {
        edge_ptr e = v1->get_shared_edge(v2);
        cons.collapse_edge(mMesh, e);
        sharedEdges++;
      } else {
        // we have to deal with two odd topologies I can't account for, yet
        // f1 and f2 somehow are equal, yet they didn't start from the same
        // valence set v2 is somehow in f1 after the vertex deletion, so we
        // abort and let triangulation reclaim the the new face
        face_ptr f1 = cons.delete_vertex_primitive(mMesh, v1);
        if (!f1)
          continue;
        if (f1->has_vertex(v2))
          continue;
        face_ptr f2 = cons.delete_vertex_primitive(mMesh, v2);
        if (!f2)
          continue;
        if (f1 != f2)
          stitchFaces(mMesh, f1, f2);
        sepEdges++;
      }
    }

    std::cout << "  -collapse: " << sharedEdges
              << " edges, and piped: " << sepEdges << " edges" << std::endl;
#endif

    if (topology_change) {
      m2::remesh<space3> rem;
      rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    return topology_change;
  }

  T tol;
  control_ptr mMesh;
  dynamic_octree<SPACE, triangle_type> exteriorTree;
  dynamic_octree<SPACE, triangle_type> interiorTree;
  vector<triangle_type> oldFaces;
  vector<triangle_type> newFaces;
  vector<coordinate_type> velocities;
  T dt;
  T minLength;
}; // join_mesh;

template <typename SPACE> class calc_sdf {
  M2_TYPEDEFS
public:
  typedef aabb_tree<SPACE, triangle_type> triangle_tree;
  typedef typename triangle_tree::node_type edge_node_type;
  control_ptr mMesh;

  vector<coordinate_type> normals;
  vector<coordinate_type> vertices;
  vector<triangle_type> tris;
  vector<int> triIndices;
  triangle_tree tree = triangle_tree(this->tris);

  calc_sdf(control_ptr targetMesh) { this->mMesh = targetMesh; };

  void initialize() {
    std::cout << " in init sequence" << std::endl;
    this->normals = mMesh->get_normals();
    this->vertices = mMesh->get_coordinates();

    this->tris = this->getTris(mMesh);
    this->triIndices = this->getTriIndices(mMesh);
    this->tree = triangle_tree(this->tris);
  }

  vector<triangle_type> getTris(control_ptr targetMesh) {
    vector<face_ptr> &faces = targetMesh->get_faces();
    vector<triangle_type> ltris;
    ltris.resize(faces.size());
    for (int i = 0; i < faces.size(); i++) {
      ltris[i][0] = faces[i]->fbegin()->coordinate();         // 0
      ltris[i][1] = faces[i]->fbegin()->next()->coordinate(); // 1
      ltris[i][2] =
          faces[i]->fbegin()->prev()->coordinate(); // 2 {next->next == prev}
    }
    return ltris;
  }

  vector<int> getTriIndices(control_ptr targetMesh) {
    vector<face_ptr> &faces = targetMesh->get_faces();
    vector<int> ltriIndex;
    ltriIndex.resize(3 * faces.size());
    for (int i = 0; i < faces.size(); i++) {
      ltriIndex[3 * i + 0] = faces[i]->fbegin()->vertex()->position_in_set();
      ltriIndex[3 * i + 1] =
          faces[i]->fbegin()->next()->vertex()->position_in_set();
      ltriIndex[3 * i + 2] =
          faces[i]->fbegin()->prev()->vertex()->position_in_set();
    }
    return ltriIndex;
  }

  inline T clamp(T x, T a, T b) { return x < a ? a : (x > b ? b : x); };

  inline void distanceFromTriangle(coordinate_type p, int tri, T &d,
                                   coordinate_type &Ni) {
    int i0 = this->triIndices[3 * tri + 0];
    int i1 = this->triIndices[3 * tri + 1];
    int i2 = this->triIndices[3 * tri + 2];
    const coordinate_type &v0 = this->vertices[i0];
    const coordinate_type &v1 = this->vertices[i1];
    const coordinate_type &v2 = this->vertices[i2];

    coordinate_type u = v1 - v0;
    coordinate_type v = v2 - v0;
    coordinate_type N = va::cross(u, v);
    T iN2 = 1.0 / (va::dot(N, N));
    coordinate_type w = p - v0;
    T b10 = va::dot(va::cross(u, w), N) * iN2;
    T b20 = va::dot(va::cross(w, v), N) * iN2;
    T b12 = 1.0 - b10 - b20;

    if (b10 <= 0) {
      // line - v0, v1
      b20 = va::dot((p - v0), u) / (va::dot(u, u));
      b20 = clamp(b20, 0.0, 1.0);
      b12 = 1.0 - b20;
      b10 = 0;
    } else if (b20 <= 0) {
      // line - v0, v2
      b10 = va::dot((p - v0), v) / (va::dot(v, v));
      b10 = clamp(b10, 0.0, 1.0);
      b12 = 1.0 - b10;
      b20 = 0;
    } else if (b12 <= 0) {
      // line - v1, v2
      coordinate_type x = v2 - v1;
      b10 = va::dot((p - v1), x) / (va::dot(x, x));
      b10 = clamp(b10, 0.0, 1.0);
      b20 = 1.0 - b10;
      b12 = 0;
    }

    coordinate_type c = b10 * v2 + b20 * v1 + b12 * v0;
    // std::cout << tri << " | " << c << " | " << b10 << " " << b20 << " " <<
    // b12 << " " << b10 + b20 + b12 << std::endl;

    m2::Debugger &debug = m2::Debugger::get_instance();
    // debug.DebugLines.push_back(p);
    // debug.DebugLines.push_back(c);
    Ni = b12 * this->normals[i0] + b20 * this->normals[i1] +
         b10 * this->normals[i2];
    Ni.normalize();
    d = va::dot(p - c, Ni);
  }

  vector<pair<T, coordinate_type>>
  getSdf(vector<coordinate_type> &coordinates) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &vertices = mMesh->get_vertices();
    vector<coordinate_type> points;
    points.resize(coordinates.size());
    vector<pair<T, coordinate_type>> out;
    out.resize(coordinates.size());

    //#pragma omp parallal for
    for (int i = 0; i < coordinates.size(); i++) {
      int cId = getNearest<SPACE, triangle_type>(coordinates[i], this->tree,
                                                 this->tris);
      coordinate_type N;
      T d;
      this->distanceFromTriangle(coordinates[i], cId, d, N);
      pair<T, coordinate_type> dN;

      dN.first = d;
      dN.second = N;
      out[i] = dN;
    }
    return out;
  }

}; // calc sdf;

template <typename SPACE> class moving_mesh {
  M2_TYPEDEFS;

  typedef coordinate_type point_type;
  typedef typename SPACE::triangle_type triangle;
  typedef typename list<coordinate_type>::iterator pl_iterator;
  typedef typename m2::Debugger Debugger;

public:
  moving_mesh(control_ptr obj_in) {
    mMesh = obj_in;
    maxCurvature = 3.0;
    minCurvature = 0.001;
    minLength = 0.01;
    minCollapseLength = 0.00025;
    maxLength = 0.06;

    coordinate_type cen = mMesh->calc_center();
    coordinate_type min = mMesh->calc_min();
    coordinate_type max = mMesh->calc_max();
    coordinate_type lengths = max - min;
    int minRes = 64;
    T minLength = lengths[0];
    minLength = minLength < lengths[1] ? minLength : lengths[1];
    minLength = minLength < lengths[2] ? minLength : lengths[2];

    T dx = minLength / (T)(minRes - 1);
    remesh<SPACE> rem;
    rem.triangulate(mMesh);
    // construct<SPACE> cons;
    // face_ptr nf = cons.delete_vertex(mMesh,mMesh->get_vertices()[0]);
    // vector<T> edgeWeights;
    // vector<T> vertexWeights;
    // calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);
  }

  face_bin<SPACE> &getHash() { return mFaceHash; }

  void pin_half() {
    mMesh->pack();
    vertex_array &vl = mMesh->get_vertices();
    coordinate_type cen(0, 0, 0);
    for (int j = 0; j < vl.size(); j++) {
      cen += vl[j]->coordinate();
      vl[j]->pinned = false;
    }

    cen /= (T)vl.size();
    std::cout << " cen: " << cen << std::endl;
    for (int j = 0; j < vl.size(); j++) {
      if (vl[j]->coordinate()[1] < cen[1]) {
        vl[j]->pinned = true;
      }
    }
  }

  void drawField() {
    vector<vertex_ptr> &tverts = mMesh->get_vertices();
    if (tverts.size() > 0 && mDirectionWeights.size() > 0 &&
        mDirectionField.size() > 0) {
      for (long i = 0; i < tverts.size(); i++) {

        coordinate_type p0 = tverts[i]->coordinate();
        coordinate_type c = 0.0035 * mDirectionWeights[i];
        coordinate_type dp0 = p0 - c[0] * mDirectionField[i][0];
        coordinate_type dp1 = p0 + c[0] * mDirectionField[i][0];
        glBegin(GL_LINES);

        glColor4f(0.25f, 0.25f, 0.8f, 0.5f);
        glVertex3f(dp0[0], dp0[1], dp0[2]);
        glVertex3f(dp1[0], dp1[1], dp1[2]);
        glEnd();

        dp0 = p0 - c[1] * mDirectionField[i][1];
        dp1 = p0 + c[1] * mDirectionField[i][1];
        glBegin(GL_LINES);
        glColor4f(0.25f, 0.8f, 0.25, 0.5f);
        glVertex3f(dp0[0], dp0[1], dp0[2]);
        glVertex3f(dp1[0], dp1[1], dp1[2]);
        glEnd();

        dp0 = p0 - c[2] * mDirectionField[i][2];
        dp1 = p0 + c[2] * mDirectionField[i][2];
        glBegin(GL_LINES);
        glColor4f(0.8f, 0.25f, 0.25, 0.5f);
        glVertex3f(dp0[0], dp0[1], dp0[2]);
        glVertex3f(dp1[0], dp1[1], dp1[2]);
        glEnd();
      }
    }
  }

  coordinate_type getSvdNormal(vertex_ptr v) {
    if (v->size() == 0)
      return coordinate_type(0, 0, 0);
    face_vertex_ptr itb = v->fbegin();
    face_vertex_ptr ite = v->fend();
    bool at_head = false;
    coordinate_type kA(0, 0, 0);
    T wTotal = 0;
    coordinate_type c0 = v->coordinate();
    int k = 0;
    coordinate_type cov[3];
    while (!at_head) {
      at_head = itb == ite;
      T wij = itb->face()->area();
      coordinate_type c1 = itb->next()->coordinate();
      wTotal += wij;
      kA += wij * (c1 - c0);

      coordinate_type dc = itb->face()->normal();
      for (int l = 0; l < 3; l++)
        for (int m = 0; m < 3; m++) {
          cov[l][m] += dc[l] * dc[m];
        }
      itb = itb->vnext();
      k++;
    }

    if (wTotal > 1e-10) {
      v->update_normal();
      coordinate_type w;
      calcSVD<SPACE>(cov, w);
      coordinate_type Nv = v->normal();
      coordinate_type N = cov[0];
      if (va::dot(Nv, N) < 0.0)
        return -N;
      else
        return N;
    } else
      return v->normal();
  }

  vector<coordinate_type> getFaceOffsets(m2::control<SPACE> &in) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<T> vertexWeights;
    vector<T> edgeWeights;
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> normals;
    normals.resize(tverts.size());
    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;
      normals[i] = getSvdNormal(v);
    }
    return normals;
  }

  coordinate_type getButterflyVertexWeight(face_vertex_ptr v) {
    // \/\./
    // vn = next
    // vnn = vertext_next next
    // pvpp = prev vertex_prev next
    // coordinate_type c = v->coordinate();
    // coordinate_type cn = v->next()->coordinate();
    // coordinate_type cvnn = v->vnext()->next()->coordinate();
    // coordinate_type cpvpp = v->prev()->vprev()->prev()->coordinate();
    // return (8.0*c + 2.0*cn - 1.0*(cvnn + cpvpp));
    face_vertex_ptr itb = v;
    face_vertex_ptr ite = v->vprev();
    bool iterating = true;
    int K = 0, j = 0;
    while (iterating) {
      iterating = itb != ite;
      itb = itb->vnext();
      K++;
    }
    itb = v;
    iterating = true;
    T s = 0;
    coordinate_type cs(0, 0, 0);
    while (iterating) {
      iterating = itb != ite;
      itb = itb->vnext();
      T frac = (T)j / (T)K;
      T sj =
          (0.25 + cos(2.0 * M_PI * frac) + 0.5 * cos(4.0 * M_PI * frac)) / (T)K;
      coordinate_type ci = itb->next()->coordinate();
      cs += sj * ci;
      s += sj;
      j++;
    }
    coordinate_type c = v->coordinate();
    return 0.75 * c + cs;
  }

  coordinate_type getButterflyWeight(edge_ptr e) {
    face_vertex_ptr v1 = e->v1();
    face_vertex_ptr v2 = e->v2();
    coordinate_type c1 = getButterflyVertexWeight(v1);
    coordinate_type c2 = getButterflyVertexWeight(v2);
    return 0.5 * (c1 + c2);
  }

  coordinate_type getAverageWeight(edge_ptr e) {
    face_vertex_ptr v1 = e->v1();
    face_vertex_ptr v2 = e->v2();
    coordinate_type c1 = v1->coordinate();
    coordinate_type c2 = v2->coordinate();
    return 0.5 * (c1 + c2);
  }

  void checkEdge(int edgeNum) {
    edge_array &edges = mMesh->get_edges();
    if (edgeNum < edges.size()) {
      edge_ptr ei = mMesh->get_edges()[edgeNum];
      std::cout << "edge check: " << ei->v1()->vertex() << " "
                << ei->v2()->vertex() << " " << ei->v1()->vertex()->size()
                << " " << ei->v2()->vertex()->size() << " " << ei->length()
                << std::endl;
    }
  }

  void relax_mesh() {
    bool relaxing = true;
    int k = 0;
    // TIMER function//TIMER(__FUNCTION__);
    std::cout << " mesh size: " << mMesh->get_vertices().size() << std::endl;
    // relaxing = this->delete_vertices_curve();

    // relaxing = this->insert_edges();

    mMesh->colorVertices();
    std::cout << " max graph color: " << mMesh->maxGraphColor << std::endl;
    while (relaxing) {
      relaxing = this->collapse_edges();
      delete_degenerate(); //
    }

    relaxing = true;
    while (relaxing)
      relaxing = this->insert_edges();

    k = 0;
    relaxing = true;
    while (relaxing && k < 2) {
      set_operations<SPACE> setOps;
      relaxing = setOps.flip_edges(mMesh);
      k++;
    }
    for (int i = 0; i < 5; i++)
      surface_filter<SPACE>::filter(*mMesh, 0.01);

    this->joinEdges();
    // surface_calculator<SPACE>::shadeVertices(*mMesh);
    // 23765
  }

  void cache_mesh() {
    face_array &faces = mMesh->get_faces();
    Debugger &debug = Debugger::get_instance();
    for (int i = 0; i < faces.size(); i++) {
      debug.CachedTriangles.push_back(
          faces[i]->fbegin()->vertex()->coordinate());
      debug.CachedTriangles.push_back(
          faces[i]->fbegin()->next()->vertex()->coordinate());
      debug.CachedTriangles.push_back(
          faces[i]->fbegin()->next()->next()->vertex()->coordinate());
    }
  }

  bool join(vector<coordinate_type> velocities, T dt) {
    join_mesh<SPACE> joiner(mMesh, 0.01 * minLength);
    joiner.init();
    joiner.testCollision();
    joiner.velocities = velocities;
    joiner.minLength = this->minLength;
    joiner.dt = dt;
    bool hit = joiner.join();
    int i = 0;
    bool relaxing = true;
    while (relaxing && i < 20) {
      relaxing = delete_degenerate();
      i++;
    }
    return hit;
  }

  bool join() {
    join_mesh<SPACE> joiner(mMesh, 0.02 * minLength);
    joiner.init();
    bool hit = joiner.join();
    int i = 0;
    bool relaxing = true;
    while (relaxing && i < 20) {
      relaxing = delete_degenerate();
      i++;
    }
    return hit;
  }

  edge_ptr subdivideFace(face_vertex_ptr fv) {
    // TIMER function//TIMER(__FUNCTION__);
    // assuming that the face vertex is the newly inserted one.
    face_vertex_ptr fv1 = fv; // this is the
    face_vertex_ptr fv2 = fv->next();
    face_vertex_ptr fv3 = fv2->next();
    face_vertex_ptr fv4 = fv3->next();
    m2::construct<SPACE> cons;
    edge_ptr enew = cons.insert_edge(mMesh, fv1, fv3);
    return enew;
  }

  void split_edge(edge_ptr e) {
    // TIMER function//TIMER(__FUNCTION__);
    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();
    coordinate_type dp =
        0.5 * (e->v1()->vertex()->data + e->v2()->vertex()->data);
    T dp2 = 0.5 * (e->v1()->vertex()->data2 + e->v2()->vertex()->data2);
    coordinate_type c = getButterflyWeight(e);
    construct<SPACE> cons;
    vertex_ptr nv = cons.subdivide_edge(mMesh, e);
    nv->coordinate() = c;
    nv->data = dp;
    nv->data2 = dp2;
    subdivideFace(fv1->next());
    subdivideFace(fv2->next());
  }

  bool insert_edges() {
    // TIMER function//TIMER(__FUNCTION__);
    edge_array &edges = mMesh->get_edges();
    vector<vertex_ptr> &verts = mMesh->get_vertices();

    face_array &faces = mMesh->get_faces();
    int N = faces.size();

    bool topology_change = false;
    vector<T> edgeWeights;
    vector<T> vertexWeights;
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();
    surface_calculator<SPACE> calc;

    vector<edge_ptr> edgesToSplit;
#if 1
    for (int i = 0; i < edges.size(); i++) {
      if (!mMesh->has_edge(i))
        continue;
      edge_ptr ei = edges[i];
      bool pinned = ei->v1()->vertex()->pinned == true &&
                    ei->v2()->vertex()->pinned == true;
      if (pinned)
        continue;
      T l = ei->length();
      if (l > 1.75 * minLength) {
        edgesToSplit.push_back(ei);
        continue;
      }
    }
#endif

    std::cout << "sorting " << edgesToSplit.size() << " edges" << std::endl;
    std::sort(edgesToSplit.begin(), edgesToSplit.end(), mEdgeSorter);
    std::cout << "splitting " << edgesToSplit.size() << " edges" << std::endl;
    for (int i = edgesToSplit.size(); i > 0; i--) {
      this->split_edge(edgesToSplit[i - 1]);
    }

    return topology_change;
  }

  bool collapse_edges() {
    // TIMER function//TIMER(__FUNCTION__);
    bool topology_change = false;
    edge_array collectedEdges;
    mMesh->reset_flags();
    vector<edge_ptr> &edges = mMesh->get_edges();
    for (int i = 0; i < edges.size(); i++) {

      if (mMesh->has_edge(i) > 0) {
        edge_ptr e = edges[i];
        if (e->v1()->vertex()->flag == 1 || e->v2()->vertex()->flag == 1
            // e->v1()->vertex()->size() < 4 ||
            // e->v2()->vertex()->size() < 4
        ) {
          continue;
        }
        T dist = e->dist();
        T distNext = e->v1()->next()->edge()->dist();
        vertex_ptr v = e->v1()->vertex();
        if (dist < minCollapseLength) {
          e->v1()->vertex()->flag = 1;
          e->v2()->vertex()->flag = 1;
          collectedEdges.push_back(e);
        }
      }
    }

    std::cout << "deleting: " << collectedEdges.size() << " Tiny edges"
              << std::endl;

    for (int i = 0; i < collectedEdges.size(); i++) {
      if (mMesh->has_edge(i)) {
        construct<SPACE> cons;
        edge_ptr e = collectedEdges[i];
        coordinate_type dp =
            0.5 * (e->v1()->vertex()->data + e->v2()->vertex()->data);
        T dp2 = 0.5 * (e->v1()->vertex()->data2 + e->v2()->vertex()->data2);

        vertex_ptr vi = cons.collapse_edge(mMesh, e);
        vi->data = dp;
        vi->data2 = dp2;
        topology_change = true;
      }
    }
    if (topology_change) {
      m2::remesh<SPACE> rem;
      rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    return topology_change;
  }

  bool collapse_low_curve() {
    // TIMER function//TIMER(__FUNCTION__);
    vector<T> edgeWeights;
    vector<T> vertexWeights;
    m2::surface_calculator<SPACE> calc;
    calc.calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);
    bool topology_change = false;
    edge_array collectedEdges;
    mMesh->reset_flags();
    vector<edge_ptr> &edges = mMesh->get_edges();
    for (int i = 0; i < edges.size(); i++) {

      if (mMesh->has_edge(i) > 0) {
        edge_ptr e = edges[i];
        if (e->v1()->vertex()->flag == 1 || e->v2()->vertex()->flag == 1) {
          continue;
        }
        T dist = e->dist();
        T distNext = e->v1()->next()->edge()->dist();
        vertex_ptr v = e->v1()->vertex();
        int iv0 = e->v1()->vertex()->position_in_set();
        int iv1 = e->v2()->vertex()->position_in_set();
        T kAvg = 0.5 * (vertexWeights[iv0] + vertexWeights[iv1]);
        if (kAvg < 200.0) {
          e->v1()->vertex()->flag = 1;
          e->v2()->vertex()->flag = 1;
          collectedEdges.push_back(e);
        }
      }
    }

    std::cout << "collapsing: " << collectedEdges.size()
              << " low curvature edges" << std::endl;

    for (int i = 0; i < collectedEdges.size(); i++) {
      if (mMesh->has_edge(i)) {
        construct<SPACE> cons;
        edge_ptr e = collectedEdges[i];
        coordinate_type dp =
            0.5 * (e->v1()->vertex()->data + e->v2()->vertex()->data);
        T dp2 = 0.5 * (e->v1()->vertex()->data2 + e->v2()->vertex()->data2);

        vertex_ptr vi = cons.collapse_edge(mMesh, e);
        vi->data = dp;
        vi->data2 = dp2;
        topology_change = true;
      }
    }
    if (topology_change) {
      m2::remesh<SPACE> rem;
      rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    return topology_change;
  }

  bool prune_edges() {
    // TIMER function//TIMER(__FUNCTION__);
    bool topology_change = false;
    // vertex_array  collectedEdges;
    // vector<T> edgeWeights;
    // vector<T> vertexWeights;
    surface_calculator<SPACE> calc;
    // calc.calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);
    mMesh->reset_flags();
    vector<edge_ptr> &edges = mMesh->get_edges();
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<edge_ptr> edgesToCollapse;
    vector<T> vertThinness;
    vertThinness.resize(verts.size(), 0.0);
    vector<coordinate_type> newCoord;
    newCoord.resize(verts.size());
    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      vertex_ptr v = verts[i];
      T vthin = v->thinness();
      vertThinness[i] = vthin;
    }

    calc.calcDiffuseQuantity(*mMesh, vertThinness, 0.2);
    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      // std::cout << vertThinness[i] << std::endl;;
    }
#if 1
    for (int i = 0; i < edges.size(); i++) {
      if (mMesh->has_edge(i) > 0) {
        edge_ptr e = edges[i];
        T dist = e->dist();
        T distNext = e->v1()->next()->edge()->dist();
        vertex_ptr v1 = e->v1()->vertex();
        vertex_ptr v2 = e->v2()->vertex();
        if (v1->pinned)
          continue;
        if (v2->pinned)
          continue;
        coordinate_type n1 = e->v1()->face()->normal();
        coordinate_type n2 = e->v2()->face()->normal();
        // T k1 = vertexWeights[e->v1()->vertex()->position_in_set()];
        // T k2 = vertexWeights[e->v2()->vertex()->position_in_set()];
        // T k1p = vertexWeights[e->v1()->prev()->vertex()->position_in_set()];
        // T k2p = vertexWeights[e->v2()->prev()->vertex()->position_in_set()];
        // T rat = (k1p+k2p)/(k1+k2);
        // if (rat < 2.0e-1){
        if (0.5 * (e->v1()->face()->thinness() + e->v2()->face()->thinness()) <
            0.98) {
          e->v1()->vertex()->flag = 1;
          e->v2()->vertex()->flag = 1;
        }
        // if (dot(n1,n2) < -0.85){
        //   // edgesToCollapse.push_back(e->v2()->prev()->edge());

        // }
      }
    }
#endif

    for (int i = 0; i < verts.size(); i++) {
      if (!mMesh->has_vertex(i))
        continue;
      if (verts[i]->size() == 1)
        continue;
      if (verts[i]->flag == 1) {
        topology_change = true;
        vertex_ptr v = verts[i];
        m2::surface_filter<SPACE> filter;
        newCoord[i] = 0.5 * (1.0 - vertThinness[i]) *
                      filter.laplacianFilterVertex(verts[i]);
      }
    }
    for (int i = 0; i < verts.size(); i++) {
      if (!mMesh->has_vertex(i))
        continue;
      if (verts[i]->size() == 1)
        continue;
      if (verts[i]->flag == 1) {
        topology_change = true;
        vertex_ptr v = verts[i];
        m2::surface_filter<SPACE> filter;
        verts[i]->coordinate() += newCoord[i];
      }
    }

    return topology_change;
  }

  bool delete_degenerate() {
    // TIMER function//TIMER(__FUNCTION__);
    bool topology_change = false;
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<edge_ptr> &edges = mMesh->get_edges();
    int zeroEdges = 0;
    int twoFaces = 0;
    int flatVol = 0;
    vector<edge_ptr> edgesToDelete;
    for (int i = 0; i < edges.size(); i++) {
      if (mMesh->has_edge(i)) {
        edge_ptr e = edges[i];
        if (e->flag == 1)
          continue;
        int s1 = e->v1()->face()->size();
        int s2 = e->v2()->face()->size();
        T a1 = e->v1()->face()->calc_area();
        T a2 = e->v2()->face()->calc_area();

        if (a1 < 1e-12 || a2 < 1e-12) {
          e->flag = 1;
          edgesToDelete.push_back(e);
          continue;
        }

        if (s1 < 3 || s2 < 3) {
          twoFaces++;
          e->flag = 1;
          edgesToDelete.push_back(e);
          continue;
        }

        if (e->v1()->face() == e->v2()->face()) {
          edgesToDelete.push_back(e);
          e->flag = 1;
          continue;
        }

        if (e->v1()->vertex() == e->v2()->vertex()) {
          zeroEdges++;
          e->flag = 1;
          edgesToDelete.push_back(e);
          continue;
        };

        if (e->v1()->prev()->vertex() == e->v2()->prev()->vertex()) {
          flatVol++;
          edgesToDelete.push_back(e);
          e->flag = 1;
          continue;
        };
      }
    }
    topology_change = edgesToDelete.size() > 0;
    construct<SPACE> cons;
    for (int i = 0; i < edgesToDelete.size(); i++) {
      edge_ptr e = edgesToDelete[i];
      face_ptr nf = cons.delete_edge(mMesh, e);
    }
#if 0
      std::cout << " - deleted: " 
		<< zeroEdges << " zero edges and " 
		<< twoFaces << " two faces and "
		<< flatVol << " flat volumes."
		<< std::endl;
#endif
#if 1
    int lowVerts = 0;
    for (int i = 0; i < verts.size(); i++) {
      if (!mMesh->has_vertex(i))
        continue;
      if (verts[i]->pinned)
        continue;
      vertex_ptr v = verts[i];

      if (verts[i]->size() < 3) {
        construct<SPACE> cons;
        face_ptr nf = cons.delete_vertex_primitive(mMesh, v);
        topology_change = true;
        continue;
      }

      // if(verts[i]->size() == 2){
      //   //std::cout << "degenerate vert: " << verts[i]->size() << std::endl;
      //   int sz = verts[i]->size();
      //   construct<SPACE> cons;

      //   edge_ptr e = v->fbegin()->next()->edge();
      //   cons.delete_edge(mMesh,e);

      //   face_ptr nf = cons.delete_vertex_primitive(mMesh,v);
      //   //numVerts++;
      //   topology_change = true;
      //   continue;
      // }
    }
#endif
    if (topology_change) {
      m2::remesh<space3> rem;
      // rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    return topology_change;
  }

  bool checkCompatibility(edge_ptr ei, edge_ptr ej) {
    int vertIds[4];
    vertIds[0] = ei->v1()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[2] =
    // ei->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

    vertIds[1] = ei->v2()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[5] =
    // ei->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

    vertIds[2] = ej->v1()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[8] =
    // ej->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

    vertIds[3] = ej->v2()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[11] =
    // ej->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        if (i != j)
          if (vertIds[i] == vertIds[j])
            return false;
      }

    return true;
  };

  //---------------------------------------------------------------------------
  // sorting operators
  //---------------------------------------------------------------------------

  struct edge_sort {
    bool operator()(edge_ptr ei, edge_ptr ej) {
      return (ei->length() < ej->length());
    }
  } mEdgeSorter;

  vector<T> calcAvgSdf(vector<coordinate_type> &coordinates,
                       vector<coordinate_type> &norms) {
    vector<T> sdf;
    vector<T> sdfWeights;

    sdf.resize(coordinates.size(), 0.0);
    sdfWeights.resize(coordinates.size(), 0.0);

    vector<face_ptr> &faces = mMesh->get_faces();
    vector<coordinate_type> chargeCenters;
    vector<T> charges;
    vector<T> chargeMags;
    vector<T> chargeK;
    vector<coordinate_type> chargeNorms;
    for (int i = 0; i < faces.size(); i++) {
      T area = faces[i]->area();
      coordinate_type norm = faces[i]->normal();
      if (area < 1e-9)
        continue;
      T charge = area;

      chargeCenters.push_back(faces[i]->center());
      charges.push_back(charge);
      chargeMags.push_back(area);
      chargeNorms.push_back(norm);
    }

    typedef pole_tree<SPACE> Tree;
    typedef pole_node<SPACE> Node;
    Tree octree(chargeCenters);

    vector<T> nodeCharges;
    vector<coordinate_type> nodeChargeCenters;
    vector<coordinate_type> nodeChargeNorms;
    nodeCharges.resize(octree.nodes.size());
    nodeChargeCenters.resize(octree.nodes.size());
    nodeChargeNorms.resize(octree.nodes.size());
    std::stack<int> stack;
    stack.push(0);

    while (stack.size() > 0) {
      int pId = stack.top();
      stack.pop();
      Node &pNode = octree.nodes[pId];
      coordinate_type chargeCenter(0, 0, 0);
      T chargeNet = 0;
      T netChargeMag = 0;
      coordinate_type normalNet(0, 0, 0);
      int N = pNode.size;
      T iN = 1.0 / (T)N;
      int beg = pNode.begin;
      for (int i = beg; i < beg + N; i++) {
        int ii = octree.permutation[i];
        T chargeMag = chargeMags[ii];
        // T chargeMag = 1.0;
        coordinate_type chargeLoc = chargeCenters[ii];
        coordinate_type chargeNorm = chargeNorms[ii];
        T charge = charges[ii];
        chargeNet += charges[ii];
        netChargeMag += chargeMag;
        chargeCenter += chargeMag * chargeLoc;
      }
      chargeCenter /= netChargeMag;
      T netDistMag = 0;
      for (int i = beg; i < beg + N; i++) {
        int ii = octree.permutation[i];
        T chargeMag = chargeMags[ii];
        coordinate_type chargeLoc = chargeCenters[ii];
        coordinate_type chargeNorm = chargeNorms[ii];
        coordinate_type dp = chargeCenters[ii] - chargeCenter;
        T dist = norm(dp);
        T regLength = 1.5 * minLength;
        T w = 1.0 / powf(dist * dist + regLength * regLength, 1.0);
        T charge = charges[ii];
        netDistMag += w * chargeMag;
        normalNet += w * chargeMag * chargeNorm;
      }

      normalNet /= netDistMag;
      T netWeight = 0;

      normalNet.normalize();
      nodeCharges[pId] = chargeNet;
      nodeChargeCenters[pId] = chargeCenter;
      nodeChargeNorms[pId] = normalNet;

      for (int j = 0; j < 8; j++) {
        if (pNode.children[j] != -1)
          stack.push(pNode.children[j]);
      }
    }
    T thresh = 0.5;

    //#pragma omp parallel for
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    for (int i = 0; i < coordinates.size(); i++) {
      int count = 0;
      // if(!verts[i]) continue;
      coordinate_type pi = coordinates[i];
      std::stack<int> stack1;
      stack1.push(0);

      double sdfi = sdf[i];
      double sdfWeightsi = sdfWeights[i];
      while (stack1.size() > 0) {

        int pId = stack1.top();
        stack1.pop();
        Node &pNode = octree.nodes[pId];
        coordinate_type pj = nodeChargeCenters[pId];
        coordinate_type Nj = nodeChargeNorms[pId];
        Nj.normalize();
        Debugger &debug = Debugger::get_instance();
        // debug.DebugLines.push_back(pj);
        // debug.DebugLines.push_back(pj+0.025*Nj);

        coordinate_type dp = pi - pj;
        T dc = norm(dp);
        T sc = norm(pNode.half);

        if (sc / dc > thresh || pNode.level < 2) {
          for (int j = 0; j < 8; j++) {
            if (pNode.children[j] != -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }

        else {

          int ii = octree.permutation[pNode.begin];
          T ci = nodeCharges[pId];

          // coordinate_type pj = pNode.centerOfMass;
          if (norm(dp) < 1e-12)
            continue;

          // Debugger& debug = Debugger::get_instance();
          // debug.DebugBoxes.push_back(pNode.center);
          // debug.DebugBoxes.push_back(pNode.half);
          // debug.DebugLines.push_back(pi);
          // debug.DebugLines.push_back(pj);

          T dist = norm(dp);
          T Px = va::dot(dp, Nj);

          T i4pi = 0.25 / M_PI;
          // T regLength = 2.0*minLength;
          // T expon = -(dist*dist)/(regLength*regLength);
          // T w = exp(expon);
          T regLength = 0.5 * minLength;
          T w = 1.0 / powf(dist * dist + regLength * regLength, 1.0);
          // T w = 1.0/powf(dist*dist + regLength*regLength,2.0);
          /*
std::cout << dp << " " << Nj
<< " Px: " << Px << " w: " << w
<< " " << dist << std::endl;
std::cout << i << " " << sdfWeightsi << " " << sdfi << std::endl;
*/
          sdfWeightsi += w;
          if (pNode.size < 512)
            sdfi += w * Px;
          else
            sdfi += w * dist;
        }
      }
      sdfWeights[i] = sdfWeightsi;
      sdf[i] = sdfi;
    }

    for (int i = 0; i < sdf.size(); i++) {
      // std::cout << " sdf: " << sdf[i] << " " << sdfWeights[i];
      sdf[i] *= 1.0 / sdfWeights[i];
    }
    return sdf;
  }

  vector<coordinate_type>
  integrateChargesBarnesHut(vector<coordinate_type> &coordinates) {
    vector<coordinate_type> u;
    u.resize(coordinates.size(), coordinate_type(0, 0, 0));

    vector<face_ptr> &faces = mMesh->get_faces();
    vector<coordinate_type> chargeCenters;
    vector<T> charges;
    vector<T> chargeMags;
    for (int i = 0; i < faces.size(); i++) {
      T area = faces[i]->area();
      coordinate_type norm = faces[i]->normal();
      if (area < 1e-9)
        continue;
      T charge = area;
      chargeCenters.push_back(faces[i]->center());
      charges.push_back(charge);
      chargeMags.push_back(area);
    }

    typedef pole_tree<SPACE> Tree;
    typedef pole_node<SPACE> Node;
    Tree octree(chargeCenters);

    vector<T> nodeCharges;
    vector<coordinate_type> nodeChargeCenters;
    nodeCharges.resize(octree.nodes.size());
    nodeChargeCenters.resize(octree.nodes.size());
    std::stack<int> stack;
    stack.push(0);

    while (stack.size() > 0) {
      int pId = stack.top();
      stack.pop();
      Node &pNode = octree.nodes[pId];
      coordinate_type chargeCenter(0, 0, 0);
      T chargeNet = 0;
      T netChargeMag = 0;
      int N = pNode.size;
      T iN = 1.0 / (T)N;
      int beg = pNode.begin;
      for (int i = beg; i < beg + N; i++) {
        int ii = octree.permutation[i];
        T chargeMag = chargeMags[ii];
        // T chargeMag = 1.0;
        coordinate_type chargeLoc = chargeCenters[ii];
        T charge = charges[ii];
        chargeNet += charges[ii];
        netChargeMag += chargeMag;
        chargeCenter += chargeMag * chargeLoc;
      }
      chargeCenter /= netChargeMag;
      T netWeight = 0;

      nodeCharges[pId] = chargeNet;
      nodeChargeCenters[pId] = chargeCenter;
      for (int j = 0; j < 8; j++) {
        if (pNode.children[j] != -1)
          stack.push(pNode.children[j]);
      }
    }
    T thresh = 0.5;
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    for (int i = 0; i < coordinates.size(); i++) {
      int count = 0;
      // if(!verts[i]) continue;
      coordinate_type pi = coordinates[i];

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int pId = stack1.top();
        stack1.pop();
        Node &pNode = octree.nodes[pId];
        coordinate_type pj = nodeChargeCenters[pId];
        coordinate_type dp = pi - pj;
        T dc = norm(dp);
        T sc = norm(pNode.half);

        if (sc / dc > thresh) {
          for (int j = 0; j < 8; j++) {
            if (pNode.children[j] != -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }

        else {

          int ii = octree.permutation[pNode.begin];
          T ci = nodeCharges[pId];

          // coordinate_type pj = pNode.centerOfMass;
          if (norm(dp) < 1e-12)
            continue;
          T dist = norm(dp);
          dp.normalize();
          T i4pi = 0.25 / M_PI;
          T regLength = 4.0 * minLength;
          T w = 1.0 / powf(dist * dist + regLength * regLength, 1.0);
          coordinate_type ui = i4pi * ci * w * dp;
          // std::cout << ui << std::endl;
          u[i] += ui;
        }
      }
    }
    return u;
  }

  vector<coordinate_type>
  integrateOrientedChargesBarnesHut(vector<coordinate_type> &coordinates) {
    vector<coordinate_type> u = integrateChargesBarnesHut(coordinates);
    vector<coordinate_type> N = integrateNormalsBarnesHut(coordinates);
    for (int i = 0; i < u.size(); i++) {
      T dotU = va::dot(u[i], N[i]);
      if (dotU < 0)
        u[i] *= -1.0;
      // u[i] *= 1.0/uSum[i];
      // u[i].normalize();
    }
    return u;
  }

  vector<coordinate_type>
  integrateNormalsBarnesHut(const vector<coordinate_type> &coordinates) {
    mMesh->update_all();
    vector<face_ptr> &faces = mMesh->get_faces();

    vector<coordinate_type> chargeCenters;
    vector<coordinate_type> charges;
    vector<T> chargeMags;
    for (int i = 0; i < faces.size(); i++) {
      T area = faces[i]->area();

      coordinate_type norm = faces[i]->normal();
      norm.normalize();
      if (area < 1e-9)
        continue;
      coordinate_type charge = area * norm;
      chargeCenters.push_back(faces[i]->center());
      charges.push_back(charge);
      chargeMags.push_back(area);
    }

    typedef pole_tree<SPACE> Tree;
    typedef pole_node<SPACE> Node;
    Tree octree(chargeCenters);

    vector<coordinate_type> nodeCharges;
    vector<coordinate_type> nodeChargeCenters;
    nodeCharges.resize(octree.nodes.size());
    nodeChargeCenters.resize(octree.nodes.size());
    std::stack<int> stack;
    stack.push(0);
    T netWeight = 0;
    while (stack.size() > 0) {
      int pId = stack.top();
      stack.pop();
      Node &pNode = octree.nodes[pId];
      coordinate_type chargeCenter(0, 0, 0);
      coordinate_type chargeNet = 0;
      T netChargeMag = 0;
      int N = pNode.size;
      T iN = 1.0 / (T)N;
      int beg = pNode.begin;
      for (int i = beg; i < beg + N; i++) {
        int ii = octree.permutation[i];
        T chargeMag = chargeMags[ii];
        // T chargeMag = 1.0;
        coordinate_type chargeLoc = chargeCenters[ii];
        coordinate_type charge = charges[ii];
        netChargeMag += chargeMag;
        chargeCenter += chargeMag * chargeLoc;
        T w = 1.0;
        chargeNet += w * charges[ii];

        netWeight += w;
      }
      chargeCenter /= netChargeMag;

      nodeCharges[pId] = chargeNet;
      nodeChargeCenters[pId] = chargeCenter;
      for (int j = 0; j < 8; j++) {
        if (pNode.children[j] != -1)
          stack.push(pNode.children[j]);
      }
    }
    T thresh = 0.5;
    vector<coordinate_type> u;
    vector<T> uSum;
    u.resize(coordinates.size(), coordinate_type(0, 0, 0));
    uSum.resize(coordinates.size(), 0.0);
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    for (int i = 0; i < coordinates.size(); i++) {
      int count = 0;
      // if(!verts[i]) continue;
      coordinate_type pi = coordinates[i];

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int pId = stack1.top();
        stack1.pop();
        Node &pNode = octree.nodes[pId];
        coordinate_type pj = nodeChargeCenters[pId];
        coordinate_type dp = pi - pj;
        T dc = norm(dp);
        T sc = norm(pNode.half);

        if (sc / dc > thresh) {
          for (int j = 0; j < 8; j++) {
            if (pNode.children[j] != -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }

        else {

          int ii = octree.permutation[pNode.begin];
          coordinate_type ci = nodeCharges[pId];
          if (norm(dp) < 1e-16)
            continue;

          if (norm(ci) < 1e-16)
            continue;
          ci.normalize();
          // coordinate_type pj = pNode.centerOfMass;
          // coordinate_type ndp = dp; ndp.normalize();
          // if(va::dot(ci,ndp) < 0.0) continue;
          // coordinate_type ciN = ci; ciN.normalize();
          // coordinate_type dpp = va::dot(dp,ciN)/va::dot(dp,dp)*dp;
          T dist = norm(dp);

          T i4pi = 0.25 / M_PI;
          T regLength = 1.0 * minLength;
          T m = 2.0;
          T w = 1.0 / powf(dist * dist + regLength * regLength, m);
          // T w = exp(-dist*dist/regLength*regLength);
          u[i] += w * ci;
          // u[i] += w*coordinate_type(0.01,0.01,0.01);
          uSum[i] += w;
        }
      }
    }

    for (int i = 0; i < u.size(); i++) {
      if (uSum[i] < 1e-9)
        continue;
      u[i] *= 1.0 / uSum[i];
      u[i].normalize();
    }
    return u;
  }

  vector<coordinate_type>
  integrateNormalDerivativeBarnesHut(vector<coordinate_type> &coordinates) {
    mMesh->update_all();
    vector<face_ptr> &faces = mMesh->get_faces();

    vector<coordinate_type> chargeCenters;
    vector<coordinate_type> charges;
    vector<T> chargeMags;
    for (int i = 0; i < faces.size(); i++) {
      T area = faces[i]->area();

      coordinate_type norm = faces[i]->normal();
      norm.normalize();
      // if(area < 1e-9) continue;
      coordinate_type charge = area * norm;
      chargeCenters.push_back(faces[i]->center());
      charges.push_back(charge);
      chargeMags.push_back(area);
    }

    typedef pole_tree<SPACE> Tree;
    typedef pole_node<SPACE> Node;
    Tree octree(chargeCenters);

    vector<coordinate_type> nodeCharges;
    vector<coordinate_type> nodeChargeCenters;
    nodeCharges.resize(octree.nodes.size());
    nodeChargeCenters.resize(octree.nodes.size());
    std::stack<int> stack;
    stack.push(0);

    while (stack.size() > 0) {
      int pId = stack.top();
      stack.pop();
      Node &pNode = octree.nodes[pId];
      coordinate_type chargeCenter(0, 0, 0);
      coordinate_type chargeNet = 0;
      T netChargeMag = 0;
      int N = pNode.size;
      T iN = 1.0 / (T)N;
      int beg = pNode.begin;
      for (int i = beg; i < beg + N; i++) {
        int ii = octree.permutation[i];
        T chargeMag = chargeMags[ii];
        // T chargeMag = 1.0;
        coordinate_type chargeLoc = chargeCenters[ii];
        coordinate_type charge = charges[ii];
        netChargeMag += chargeMag;
        chargeCenter += chargeMag * chargeLoc;
        chargeNet += charges[ii];
      }
      chargeCenter /= netChargeMag;
      nodeCharges[pId] = chargeNet;
      nodeChargeCenters[pId] = chargeCenter;
      for (int j = 0; j < 8; j++) {
        if (pNode.children[j] != -1)
          stack.push(pNode.children[j]);
      }
    }
    T thresh = 0.5;
    vector<coordinate_type> u;
    vector<T> uSum;
    u.resize(coordinates.size(), coordinate_type(0, 0, 0));
    uSum.resize(coordinates.size(), 0.0);
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    for (int i = 0; i < coordinates.size(); i++) {
      int count = 0;
      // if(!verts[i]) continue;
      coordinate_type pi = coordinates[i];

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int pId = stack1.top();
        stack1.pop();
        Node &pNode = octree.nodes[pId];
        coordinate_type pij = nodeChargeCenters[pId];
        coordinate_type dpi = pi - pij;
        T dc = norm(dpi);
        T sc = norm(pNode.half);

        if (sc / dc > 0.5) {
          for (int j = 0; j < 8; j++) {
            if (pNode.children[j] != -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }

        else {

          int ii = octree.permutation[pNode.begin];
          coordinate_type ci = nodeCharges[pId];

          // coordinate_type pj = pNode.centerOfMass;
          if (norm(dpi) < 1e-16)
            continue;
          // coordinate_type ciN = ci; ciN.normalize();
          // coordinate_type dpp = va::dot(dp,ciN)/va::dot(dp,dp)*dp;
          T dist = norm(dpi);

          T i4pi = 0.25 / M_PI;
          T regLength = 1.0 * minLength;
          T m = 2.0;
          T wi = 1.0 / powf(dist * dist + regLength * regLength, m);

          int count = 0;
          // if(!verts[i]) continue;
          coordinate_type pi = coordinates[i];

          std::stack<int> stack2;
          stack2.push(0);
          while (stack2.size() > 0) {

            int pJd = stack2.top();
            stack2.pop();
            Node &pNodej = octree.nodes[pJd];
            coordinate_type pjj = nodeChargeCenters[pJd];
            coordinate_type dpj = pi - pjj;
            T dcj = norm(dpj);
            T scj = norm(pNodej.half);

            if (scj / dcj > 0.5) {
              for (int j = 0; j < 8; j++) {
                if (pNodej.children[j] != -1) {
                  stack2.push(pNodej.children[j]);
                }
              }
            }

            else {
              T distj = norm(dpj);

              T i4pi = 0.25 / M_PI;
              T m = 2.0;
              T wj = 1.0 / powf(distj * distj + regLength * regLength, m);
              u[i] += wi * wj * ci;
            }
          }
          uSum[i] += wi;
        }
      }
    }
    for (int i = 0; i < u.size(); i++) {
      u[i] *= 1.0 / uSum[i] / uSum[i];
    }
    return u;
  }

  void integrateNormalsRK2(T dt) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<coordinate_type> y0;
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();
    y0.resize(verts.size());
    vector<T> edgeWeights;
    vector<T> vertexWeights;
    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;
      y0[i] = verts[i]->coordinate();
    }

    vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
    T norm2 = 0;
    vector<coordinate_type> y1(y0);
    for (int i = 0; i < u0.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;

      y1[i] += 0.5 * dt * u0[i];
      norm2 += va::dot(u0[i], u0[i]);
    }
    std::cout << " first stage norm: " << sqrt(norm2) << std::endl;
    vector<coordinate_type> u1 = integrateNormalsBarnesHut(y1);
    norm2 = 0;
    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;
      verts[i]->coordinate() += dt * u1[i];
      norm2 += va::dot(u1[i], u1[i]);
    }
    std::cout << " second stage norm: " << sqrt(norm2) << std::endl;
  }

  vector<coordinate_type>
  integrateNormals(T dt, vector<coordinate_type> &coordinates) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<coordinate_type> y0(coordinates);
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();

    vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
    T norm2 = 0;
    for (int i = 0; i < u0.size(); i++) {
      y0[i] += dt * u0[i];
      norm2 += va::dot(u0[i], u0[i]);
    }
    std::cout << " first stage norm: " << sqrt(norm2) << std::endl;
    return y0;
  }

  vector<coordinate_type>
  integrateNormals(vector<T> dt, vector<coordinate_type> &coordinates) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<coordinate_type> y0(coordinates);
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();

    vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
    T norm2 = 0;
    for (int i = 0; i < u0.size(); i++) {
      y0[i] += dt[i] * u0[i];
      norm2 += va::dot(u0[i], u0[i]);
    }
    std::cout << " first stage norm: " << sqrt(norm2) << std::endl;
    return y0;
  }

  vector<coordinate_type>
  integrateCharges(T dt, vector<coordinate_type> &coordinates) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<coordinate_type> y0(coordinates);
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();

    vector<coordinate_type> u0 = integrateOrientedChargesBarnesHut(y0);
    T norm2 = 0;
    for (int i = 0; i < u0.size(); i++) {
      y0[i] += dt * u0[i];
      norm2 += va::dot(u0[i], u0[i]);
    }
    std::cout << " first stage norm: " << sqrt(norm2) << std::endl;
    return y0;
  }

  vector<coordinate_type>
  integrateMorphing(T dt, vector<coordinate_type> &coordinates,
                    vector<coordinate_type> &normals) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<coordinate_type> n0 = mMesh->get_normals();

    vector<coordinate_type> y0(coordinates);
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();

    vector<coordinate_type> u0 = integrateMorphingCharges(y0, normals);
    // vector<coordinate_type> u1 = integrateChargesBarnesHut(y0);
    T norm2 = 0;
    for (int i = 0; i < u0.size(); i++) {
      y0[i] += dt * u0[i];
      norm2 += va::dot(u0[i], u0[i]);
    }
    std::cout << " first stage norm: " << sqrt(norm2) << std::endl;
    return y0;
  }

  vector<coordinate_type>
  integrateNormalsRK2(T dt, vector<coordinate_type> &coordinates) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<coordinate_type> y0(coordinates);
    vector<coordinate_type> y1(y0);
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();

    vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
    T norm2 = 0;
    for (int i = 0; i < u0.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;

      y0[i] += 0.5 * dt * u0[i];
      norm2 += va::dot(u0[i], u0[i]);
    }
    std::cout << " first stage norm: " << sqrt(norm2) << std::endl;
    vector<coordinate_type> u1 = integrateNormalsBarnesHut(y0);
    norm2 = 0;
    for (int i = 0; i < coordinates.size(); i++) {
      y1[i] += dt * u1[i];
      norm2 += va::dot(u1[i], u1[i]);
    }
    std::cout << " second stage norm: " << sqrt(norm2) << std::endl;
    return y1;
  }

  void hashFaces() {
    face_bin<SPACE> newHash(*mMesh, 2.0 * minLength);
    std::cout << "hashed faces: " << newHash.faces().size() << std::endl;
    mFaceHash = newHash;
  }

  // void staticAddHandles(){
  //   vector<vertex_ptr>& verts = mMesh->get_vertices();
  //   for (int i = 0; i < verts.size(); i++){
  //   }
  // }

  void randomize() {
    vertex_array &vl = mMesh->get_vertices();
    va_iterator itb = vl.begin();
    va_iterator ite = vl.end();
    while (itb != ite) {
      (*itb)->coordinate()[0] = randd(0.1);
      (*itb)->coordinate()[1] = randd(0.1);
      (*itb)->coordinate()[2] = randd(0.1);

      itb++;
    }
  }

public:
  // aabb_tree<SPACE,2> & aabb(){return mAaBb;}
  face_bin<SPACE> &faceHash() { return mFaceHash; }
  face_bin<SPACE> mFaceHash;
  // aabb_tree<SPACE,2> mAaBb;
  vector<coordinate_type> mDirectionWeights;
  vector<coordinate_type> mDebug;
  vector<coordinate_type *> mDirectionField;
  T maxCurvature, minCurvature, minLength, minCollapseLength, maxLength;
  control_ptr mMesh;
};

template <typename SPACE> class vortex_sheet {
  M2_TYPEDEFS;

public:
  vortex_sheet(control_ptr obj_in) {

    mMesh = obj_in;
    maxCurvature = 3.0;
    minCurvature = 0.01;
    minLength = 0.03;
    regLength = 0.03;
    minCollapseLength = 0.0001;
    maxLength = 0.0005;
    edgeJoinThresh = 0.00025;
    m2::remesh<SPACE> rem;
    coordinate_type cen = mMesh->calc_center();
    coordinate_type min = mMesh->calc_min();
    coordinate_type max = mMesh->calc_max();
    coordinate_type lengths = max - min;
    int minRes = 64;
    T minLength = lengths[0];
    minLength = minLength < lengths[1] ? minLength : lengths[1];
    minLength = minLength < lengths[2] ? minLength : lengths[2];

    T dx = minLength / (T)(minRes - 1);

    for (int i = 0; i < 1; i++) {
      this->remesh();
      set_operations<SPACE> setOps;
      bool relaxing = true;
      while (relaxing)
        relaxing = setOps.flip_edges(mMesh);
    }

    vector<face_ptr> &faces = mMesh->get_faces();
    for (int i = 0; i < faces.size(); i++) {
      faces[i]->data = coordinate_type(0, 0, 0);
      faces[i]->data2 = coordinate_type(0, 0, 0);
    }

    // pin_bottom();
    // pin_half();

    this->hashFaces();
    // rem.triangulate(mMesh);
    // construct<SPACE> cons;
    // face_ptr nf = cons.delete_vertex(mMesh,mMesh->get_vertices()[0]);
    // vector<T> edgeWeights;
    // vector<T> vertexWeights;
    // calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);
  }

  coordinate_type getButterflyVertexWeight(face_vertex_ptr v) {
    // \/\./
    // vn = next
    // vnn = vertext_next next
    // pvpp = prev vertex_prev next
    // coordinate_type c = v->coordinate();
    // coordinate_type cn = v->next()->coordinate();
    // coordinate_type cvnn = v->vnext()->next()->coordinate();
    // coordinate_type cpvpp = v->prev()->vprev()->prev()->coordinate();
    // return (8.0*c + 2.0*cn - 1.0*(cvnn + cpvpp));
    face_vertex_ptr itb = v;
    face_vertex_ptr ite = v->vprev();
    bool iterating = true;
    int K = 0, j = 0;
    while (iterating) {
      iterating = itb != ite;
      itb = itb->vnext();
      K++;
    }
    itb = v;
    iterating = true;
    T s = 0;
    coordinate_type cs(0, 0, 0);
    while (iterating) {
      iterating = itb != ite;
      itb = itb->vnext();
      T frac = (T)j / (T)K;
      T sj =
          (0.25 + cos(2.0 * M_PI * frac) + 0.5 * cos(4.0 * M_PI * frac)) / (T)K;
      coordinate_type ci = itb->next()->coordinate();
      cs += sj * ci;
      s += sj;
      j++;
    }
    coordinate_type c = v->coordinate();
    return 0.75 * c + cs;
  }

  coordinate_type getButterflyWeight(edge_ptr e) {
    face_vertex_ptr v1 = e->v1();
    face_vertex_ptr v2 = e->v2();
    coordinate_type c1 = getButterflyVertexWeight(v1);
    coordinate_type c2 = getButterflyVertexWeight(v2);
    return 0.5 * (c1 + c2);
  }

  void flip_edges() {
    int k = 0;
    bool relaxing = true;
    while (relaxing && k < 3) {
      relaxing = set_operations<SPACE>::flip_edges(mMesh);
      k++;
    }
  }

  void remesh() {
    bool relaxing = true;
    int k = 0;
    // TIMER function//TIMER(__FUNCTION__);
    std::cout << " mesh size: " << mMesh->get_vertices().size() << std::endl;
    // relaxing = this->delete_vertices_curve();

    mMesh->colorVertices();
    std::cout << " max graph color: " << mMesh->maxGraphColor << std::endl;

    relaxing = true;
    while (relaxing) {
      relaxing = delete_degenerate();
    }

    relaxing = this->collapse_edges();
    relaxing = true;

    while (relaxing) {
      relaxing = delete_degenerate();
    }

    relaxing = true;
    while (relaxing) {
      relaxing = this->insert_edges_curve();
    }

    relaxing = true;
    while (relaxing) {
      relaxing = delete_degenerate();
    }
    // surface_filter<SPACE> filter;
    // for (int i = 0; i < 5; i++) {
    //  mMesh->update_all();
    //  surface_filter<SPACE>::filter(*mMesh, 0.01);
    //}

    // for (int i = 0; i < 1; i++)
    //	surface_filter<SPACE>::filter(*mMesh, 0.01);

    // this->joinEdges();
    this->flip_edges();
    // this->hashFaces();
    surface_calculator<SPACE> calc;
    // calc.shadeVerticesWinding(*mMesh);
    calc.shadeVertices(*mMesh);
  }

  // face_bin<SPACE> & getHash(){return mFaceHash;}

  pole_tree<SPACE> &octree() { return mOctree; }
  // aabb_tree<SPACE,3> & aabb(){return mAaBb;}
  vertex_bin<SPACE> &vertexHash() { return mVertexHash; }
  edge_bin<SPACE> &edgeHash() { return mEdgeHash; }
  face_bin<SPACE> &faceHash() { return mFaceHash; }

  void hashVertices() {
    vertex_bin<SPACE> newHash(*mMesh, 1.0 * regLength);
    mVertexHash = newHash;
  }

  void hashEdges() {
    edge_bin<SPACE> newHash(*mMesh, 1.0 * regLength);
    mEdgeHash = newHash;
  }

  void hashFaces() {
    face_bin<SPACE> newHash(*mMesh, 2.0 * minLength);
    mFaceHash = newHash;
  }

  void pin_half() {
    mMesh->pack();
    vertex_array &vl = mMesh->get_vertices();
    coordinate_type cen = mMesh->calc_center();
    coordinate_type min = mMesh->calc_min();
    coordinate_type max = mMesh->calc_max();
    for (int j = 0; j < vl.size(); j++) {
      vl[j]->pinned = false;
    }

    std::cout << " cen: " << cen << std::endl;
    for (int j = 0; j < vl.size(); j++) {
      if (vl[j]->coordinate()[1] < min[1] + 0.25 * (max[1] - min[1])) {
        vl[j]->pinned = true;
      }
    }
  }

  void pin_bottom() {
    mMesh->pack();
    vertex_array &vl = mMesh->get_vertices();
    coordinate_type cen = mMesh->calc_center();
    coordinate_type min = mMesh->calc_min();
    coordinate_type max = mMesh->calc_max();
    for (int j = 0; j < vl.size(); j++) {
      vl[j]->pinned = false;
    }

    std::cout << " cen: " << cen << std::endl;
    for (int j = 0; j < vl.size(); j++) {
      if (vl[j]->coordinate()[1] < min[1] + 0.05 * (max[1] - min[1])) {
        vl[j]->pinned = true;
      }
    }
  }

  void pin_bottom_facing() {
    mMesh->pack();
    vertex_array &vl = mMesh->get_vertices();
    coordinate_type cen = mMesh->calc_center();
    coordinate_type min = mMesh->calc_min();
    coordinate_type max = mMesh->calc_max();
    for (int j = 0; j < vl.size(); j++) {
      vl[j]->pinned = false;
    }

    std::cout << " cen: " << cen << std::endl;
    for (int j = 0; j < vl.size(); j++) {
      coordinate_type N = vl[j]->normal();
      T downAngle = va::dot(N, coordinate_type(0, -1.0, 0));
      T bottomThresh = min[1] + 0.2 * (max[1] - min[1]);
      if (vl[j]->coordinate()[1] < bottomThresh && downAngle > -0.5) {
        vl[j]->pinned = true;
      }
    }
  }

  void updateVorticity() {
    T normt = 0;
    T normx = 0;
    vector<face_ptr> &faces = mMesh->get_faces();

    for (int i = 0; i < faces.size(); i++) {
      if (!faces[i])
        continue;
      circulationToVorticity(faces[i]);
      T mag = va::norm2(faces[i]->data);
      normt += mag;
      normx = max(normx, sqrt(mag));
    }
    mMesh->verify();
    std::cout << " vorticity 1 norm: " << sqrt(normt)
              << " max magnitude: " << normx << std::endl;
  }

  void updateCirculation() {
    // TIMER function//TIMER(__FUNCTION__);
    T norm = 0;
    vector<face_ptr> &faces = mMesh->get_faces();
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    for (int i = 0; i < faces.size(); i++) {
      vorticityToCirculation(faces[i]);
      norm += va::dot(faces[i]->data, faces[i]->data);
    }
    std::cout << " circulation update norm: " << sqrt(norm) << std::endl;
  }

  void integrateBaroclinity(T dt, T C) {
    // TIMER function//TIMER(__FUNCTION__);
    mMesh->update_all();
    vector<face_ptr> &faces = mMesh->get_faces();
    T normt = 0;
    T normx = 0;
    coordinate_type g(0, -9.8, 0);
    for (int i = 0; i < faces.size(); i++) {
      if (!faces[i])
        continue;
      if (faces[i]->fbegin()->vertex()->pinned)
        continue;
      if (faces[i]->fbegin()->next()->vertex()->pinned)
        continue;
      if (faces[i]->fbegin()->next()->next()->vertex()->pinned)
        continue;
      T a = faces[i]->calc_area();
      if (a > 1e-16) {
        // T Ba  = 2.0*faces[i]->data3*0.0005;
        T Ba = 2.0 * C;
        faces[i]->update_normal();
        coordinate_type normal = faces[i]->normal();
        coordinate_type dB = -dt * Ba * va::cross(g, normal);
        faces[i]->data += dB;
        faces[i]->data2 = dB;
        T mag = va::norm2(faces[i]->data);
        normx = max(mag, sqrt(mag));
        normt += mag;
      }
    }
    std::cout << " vorticity 0 norm: " << sqrt(normt)
              << " max magnitude: " << normx << std::endl;
  }

  void integrateVelocityEuler(T dt) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    // vector<coordinate_type> u0 = integrateVelocity(0.5*dt);
    vector<coordinate_type> u0 = integrateVelocityTreeCode();
    // vector<coordinate_type> u0 = integrateVelocityBruteForce(0.5*dt);
    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;

      verts[i]->coordinate() += dt * u0[i];
      verts[i]->data = u0[i];
    }
  }

  void integrateVelocityRK2(T dt) {
    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<coordinate_type> y0;
    y0.resize(verts.size());
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();

    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;
      y0[i] = verts[i]->coordinate();
    }
    // vector<coordinate_type> uheat = integrateHeat(0.5*dt);

    vector<coordinate_type> u0 = integrateVelocityTreeCode();

    T normt = 0;
    T normx = 0;
    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;

      verts[i]->coordinate() += 0.5 * dt * (u0[i]);
      verts[i]->verify();
      // verts[i]->data = u0[i] + uheat[i];
      T mag = va::norm2(u0[i]);
      normx = max(normx, sqrt(mag));
      normt += mag;
    }
    std::cout << " velocityRK2 0 norm: " << sqrt(normt)
              << " max velocity mag: " << normx << std::endl;
    normt = 0;
    normx = 0;
    vector<coordinate_type> u1 = integrateVelocityTreeCode();
    // vector<coordinate_type> u1 = integrateVelocityBarnesHut();
    // vector<coordinate_type> u1 = integrateVelocityBruteForce(0.5*dt);

    for (int i = 0; i < verts.size(); i++) {
      if (!verts[i])
        continue;
      if (verts[i]->pinned)
        continue;

      coordinate_type du = dt * (u1[i]);
      verts[i]->coordinate() = y0[i] + du;
      verts[i]->verify();
      // verts[i]->data = u1[i] + uheat[i];
      T mag = va::norm2(u1[i]);
      normx = max(normx, sqrt(mag));
      normt += mag;
    }
    std::cout << " velocityRK2 1 norm: " << sqrt(normt)
              << " max velocity mag: " << normx << std::endl;
  }

  vector<coordinate_type> integrateVelocityBruteForce(T dt) {
    vector<coordinate_type> u;
    u.resize(mMesh->get_vertices().size(), coordinate_type(0, 0, 0));

    vector<face_ptr> &faces = mMesh->get_faces();
    vector<coordinate_type> chargeCenters;
    vector<coordinate_type> charges;
    vector<T> chargeMags;

    for (int i = 0; i < faces.size(); i++) {
      T area = faces[i]->area();
      coordinate_type vort = faces[i]->data;
      coordinate_type charge = area * vort;
      chargeCenters.push_back(faces[i]->center());
      charges.push_back(charge);
      chargeMags.push_back(norm(charge));
    }

    vector<vertex_ptr> &verts = mMesh->get_vertices();
    for (int i = 0; i < verts.size(); i++) {
      for (int j = 0; j < charges.size(); j++) {
        if (!verts[i])
          continue;
        coordinate_type pi = verts[i]->coordinate();
        coordinate_type pj = chargeCenters[j];

        coordinate_type ci = charges[j];
        // coordinate_type pi = chargeCenters[ii];
        coordinate_type dp = pi - pj;
        if (norm(dp) < 1e-12)
          continue;
        T dist = norm(dp);
        T i4pi = 0.25 / M_PI;

        // T dist3 = dist*dist*dist;
        // T l3   = regLength*regLength*regLength;
        // T kappa = (1.0-exp(-dist3/l3))/dist3;

        T l2 = regLength * regLength;
        T denom = powf((dist * dist + l2), 1.5);
        coordinate_type fR = dp / denom;
        // coordinate_type fR = dp*kappa;
        coordinate_type ui = i4pi * va::cross(ci, fR);
        // std::cout << ui << std::endl;
        u[i] += ui;
      }
    }
    return u;
  }

  vector<coordinate_type> integrateVelocityTreeCode() {

    auto pre = [](const vector<coordinate_type> &charges,
                  const vector<T> &weights,
                  const vector<coordinate_type> &points, int begin, int N,
                  const vector<int> &permutation, coordinate_type &netCharge,
                  coordinate_type &avgPoint) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      netCharge = coordinate_type(0, 0, 0);

      T netWeight = 0;

      for (int i = begin; i < begin + N; i++) {
        int ii = permutation[i];
        T w = weights[ii] + .001;
        // T chargeMag = 1.0;

        avgPoint += w * points[ii];
        netWeight += w;

        netCharge += charges[ii];
      }

      avgPoint /= netWeight;
    };

    auto compute = [this](const coordinate_type &charge,
                          const coordinate_type &chargePoint,
                          const coordinate_type &evalPoint) -> coordinate_type {
      coordinate_type dp = evalPoint - chargePoint;
      T dist = va::norm(dp);

      if (dist < regLength)
        return coordinate_type(0.0, 0.0, 0.0);
      T i4pi = 0.25 / M_PI;
#if 1
      T dist3 = dist * dist * dist;
      T l3 = regLength * regLength * regLength;
      T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
      coordinate_type fR = dp * kappa;
#else
      T dist2 = dist * dist;
      T l2 = regLength * regLength;
      T denom = powf((dist2 + l2), 1.5);
      coordinate_type fR = dp / denom;
#endif
      coordinate_type u = i4pi * va::cross(charge, fR);

      return dist3 > 0 ? u : coordinate_type(0, 0, 0);
    };

    vector<face_ptr> &faces = mMesh->get_faces();
    vector<coordinate_type> vorticityPositions;
    vector<coordinate_type> vorticity;
    vector<T> weights;
    // TIMER build//TIMER("building tree hierarchy");
    for (int i = 0; i < faces.size(); i++) {
      T area = faces[i]->area();
      // assert(area != 0);
      coordinate_type vort = faces[i]->data;
      coordinate_type charge = area * vort;
      vorticityPositions.push_back(faces[i]->center());
      vorticity.push_back(charge);
      weights.push_back(va::norm(charge));
    }
    vector<coordinate_type> vertexCoordinates = mMesh->get_coordinates();
    Simple_BarnesHutt<SPACE, coordinate_type, coordinate_type> integrator;

    std::vector<coordinate_type> u =
        integrator.integrate(vorticity, weights, vorticityPositions,
                             vertexCoordinates, pre, compute);

    return u;
  }

  void circulationToVorticity(face_ptr f) {
    // TIMER function//TIMER(__FUNCTION__);
    T a = f->calc_area();
    if (a > 1e-16) {
      face_vertex_ptr fv1 = f->fbegin();
      face_vertex_ptr fv2 = fv1->next();
      face_vertex_ptr fv3 = fv2->next();

      coordinate_type e13 = fv1->edge()->vec();
      coordinate_type e21 = fv2->edge()->vec();
      coordinate_type e32 = fv3->edge()->vec();

      T e1 = fv1->data;
      T e2 = fv2->data;
      T e3 = fv3->data;
      f->data = 1.0 / a * (e1 * e13 + e2 * e21 + e3 * e32);
      // Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
      // Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";"); std::cout
      // << f->data.format(CommaInitFmt) << std::endl;
      assert(a > 0);
    } else {
      f->data = coordinate_type(0, 0, 0);
    }
  }

  void vorticityToCirculation(face_ptr f) {
    // TIMER function//TIMER(__FUNCTION__);
    T a = f->calc_area();
    if (a > 1e-16) {
      face_vertex_ptr fv1 = f->fbegin();
      face_vertex_ptr fv2 = fv1->next();
      face_vertex_ptr fv3 = fv2->next();

      coordinate_type e13 = fv1->edge()->vec();
      coordinate_type e21 = fv2->edge()->vec();
      coordinate_type e32 = fv3->edge()->vec();

      coordinate_type gamma = a * f->data2;
      if (va::norm(e13) < 1e-12 || va::norm(e21) < 1e-12 ||
          va::norm(e32) < 1e-12) {
        fv1->data = 0.0;
        fv2->data = 0.0;
        fv3->data = 0.0;
        return;
      }

      Eigen::MatrixXd E(4, 3);
      E = Eigen::MatrixXd::Zero(4, 3);
      for (int i = 0; i < 3; i++) {
        E(i, 0) = e13[i];
        E(i, 1) = e21[i];
        E(i, 2) = e32[i];
      }
      E(3, 0) = 1;
      E(3, 1) = 1;
      E(3, 2) = 1;
      // Eigen::MatrixXd EtE(3, 3);
      // EtE = E.transpose() * E;
      Eigen::VectorXd b(4);
      b(0) = gamma[0];
      b(1) = gamma[1];
      b(2) = gamma[2];
      b(3) = 0;
      // b = E.transpose() * b;

      // Eigen::VectorXd c = EtE.ldlt().solve(b);
      Eigen::VectorXd c =
          E.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
      /*
     if(c.norm() > 0.01){
       std::cout << " a: " << a << std::endl;
       std::cout << " b: " << b.transpose() << std::endl;
       std::cout << " e: " << std::endl;
       std::cout << E << std::endl;
       std::cout << " c: " << c.norm() << std::endl;
       std::cout << " =======" << std::endl;
     }
     */
      m2::construct<SPACE> cons;
      if (c.squaredNorm() > 100.00) {
        fv1->data = 0.0;
        fv2->data = 0.0;
        fv3->data = 0.0;
        return;
      }

      fv1->data += c(0);
      fv2->data += c(1);
      fv3->data += c(2);
    }
  }

  edge_ptr subdivideFace(face_vertex_ptr fv) {

    // coordinate_type data = fv->face()->data;
    // coordinate_type data2 = fv->face()->data2;
    T data3 = fv->face()->data3;

    // TIMER function//TIMER(__FUNCTION__);
    // assuming that the face vertex is the newly inserted one.
    face_vertex_ptr fv1 = fv; // this is the
    face_vertex_ptr fv2 = fv->next();
    face_vertex_ptr fv3 = fv2->next();
    face_vertex_ptr fv4 = fv3->next();
#if 1
    coordinate_type e21 =
        fv2->vertex()->coordinate() - fv1->vertex()->coordinate();
    coordinate_type e32 =
        fv3->vertex()->coordinate() - fv2->vertex()->coordinate();
    coordinate_type e43 =
        fv4->vertex()->coordinate() - fv3->vertex()->coordinate();
    coordinate_type e14 =
        fv1->vertex()->coordinate() - fv4->vertex()->coordinate();
    coordinate_type e13 =
        fv1->vertex()->coordinate() - fv3->vertex()->coordinate();

    T e1 = fv1->data;
    T e2 = fv2->data;
    T e3 = fv3->data;
    T e4 = fv4->data;
    coordinate_type gamma0 = e2 * e21 + e3 * e32 + e4 * e43 + e1 * e14;
    Eigen::MatrixXd E(5, 6);
    E = Eigen::MatrixXd::Zero(5, 6);
    for (int i = 0; i < 3; i++) {
      E(i, 0) = e21[i];
      E(i, 1) = e13[i];
      E(i, 2) = e32[i];
      E(i, 3) = e14[i];
      E(i, 4) = e43[i];
      E(i, 5) = -e13[i];
    }
    E(3, 0) = 1;
    E(3, 1) = 1;
    E(3, 2) = 1;
    E(4, 3) = 1;
    E(4, 4) = 1;
    E(4, 5) = 1;
    Eigen::VectorXd b(5);
    b(0) = gamma0[0];
    b(1) = gamma0[1];
    b(2) = gamma0[2];
    b(3) = 0;
    b(4) = 0;
    // b = E.transpose()*b;
    // Eigen::MatrixXd EtE(3,3);
    // EtE = E.transpose()*E;
    // Eigen::VectorXd g = EtE.ldlt().solve(b);
    Eigen::VectorXd g =
        E.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

    m2::construct<SPACE> cons;
    edge_ptr enew = cons.insert_edge(mMesh, fv1, fv3);
    enew->v1()->data = g(1);
    enew->v2()->data = g(5);
    fv3->data = g(0);
    fv2->data = g(2);
    fv1->data = g(3);
    fv4->data = g(4);
#else
    m2::construct<SPACE> cons;
    edge_ptr enew = cons.insert_edge(mMesh, fv1, fv3);

    // enew->v1()->face()->data = e1;
    // enew->v2()->face()->data = e2;

    // enew->v1()->face()->data = e3;
    // enew->v2()->face()->data = e4;
    enew->v1()->data = 0.0;
    enew->v2()->data = 0.0;
#endif
    return enew;
  }

  void split_edge(edge_ptr e) {
    // TIMER function//TIMER(__FUNCTION__);
    face_vertex_ptr fv1 = e->v1();
    face_vertex_ptr fv2 = e->v2();
    T circ1 = fv1->data;
    T circ2 = fv2->data;
    construct<SPACE> cons;
    // coordinate_type c = getButterflyWeight(e);
    vertex_ptr nv = cons.subdivide_edge(mMesh, e);

    // nv->coordinate() = c;
    fv1->data = circ2;
    fv2->data = circ1;
    fv1->next()->data = circ2;
    fv2->next()->data = circ1;

    nv->winding =
        0.5 * (e->v1()->vertex()->winding + e->v2()->vertex()->winding);
    subdivideFace(fv1->next());
    subdivideFace(fv2->next());
  }

  struct edge_sort {
    bool operator()(edge_ptr ei, edge_ptr ej) {
      return (ei->length() < ej->length());
    }
  } mEdgeSorter;

  bool insert_edges_curve() {
    // TIMER function//TIMER(__FUNCTION__);
    edge_array &edges = mMesh->get_edges();
    vector<vertex_ptr> &verts = mMesh->get_vertices();

    face_array &faces = mMesh->get_faces();
    int N = faces.size();

    bool topology_change = false;
    vector<T> edgeWeights;
    vector<T> vertexWeights;
    mMesh->update_all();
    mMesh->reset_flags();
    mMesh->pack();
    surface_calculator<SPACE> calc;

    vector<edge_ptr> edgesToSplit;
#if 1
    for (int i = 0; i < edges.size(); i++) {
      if (!mMesh->has_edge(i))
        continue;
      edge_ptr ei = edges[i];
      bool pinned = ei->v1()->vertex()->pinned == true &&
                    ei->v2()->vertex()->pinned == true;
      if (pinned)
        continue;
      T l = ei->length();
      if (l > 1.75 * minLength) {
        edgesToSplit.push_back(ei);
        continue;
      }
    }
#endif
#if 0
      calc.calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);

      for (int i = 0; i < verts.size(); i++){
          if (!mMesh->has_vertex(i)) continue;
          if (verts[i]->size() == 0) continue;
          vertex_ptr v = verts[i];
          T k = fabs(vertexWeights[i]);
          
          face_vertex_ptr itb = v->fbegin();
          face_vertex_ptr ite = v->fend();
          bool at_head = false;  int i = 0;
          while (!at_head && i < 40) {
            at_head = itb==ite;
            edge_ptr ei = itb->edge();
            T l = ei->length();
            bool pinned =
              ei->v1()->vertex()->pinned == true &&
              ei->v2()->vertex()->pinned == true;
            if (pinned) 	  itb = itb->vnext();
            else{
              if(l*k > 0.5)
                if(l > 1.0*minLength)
                  ei->flag = 1;	    
            }
            itb = itb->vnext();
            i++;
          }
      }

      for (int i = 0; i < edges.size(); i++){
      	if(!edges[i]) continue;
        if(edges[i]->flag == 1){
          edges[i]->flag = 0;
          edgesToSplit.push_back(edges[i]);
        }
      }

#endif

    std::cout << " - sorting " << edgesToSplit.size() << " edges, ";
    std::sort(edgesToSplit.begin(), edgesToSplit.end(), mEdgeSorter);
    std::cout << "splitting " << edgesToSplit.size() << " edges" << std::endl;
    for (int i = edgesToSplit.size(); i > 0; i--) {
      this->split_edge(edgesToSplit[i - 1]);
    }

    return topology_change;
  }

  bool collapse_edges() {
    // TIMER function//TIMER(__FUNCTION__);
    bool topology_change = false;
    edge_array collectedEdges;
    vector<edge_ptr> &edges = mMesh->get_edges();
    mMesh->reset_flags();
    for (int i = 0; i < edges.size(); i++) {

      if (!mMesh->has_edge(i))
        continue;
      edge_ptr e = edges[i];
      if (e->flag == 1)
        continue;
      if (e->v1()->vertex()->flag == 1)
        continue;
      if (e->v2()->vertex()->flag == 1)
        continue;
      if (e->v2()->vertex()->pinned)
        continue;
      T dist = e->dist();
      if (dist < minCollapseLength) {
        e->v1()->vertex()->flag = 1;
        e->v2()->vertex()->flag = 1;
        collectedEdges.push_back(e);
      }
    }

    std::cout << " - deleting: " << collectedEdges.size() << " Tiny edges"
              << std::endl;

    for (int i = 0; i < collectedEdges.size(); i++) {
      construct<SPACE> cons;
      if (!collectedEdges[i])
        continue;
      edge_ptr e = collectedEdges[i];
      coordinate_type avg = 0.5 * (e->v1()->vertex()->coordinate() +
                                   e->v2()->vertex()->coordinate());
      e->v1()->vertex()->coordinate() = avg;
      // face_ptr nf = cons.delete_vertex_primitive(mMesh,e->v2()->vertex());
      // cons.collapse_edge_primitive(mMesh,e);
      cons.collapse_edge(mMesh, e);
      topology_change = true;
    }
    if (topology_change) {
      m2::remesh<space3> rem;
      rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    return topology_change;
  }

  bool delete_degenerate() {

    // TIMER function//TIMER(__FUNCTION__);
    bool topology_change = false;
    vector<vertex_ptr> &verts = mMesh->get_vertices();
    vector<edge_ptr> &edges = mMesh->get_edges();
    int zeroEdges = 0;
    int twoFaces = 0;
    int flatVol = 0;
    vector<edge_ptr> edgesToDelete;
    for (int i = 0; i < edges.size(); i++) {
      if (mMesh->has_edge(i)) {
        edge_ptr e = edges[i];
        if (e->flag == 1)
          continue;
        int s1 = e->v1()->face()->size();
        int s2 = e->v2()->face()->size();
        T a1 = e->v1()->face()->calc_area();
        T a2 = e->v2()->face()->calc_area();

        if (a1 < 1e-12 || a2 < 1e-12) {
          e->flag = 1;
          edgesToDelete.push_back(e);
          continue;
        }

        if (s1 < 3 || s2 < 3) {
          twoFaces++;
          e->flag = 1;
          edgesToDelete.push_back(e);
          continue;
        }

        if (e->v1()->face() == e->v2()->face()) {
          edgesToDelete.push_back(e);
          e->flag = 1;
          continue;
        }

        if (e->v1()->vertex() == e->v2()->vertex()) {
          zeroEdges++;
          e->flag = 1;
          edgesToDelete.push_back(e);
          continue;
        };

        if (e->v1()->prev()->vertex() == e->v2()->prev()->vertex()) {
          flatVol++;
          edgesToDelete.push_back(e);
          e->flag = 1;
          continue;
        };
      }
    }

    topology_change = edgesToDelete.size() > 0;
    construct<SPACE> cons;
    for (int i = 0; i < edgesToDelete.size(); i++) {
      edge_ptr e = edgesToDelete[i];
      face_ptr nf = cons.delete_edge(mMesh, e);
    }
    mMesh->verify();
#if 0
      std::cout << " - deleted: " 
		<< zeroEdges << " zero edges and " 
		<< twoFaces << " two faces and "
		<< flatVol << " flat volumes."
		<< std::endl;
#endif

#if 1
    int lowVerts = 0;
    for (int i = 0; i < verts.size(); i++) {
      if (!mMesh->has_vertex(i))
        continue;
      if (verts[i]->pinned)
        continue;
      vertex_ptr v = verts[i];

      if (verts[i]->size() < 3) {
        construct<SPACE> cons;
        face_ptr nf = cons.delete_vertex_primitive(mMesh, v);
        topology_change = true;
        continue;
      }

      // if(verts[i]->size() == 2){
      //   //std::cout << "degenerate vert: " << verts[i]->size() << std::endl;
      //   int sz = verts[i]->size();
      //   construct<SPACE> cons;

      //   edge_ptr e = v->fbegin()->next()->edge();
      //   cons.delete_edge(mMesh,e);

      //   face_ptr nf = cons.delete_vertex_primitive(mMesh,v);
      //   //numVerts++;
      //   topology_change = true;
      //   continue;
      // }
    }
#endif
#if 0
        // if somehow a vertex has two edges to another vertex...
    for (int i = 0; i < verts.size(); i++){
      if(!mMesh->has_vertex(i)) continue;
      if(verts[i]->pinned) continue;
      face_vertex_ptr fvb = verts[i]->fbegin();
      face_vertex_ptr fve = verts[i]->fend();

      bool iterating = true;
      bool degenerateVert = false;
      for(fvb; iterating;  fvb = fvb->vnext()){
        iterating = fvb != fve;
        if(fvb->next()->vertex()->flag == 0) fvb->next()->vertex()->flag = 1;
        else if(fvb->next()->vertex()->flag == 1){
          degenerateVert = true;
        }
      }
      iterating = true;
      for(fvb; iterating;  fvb = fvb->vnext()){
        iterating = fvb != fve;
        fvb->next()->vertex()->flag = 0;
      }
      if(degenerateVert){
        construct<SPACE> cons;
        face_ptr nf = cons.delete_vertex_primitive(mMesh,verts[i]);
        topology_change = true;
      }
    }

    m2::remesh<space3> rem;
    rem.triangulate(mMesh);
    mMesh->verify();
#endif

    if (topology_change) {
      m2::remesh<space3> rem;
      // rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    mMesh->verify();
    return topology_change;
  }

  bool checkCompatibility(edge_ptr ei, edge_ptr ej) {
    int vertIds[4];
    vertIds[0] = ei->v1()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[2] =
    // ei->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

    vertIds[1] = ei->v2()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[5] =
    // ei->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

    vertIds[2] = ej->v1()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[8] =
    // ej->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

    vertIds[3] = ej->v2()->face()->fbegin()->vertex()->position_in_set();
    // vertIds[11] =
    // ej->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        if (i != j)
          if (vertIds[i] == vertIds[j])
            return false;
      }

    return true;
  };

  vector<edge_ptr> getPairList() {
    // TIMER function//TIMER(__FUNCTION__);

    mMesh->reset_flags();
    this->hashEdges();
    T dx = mEdgeHash.dx();
    int xRes = mEdgeHash.xRes();
    int yRes = mEdgeHash.yRes();
    int zRes = mEdgeHash.zRes();
    int sz = 8;
    std::cout << "   - edge pair query sphere size: " << sz << std::endl;
    int *ip = new int[3 * sz];
    int ind = 0;
    for (int k = 0; k < 2; k++)
      for (int j = 0; j < 2; j++)
        for (int i = 0; i < 2; i++) {
          int index = i + j * 2 + k * 4;
          ip[3 * index + 0] = i;
          ip[3 * index + 1] = j;
          ip[3 * index + 2] = k;
          ind++;
        }

    T imax = dx / xRes;
    vector<int> &binStart = mEdgeHash.binStart();
    vector<int> &binnedVerts = mEdgeHash.binnedEdges();
    vector<edge_ptr> &edges = mMesh->get_edges();
    vector<edge_ptr> edgePairs;

    std::cout << "   - getting edge pairs: " << edges.size() << std::endl;

    for (int z = 0; z < zRes; z++)
      for (int y = 0; y < yRes; y++)
        for (int x = 0; x < xRes; x++) {
          int index = x + y * xRes + z * xRes * yRes;
          int bStart0 = binStart[index];
          int bEnd0 = binStart[index + 1];

          for (int i = bStart0; i < bEnd0; i++) {
            int i0 = binnedVerts[i];
            int jm = -1;

            if (!mMesh->has_edge(i0))
              continue;
            if (edges[i0]->flag == 1)
              continue;
            bool breakLoop = false;
            edge_ptr ei = edges[i0];
            // if((ei->v1()->vertex()->winding + ei->v1()->vertex()->winding)
            // < 2.0*M_PI) continue;
            for (int k = 0; k < sz; k++) {
              // for(int k = 0; k < 0; k++){
              int xo = ip[3 * k + 0];
              int yo = ip[3 * k + 1];
              int zo = ip[3 * k + 2];

              if (x + xo < xRes && x + xo > -1 && y + yo < yRes &&
                  y + yo > -1 && z + zo < zRes && z + zo > -1) {
                int offsetIndex =
                    x + xo + (y + yo) * xRes + (z + zo) * xRes * yRes;
                int bStart1 = binStart[offsetIndex];
                int bEnd1 = binStart[offsetIndex + 1];
                for (int j = bStart1; j < bEnd1; j++) {
                  int j0 = binnedVerts[j];
                  if (j0 == i0)
                    continue;
                  if (!mMesh->has_edge(j0))
                    continue;
                  edge_ptr ej = edges[j0];

                  int Numi1 = ei->v1()->face()->size();
                  int Numi2 = ei->v2()->face()->size();
                  int Numj1 = ei->v1()->face()->size();
                  int Numj2 = ei->v2()->face()->size();

                  if (Numi1 != 3)
                    continue;
                  if (Numi2 != 3)
                    continue;
                  if (Numj1 != 3)
                    continue;
                  if (Numj2 != 3)
                    continue;

                  int setPos = ej->position_in_set();
                  // battery of degeneracy tests:

                  if (ej->flag == 1 || ei->flag == 1)
                    continue;
                  bool compatible = checkCompatibility(ei, ej);
                  if (!compatible)
                    continue;
                  if (ei->v1()->face()->calc_area() <
                      0.1 * minLength * minLength)
                    continue;
                  if (ei->v2()->face()->calc_area() <
                      0.1 * minLength * minLength)
                    continue;
                  if (ej->v1()->face()->calc_area() <
                      0.1 * minLength * minLength)
                    continue;
                  if (ej->v2()->face()->calc_area() <
                      0.1 * minLength * minLength)
                    continue;

                  distance_calculator<SPACE> dcalc;
                  T d = dcalc.calcEdgeEdgeDistance(ei, ej);
                  coordinate_type Ni = 0.5 * ei->v1()->face()->normal() +
                                       0.5 * ei->v2()->face()->normal();
                  coordinate_type Nj = 0.5 * ej->v1()->face()->normal() +
                                       0.5 * ej->v2()->face()->normal();
                  Ni.normalize();
                  Nj.normalize();
                  T dNij = va::dot(Ni, Nj);
                  coordinate_type Vi = 0.5 * ei->v1()->face()->data +
                                       0.5 * ei->v2()->face()->data;
                  coordinate_type Vj = 0.5 * ej->v1()->face()->data +
                                       0.5 * ej->v2()->face()->data;
                  T dVij = (Vi + Vj).norm() / (Vi).norm();
                  T avgWinding = 0.25 * (ei->v1()->vertex()->winding +
                                         ei->v1()->vertex()->winding +
                                         ej->v1()->vertex()->winding +
                                         ej->v1()->vertex()->winding);
                  //   avgWinding > 2.0*M_PI
                  if (d < edgeJoinThresh && dNij < -0.98) {
                    jm = j0;
                  }
                }
              }
            }
            if (jm > 0) {
              edge_ptr ej = edges[jm];
              edgePairs.push_back(ei);
              edgePairs.push_back(ej);
              ej->flag = 1;
              ei->flag = 1;
              ei->v1()->next()->edge()->flag = 1;
              ei->v1()->prev()->edge()->flag = 1;
              ei->v2()->next()->edge()->flag = 1;
              ei->v2()->prev()->edge()->flag = 1;

              ej->v1()->next()->edge()->flag = 1;
              ej->v1()->prev()->edge()->flag = 1;
              ej->v2()->next()->edge()->flag = 1;
              ej->v2()->prev()->edge()->flag = 1;
            }
          }
        }
    return edgePairs;
  }

  void pipeFace(face_ptr f0, face_ptr f1, T d00, T d01, T d10, T d11) {
    std::cout << "-- pipe face --" << std::endl;

    face_vertex_ptr v0b = f0->fbegin();
    face_vertex_ptr v0e = f0->fend();
    face_vertex_ptr v1b = f1->fbegin();
    face_vertex_ptr v1e = f1->fend();
    T mine = 9999;
    face_vertex_ptr v0s = v0b;
    face_vertex_ptr v1s = v1b;
    bool iterating0 = true;
    while (iterating0) {
      iterating0 = v0b != v0e;
      bool iterating1 = true;
      while (iterating1) {
        iterating1 = v1b != v1e;
        T d0 = (v0b->coordinate() - v1b->coordinate()).norm();
        T d1 = (v0b->next()->next()->coordinate() -
                v1b->next()->next()->coordinate())
                   .norm();
        T e = d0 * d0 + d1 * d1;
        if (e < mine) {
          mine = e;
          v0s = v0b;
          v1s = v1b;
        }
        v1b = v1b->next();
      }
      v0b = v0b->next();
    }
    iterating0 = true;
    v0b = v0s;
    v0e = v0s->prev();
    vector<face_vertex_ptr> pairs;
    vector<vertex_ptr> verts;
    while (iterating0) {
      iterating0 = v0b != v0e;
      if (v0s->vertex() != v1s->vertex()) {
        pairs.push_back(v0s);
        pairs.push_back(v1s);
        verts.push_back(v0s->vertex());
        verts.push_back(v1s->vertex());
      } else
        verts.push_back(v0s->vertex());
      v0s = v0s->next();
      v1s = v1s->prev();
      v0b = v0b->next();
    }

    // for (int i = 0; i < pairs.size(); i += 2) {}
    for (int i = 0; i < pairs.size(); i += 2) {
      construct<SPACE> cons;
      std::cout << "edges: " << pairs[i]->edge() << " " << pairs[i + 1]->edge()
                << std::endl;

      edge_ptr ei = cons.insert_edge(mMesh, pairs[i], pairs[i + 1]);

      ei->v1()->face()->color.g = 1.0;
      ei->v2()->face()->color.g = 1.0;
      if (i == 0)
        ei->v1()->data = d00;
      ei->v2()->data = d01;
      if (i == 2)
        ei->v1()->data = d10;
      ei->v2()->data = d11;
    }

    // for(int k = 0; k < 10; k++){
    // 	mesh_filter<SPACE> filter;
    // 	for(int i = 0; i < verts.size(); i++){
    // 	  verts[i]->coordinate() +=
    // 	    0.1*filter.laplacianFilterVertex(verts[i]);
    // 	}
    // }
    std::cout << std::endl;
  }

  bool joinEdges() {
    bool topology_change = false;

    std::cout << " - join edges " << std::endl;
    vector<edge_ptr> edgePairs = this->getPairList();
    if (edgePairs.size() > 0)
      topology_change = true;

    std::cout << "   - joining: " << edgePairs.size() << " pairs" << std::endl;
    for (int i = 0; i < edgePairs.size(); i += 2) {
      edge_ptr e0 = edgePairs[i + 0];
      edge_ptr e1 = edgePairs[i + 1];

      face_vertex_ptr v00 = e0->v2()->next();
      face_vertex_ptr v01 = e0->v1()->next();
      face_vertex_ptr v00n = v00->next();
      face_vertex_ptr v01n = v01->next();

      face_vertex_ptr v10 = e1->v2()->next();
      face_vertex_ptr v11 = e1->v1()->next();
      face_vertex_ptr v10n = v10->next();
      face_vertex_ptr v11n = v11->next();

      construct<SPACE> cons;
      T dt00 = e0->v1()->data;
      T dt01 = e0->v2()->data;
      T dt10 = e1->v1()->data;
      T dt11 = e1->v2()->data;
      v00->vertex()->color.g = 0.5;
      v01->vertex()->color.g = 0.5;
      v10->vertex()->color.g = 0.5;
      v11->vertex()->color.g = 0.5;
      if (v00->vertex()->pinned)
        continue;
      if (v01->vertex()->pinned)
        continue;
      if (v10->vertex()->pinned)
        continue;
      if (v11->vertex()->pinned)
        continue;
      face_ptr f0 = cons.delete_edge(mMesh, e0);
      face_ptr f1 = cons.delete_edge(mMesh, e1);
      // std::cout << f0->size() << " " << f1->size() << std::endl;
      pipeFace(f0, f1, dt00, dt01, dt10, dt11);
      break;
    }

    if (topology_change) {
      m2::remesh<space3> rem;
      rem.triangulate(mMesh);
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
    }
    return topology_change;
  }

  void drawVorticity() {
    face_array faces = mMesh->get_faces();
    for (int i = 0; i < faces.size(); i++) {
      face_ptr f = faces[i];
      if (f) {
        coordinate_type w = f->data;
        coordinate_type cen = f->center();
        glBegin(GL_LINES);
        T c = 0.005;
        glColor3f(0.0, 0.0, 1.0);
        glVertex3d(cen[0], cen[1], cen[2]);
        glVertex3d(cen[0] + c * w[0], cen[1] + c * w[1], cen[2] + c * w[2]);
        glEnd();
      }
    }
  }

  void drawVelocity() {
    vertex_array verts = mMesh->get_vertices();
    for (int i = 0; i < verts.size(); i++) {
      vertex_ptr v = verts[i];
      if (v) {
        coordinate_type vel = v->data;
        coordinate_type cen = v->coordinate();
        glBegin(GL_LINES);
        T c = 0.01;
        glColor3f(0.0, 1.0, 0.0);
        glVertex3d(cen[0], cen[1], cen[2]);
        glVertex3d(cen[0] + c * vel[0], cen[1] + c * vel[1],
                   cen[2] + c * vel[2]);
        glEnd();
      }
    }
  }

  edge_bin<SPACE> mEdgeHash;
  pole_tree<SPACE> mOctree;
  // aabb_tree<SPACE,3>  mAaBb;
  face_bin<SPACE> mFaceHash;
  vertex_bin<SPACE> mVertexHash;
  T maxCurvature, minCurvature, minLength, minCollapseLength, maxLength,
      regLength, edgeJoinThresh;
  vector<coordinate_type> mEdgeCirculation;
  vector<coordinate_type> mFaceVorticity;
  control_ptr mMesh;
};
} // namespace m2
#endif
