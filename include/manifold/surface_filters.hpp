#ifndef __M2SURFACE_FILTER__
#define __M2SURFACE_FILTER__

#include <stack>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bins.hpp"
#include "conj_grad.hpp"
#include "debugger.h"
#include "geometry_types.hpp"
#include "m2Includes.h"
#include "modify.hpp"
#include "remesh.hpp"
#include "surface_calculator.hpp"

namespace m2 {
template <typename SPACE> class surface_filter {
  M2_TYPEDEFS;

public:
  void filterCutoff(m2::surf<SPACE> &in, T cutoff, vector<T> &vertexWeights,
                    T strength) {
    // TIMER functionTimer(__FUNCTION__);
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    m2::surface_calculator<SPACE> calc;

    for (long i = 0; i < tverts.size(); i++) {

      if (!in.has_vertex(i))
        continue;
      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;

      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;
      coordinate_type kA(0, 0, 0, 0);
      T wTotal = 0, aTotal = 0;
      coordinate_type c0 = tverts[i]->coordinate();
      int k = 0;
      while (!at_head) {
        at_head = itb == ite;
        // mesh_calculator<SPACE> calc;
        T wij = calc.getEdgeWeight(itb->edge());
        // T wij = 1.0;
        coordinate_type c1 = itb->next()->coordinate();
        T aij = calc.baryArea(itb);

        aTotal += aij;
        wTotal += wij;

        kA += wij * (c1 - c0);
        itb = itb->vnext();
        k++;
      }
      T kT = vertexWeights[v->position_in_set()];

      if (wTotal > 1e-6 && kT < cutoff)
        filteredCoordinates[i] = kA / wTotal;
      else
        filteredCoordinates[i] = coordinate_type(0, 0, 0);
    }

    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      coordinate_type ci = tverts[i]->coordinate();
      tverts[i]->coordinate() = ci + strength * filteredCoordinates[i];
    }
  }

  static void filter(m2::surf<SPACE> &in, T strength) {
    // TIMER functionTimer(__FUNCTION__);
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    for (long i = 0; i < tverts.size(); i++) {

      if (!in.has_vertex(i))
        continue;
      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;

      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;

      coordinate_type kA(0, 0, 0);
      T wTotal = 0, aTotal = 0;
      coordinate_type c0 = tverts[i]->coordinate();
      int k = 0;
      while (!at_head) {
        at_head = itb == ite;
        surface_calculator<SPACE> calc;
        T wij = calc.getEdgeWeight(itb->edge());
        // T wij = 1.0;
        coordinate_type c1 = itb->next()->coordinate();
        T aij = calc.baryArea(itb);

        aTotal += aij;
        wTotal += wij;

        kA += wij * (c1 - c0);
        itb = itb->vnext();
        k++;
      }
      T kT = (kA / (2.0 * aTotal)).norm();

      if (wTotal > 1e-10)
        filteredCoordinates[i] = kA / wTotal;
      else
        filteredCoordinates[i] = coordinate_type(0, 0, 0);
    }

    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      coordinate_type ci = tverts[i]->coordinate();
      tverts[i]->coordinate() = ci + strength * filteredCoordinates[i];
    }
  }

  static void cuspFilter(m2::surf<SPACE> &in, T strength) {
    // TIMER functionTimer(__FUNCTION__);
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    for (long i = 0; i < tverts.size(); i++) {

      if (!in.has_vertex(i))
        continue;
      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;

      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;
      coordinate_type kA(0, 0, 0);
      T wTotal = 0, aTotal = 0;
      coordinate_type c0 = tverts[i]->coordinate();
      int k = 0;
      while (!at_head) {
        at_head = itb == ite;
        surface_calculator<SPACE> calc;
        T wij = calc.getEdgeWeight(itb->edge());
        edge_ptr e = itb->edge();
        coordinate_type n1 = e->v1()->face()->normal();
        coordinate_type n2 = e->v2()->face()->normal();
        wij = 1.0 + dot(n1, n2);
        // wij = dot(n1,n2);
        coordinate_type c1 = itb->next()->coordinate();

        wTotal += wij;

        kA += wij * (c1 - c0);
        itb = itb->vnext();
        k++;
      }
      T kT = norm(kA / (2.0 * aTotal));

      if (wTotal > 1e-10)
        filteredCoordinates[i] = kA / wTotal;
      else
        filteredCoordinates[i] = coordinate_type(0, 0, 0);
    }

    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      coordinate_type ci = tverts[i]->coordinate();
      tverts[i]->coordinate() = ci + strength * filteredCoordinates[i];
    }
  }

  static coordinate_type laplacianFilterVertex(vertex_ptr v) {
    if (v->size() == 0)
      return coordinate_type(0, 0, 0);
    face_vertex_ptr itb = v->fbegin();
    face_vertex_ptr ite = v->fend();
    bool at_head = false;
    coordinate_type kA(0, 0, 0);
    T wTotal = 0;
    coordinate_type c0 = v->coordinate();
    int k = 0;
    while (!at_head) {
      at_head = itb == ite;
      T wij = 1.0;
      coordinate_type c1 = itb->next()->coordinate();
      if (itb->next()->vertex()->flag == 0) {
        itb->next()->vertex()->flag += 1;
        wTotal += wij;
        kA += wij * (c1 - c0);
        k++;
      }
      itb = itb->vnext();
    }

    itb = v->fbegin();
    ite = v->fend();
    at_head = false;

    while (!at_head) {
      at_head = itb == ite;
      itb->next()->vertex()->flag = 0;
      itb = itb->vnext();
    }

    if (wTotal > 1e-10) {
      v->update_normal();
      coordinate_type w;
      coordinate_type dv = kA / wTotal;
      return dv;
    } else
      return coordinate_type(0, 0, 0);
  }

  static coordinate_type projectOntoLine(const coordinate_type &v0,
                                         const coordinate_type &v1,
                                         const coordinate_type &pt) {
    coordinate_type s = v1 - v0;
    coordinate_type v = pt - v0;
    coordinate_type vp = dot(v, s) / dot(s, s) * s;
    coordinate_type ptl = v0 + vp;
    return ptl;
  }

  static void mlsFilter(m2::surf<SPACE> &in, T strength, T dx) {
    // TIMER functionTimer(__FUNCTION__);
    vector<T> vertexWeights;
    vector<T> edgeWeights;
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;

      vertex_ptr v = tverts[i];
      coordinate_type c0 = tverts[i]->coordinate();
      if (v->size() == 0)
        continue;
      coordinate_type ca(0, 0, 0);
      T wTot = 0;
      m2::surface_calculator<SPACE> calc;
      vector<vertex_ptr> lverts = calc.getLocalVertices(in, v, 2);

      for (int j = 0; j < lverts.size(); j++) {
        coordinate_type c1 = lverts[j]->coordinate();
        coordinate_type dc = c1 - c0;
        T dist = norm(dc);
        // T wij = 1.0;
        T exponent = -dist * dist / (dx * dx);
        T wij = exp(exponent);
        wTot += wij;
        ca += wij * dc;
      }
      ca /= wTot;
      ca += c0;

      mat3 m;
      coordinate_type w;
      calc.calcCovariance(in, v, w, m, dx);
      coordinate_type d(m(0, 0), m(1, 0), m(2, 0));
      coordinate_type d2(m(0, 2), m(1, 2), m(2, 2));

      coordinate_type cp = projectOntoLine(ca - 0.1 * d, ca + 0.1 * d, c0);
      filteredCoordinates[i] = cp - c0;
    }

    for (long i = 0; i < filteredCoordinates.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      tverts[i]->coordinate() += strength * filteredCoordinates[i];
    }
  }

  static void laplacianFilter(m2::surf<SPACE> &in, T strength) {
    // TIMER functionTimer(__FUNCTION__);
    vector<T> vertexWeights;
    vector<T> edgeWeights;
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;

      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;
      filteredCoordinates[i] = laplacianFilterVertex(v);
    }

    for (long i = 0; i < filteredCoordinates.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      tverts[i]->coordinate() += strength * filteredCoordinates[i];
    }
  }

  static void taubinFilter(m2::surf<SPACE> &in, T a, T b) {
    // TIMER functionTimer(__FUNCTION__);
    vector<T> vertexWeights;
    vector<T> edgeWeights;
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;

      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;
      filteredCoordinates[i] = laplacianFilterVertex(v);
    }

    for (long i = 0; i < filteredCoordinates.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      tverts[i]->coordinate() += a * filteredCoordinates[i];
    }

    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;
      filteredCoordinates[i] = laplacianFilterVertex(v);
    }

    for (long i = 0; i < filteredCoordinates.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      tverts[i]->coordinate() -= b * filteredCoordinates[i];
    }
  }

  static void cacheTensor(m2::surf<SPACE> &in,
                          vector<coordinate_type *> &tensorArray,
                          vector<coordinate_type> &singularArray) {
    vector<vertex_ptr> &tverts = in.get_vertices();
    for (long i = 0; i < tverts.size(); i++) {

      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;
      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;
      coordinate_type *cov = new coordinate_type[3];
      while (!at_head) {
        at_head = itb == ite;
        T wij = itb->face()->area();
        coordinate_type N = itb->face()->normal();
        for (int l = 0; l < 3; l++)
          for (int m = 0; m < 3; m++) {
            cov[l][m] += wij * N[l] * N[m];
          }
        itb = itb->vnext();
      }
      coordinate_type w;
      calcSVD<SPACE>(cov, w);
      tensorArray[i] = cov;
      singularArray[i] = w;
    }
  }

  static coordinate_type nullLaplacianFilterVertex(vertex_ptr v,
                                                   coordinate_type cov[3],
                                                   coordinate_type w) {
    if (v->size() == 0)
      return coordinate_type(0, 0, 0);
    face_vertex_ptr itb = v->fbegin();
    face_vertex_ptr ite = v->fend();
    bool at_head = false;
    coordinate_type kA(0, 0, 0, 0);
    T wTotal = 0;
    coordinate_type c0 = v->coordinate();
    int k = 0;
    // coordinate_type cov[3];
    while (!at_head) {
      at_head = itb == ite;
      T wij = itb->face()->area();
      coordinate_type c1 = itb->next()->coordinate();
      wTotal += wij;
      kA += wij * (c1 - c0);
      itb = itb->vnext();
      k++;
    }
    coordinate_type dv = kA / wTotal;
    if (wTotal > 1e-10) {
      v->update_normal();

      coordinate_type u = cov[1];
      u.normalize();
      coordinate_type v = cov[2];
      v.normalize();
      coordinate_type outu = dot(dv, u) * u;
      coordinate_type outv = dot(dv, v) * v;

      T wp = w[1] + w[2];
      if (w[1] < 1e-3 && w[2] < 1e-3)
        return outu + outv;
      else if (w[2] < 1e-3)
        return outv;
      else
        return coordinate_type(0, 0, 0);
    } else
      return coordinate_type(0, 0, 0);
  }

  static void nullLaplacianFilter(m2::surf<SPACE> &in, T strength) {
    // TIMER functionTimer(__FUNCTION__);
    vector<T> vertexWeights;
    vector<T> edgeWeights;
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    vector<coordinate_type *> tensorArray;
    tensorArray.resize(tverts.size());
    vector<coordinate_type> singularArray;
    singularArray.resize(tverts.size());
    cacheTensor(in, tensorArray, singularArray);
    for (int k = 0; k < 2; k++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (long i = 0; i < tverts.size(); i++) {
        if (!in.has_vertex(i))
          continue;
        if (tverts[i]->pinned)
          continue;
        vertex_ptr v = tverts[i];
        if (v->size() == 0)
          continue;
        filteredCoordinates[i] =
            nullLaplacianFilterVertex(v, tensorArray[i], singularArray[i]);
      }

      for (long i = 0; i < filteredCoordinates.size(); i++) {
        if (!in.has_vertex(i))
          continue;
        if (tverts[i]->pinned)
          continue;
        tverts[i]->coordinate() += strength * filteredCoordinates[i];
      }
    }
    for (int i = 0; i < tensorArray.size(); i++) {
      delete tensorArray[i];
    }
  }

  void flaggedFilter(m2::surf<SPACE> &in, T strength) {
    // TIMER functionTimer(__FUNCTION__);
    vector<T> vertexWeights;
    vector<T> edgeWeights;
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<coordinate_type> filteredCoordinates;
    filteredCoordinates.resize(tverts.size());
    for (long i = 0; i < tverts.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      vertex_ptr v = tverts[i];
      if (v->size() == 0)
        continue;
      if (v->flag != 1)
        continue;
      // filteredCoordinates[i] = laplacianFilterVertex(v);
      filteredCoordinates[i] = nullLaplacianFilterVertex(v);
    }

    for (long i = 0; i < filteredCoordinates.size(); i++) {
      if (!in.has_vertex(i))
        continue;
      if (tverts[i]->pinned)
        continue;
      if (tverts[i]->flag != 1)
        continue;
      tverts[i]->coordinate() += strength * filteredCoordinates[i];
    }
  }
};
} // namespace m2
#endif
