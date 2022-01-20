#ifndef __M2MOVING__
#define __M2MOVING__

#include <stack>

//#ifdef _OPENMP
//# include <omp.h>
//#endif

#include "bins.hpp"
#include "debugger.h"
#include "m2Includes.h"
#include "remesh.hpp"

#include "tree_code.hpp"

#include "TIMER.h"

#include "surface_calculator.hpp"
#include "surface_filters.hpp"
#include <cmath>

namespace m2 {

template <typename SPACE>
class set_operations : public default_interface<SPACE> {
  M2_TYPEDEFS;

public:
  bool flip_edges(surf_ptr obj) {
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

      coordinate_type c0 = this->coordinate(v0);
      coordinate_type c1 = this->coordinate(v1);
      coordinate_type c2 = this->coordinate(v2);
      coordinate_type c3 = this->coordinate(v3);

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

template <typename SPACE> class vortex_sheet : public default_interface<SPACE> {
  M2_TYPEDEFS;

public:
  vortex_sheet(surf_ptr obj_in) {

    mMesh = obj_in;
    maxCurvature = 3.0;
    minCurvature = 0.01;
    minLength = 0.03;
    regLength = 0.03;
    minCollapseLength = 0.0001;
    maxLength = 0.0005;
    edgeJoinThresh = 0.00025;
    m2::remesh<SPACE> rem;

    box_type bb = this->bound(mMesh);
    coordinate_type lengths = 2.0 * bb.half;
    int minRes = 64;
    T minLength = lengths[0];
    minLength = minLength < lengths[1] ? minLength : lengths[1];
    minLength = minLength < lengths[2] ? minLength : lengths[2];

    T dx = minLength / (T)(minRes - 1);

    set_operations<SPACE> setOps;
    for (int i = 0; i < 1; i++) {
      this->remesh();
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
      coordinate_type ci = this->coordinate(itb->next());
      cs += sj * ci;
      s += sj;
      j++;
    }
    coordinate_type c = this->coordinate(v);
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
    set_operations<SPACE> setOps;
    while (relaxing && k < 3) {
      relaxing = setOps.flip_edges(mMesh);
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
      if (this->coordinate(vl[j])[1] < min[1] + 0.25 * (max[1] - min[1])) {
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
      if (this->coordinate(vl[j])[1] < min[1] + 0.05 * (max[1] - min[1])) {
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
      if (this->coordinate(vl[j])[1] < bottomThresh && downAngle > -0.5) {
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
      T a = this->area(faces[i]);
      if (a > 1e-16) {
        // T Ba  = 2.0*faces[i]->data3*0.0005;
        T Ba = 2.0 * C;

        coordinate_type n = this->normal(faces[i]);
        coordinate_type dB = -dt * Ba * va::cross(g, n);
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

      this->coordinate(verts[i]) += dt * u0[i];
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
      y0[i] = this->coordinate(verts[i]);
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

      coordinate_type ci = this->coordinate(verts[i]);
      this->coordinate(ci + 0.5 * dt * (u0[i]), verts[i]);
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
      this->coordinate(y0[i] + du, verts[i]);
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
        coordinate_type pi = this->coordinate(verts[i]);
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
      T a = this->area(faces[i]);
      // assert(area != 0);
      coordinate_type vort = faces[i]->data;
      coordinate_type charge = a * vort;
      vorticityPositions.push_back(this->center(faces[i]));
      vorticity.push_back(charge);
      weights.push_back(va::norm(charge));
    }
    vector<coordinate_type> vertexCoordinates = ci::get_coordinates<SPACE>(mMesh);
    Simple_BarnesHutt<SPACE, coordinate_type, coordinate_type> integrator;

    std::vector<coordinate_type> u =
        integrator.integrate(vorticity, weights, vorticityPositions,
                             vertexCoordinates, pre, compute);

    return u;
  }

  void circulationToVorticity(face_ptr f) {
    // TIMER function//TIMER(__FUNCTION__);
    T a = this->area(f);
    if (a > 1e-16) {
      face_vertex_ptr fv1 = f->fbegin();
      face_vertex_ptr fv2 = fv1->next();
      face_vertex_ptr fv3 = fv2->next();

      coordinate_type e13 = this->dir(fv1->edge());
      coordinate_type e21 = this->dir(fv2->edge());
      coordinate_type e32 = this->dir(fv3->edge());

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
    T a = this->area(f);
    if (a > 1e-16) {
      face_vertex_ptr fv1 = f->fbegin();
      face_vertex_ptr fv2 = fv1->next();
      face_vertex_ptr fv3 = fv2->next();

      coordinate_type e13 = this->dir(fv1->edge());
      coordinate_type e21 = this->dir(fv2->edge());
      coordinate_type e32 = this->dir(fv3->edge());

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
    coordinate_type e21 = this->coordinate(fv2) - this->coordinate(fv1);
    coordinate_type e32 = this->coordinate(fv3) - this->coordinate(fv2);
    coordinate_type e43 = this->coordinate(fv4) - this->coordinate(fv3);
    coordinate_type e14 = this->coordinate(fv1) - this->coordinate(fv4);
    coordinate_type e13 = this->coordinate(fv1) - this->coordinate(fv3);

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
      return (ci::length<SPACE>(ei) < ci::length<SPACE>(ej));
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
      T l = this->length(ei);
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
      T dist = this->length(e);
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
      coordinate_type avg = 0.5 * (this->coordinate(e->v1()->vertex()) +
                                   this->coordinate(e->v2()->vertex()));
      this->coordinate(avg, e->v1()->vertex());
      // face_ptr nf = cons.delete_vertex_primitive(mMesh,e->v2()->vertex());
      // cons.collapse_edge_primitive(mMesh,e);
      cons.collapse_edge(mMesh, e);
      topology_change = true;
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
        T a1 = this->area(e->v1()->face());
        T a2 = this->area(e->v2()->face());

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
    vertIds[0] = ei->v1()->face()->position_in_set();
    vertIds[1] = ei->v2()->face()->position_in_set();
    vertIds[2] = ej->v1()->face()->position_in_set();
    vertIds[3] = ej->v2()->face()->position_in_set();

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        if (i != j)
          if (vertIds[i] == vertIds[j])
            return false;
      }

    return true;
  };
  

  vector<edge_ptr> getEdgePairList() {
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
                  if (this->area(ei->v1()->face()) <
                      0.1 * minLength * minLength)
                    continue;
                  if (this->area(ei->v2()->face()) <
                      0.1 * minLength * minLength)
                    continue;
                  if (this->area(ej->v1()->face()) <
                      0.1 * minLength * minLength)
                    continue;
                  if (this->area(ej->v2()->face()) <
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
        T d0 = (this->coordinate(v0b) - this->coordinate(v1b)).norm();
        T d1 = (this->coordinate(v0b->next()->next()) -
                this->coordinate(v1b->next()->next()))
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
    vector<edge_ptr> edgePairs = this->getEdgePairList();
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

  vector<face_ptr> getFacePairList() { 
    vector<face_ptr> pairs;
    
    return pairs;
  }

  bool joinFaces() {
    bool topology_change = false;
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
        coordinate_type cen = this->coordinate(v);
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
  surf_ptr mMesh;
};
} // namespace m2
#endif
