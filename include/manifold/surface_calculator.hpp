#ifndef __M2SURFACE_CALCULATOR__
#define __M2SURFACE_CALCULATOR__

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

namespace m2 {

template <typename SPACE>
inline void calcSVD(typename SPACE::coordinate_type *vec,
                    typename SPACE::coordinate_type &val) {
  M2_TYPEDEFS;
  Eigen::Matrix3f m3;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      m3(i, j) = vec[i][j];
  Eigen::JacobiSVD<Eigen::Matrix3f> svd(m3, Eigen::ComputeFullU);

  const Eigen::Matrix3f U = svd.matrixU();
  const Eigen::VectorXf S = svd.singularValues();
  val = coordinate_type(S[0], S[1], S[2]);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      vec[j][i] = U(i, j);
}

template <typename SPACE>
inline void calcSVD(typename SPACE::mat3 &mi,
                    typename SPACE::coordinate_type &val) {
  M2_TYPEDEFS;
  mat3 m3;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      m3(i, j) = mi(i, j);
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(m3, Eigen::ComputeFullU);

  //const mat3 U = svd.matrixU();
  //const vec3 S = svd.singularValues();

  const Eigen::Matrix3d U = svd.matrixU();
  const Eigen::VectorXd S = svd.singularValues();
  val = coordinate_type(S[0], S[1], S[2]);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      mi(i, j) = U(i, j);
}

template <typename SPACE> class surface_calculator {
  M2_TYPEDEFS;

public:
  T baryArea(face_vertex_ptr fv) {
    // assumes triangels
    coordinate_type c0 = fv->coordinate();
    coordinate_type c1 = fv->next()->coordinate();
    coordinate_type c2n = fv->vnext()->next()->coordinate();
    coordinate_type c2p = fv->vprev()->next()->coordinate();

    coordinate_type c1h = 0.5 * (c0 + c1);
    coordinate_type c2nh = (c0 + c1 + c2n) / 3.0;
    coordinate_type c2ph = (c0 + c1 + c2p) / 3.0;

    coordinate_type dc10 = c1h - c0;
    coordinate_type dc20n = c2nh - c0;
    coordinate_type dc20p = c2ph - c0;
    T an = 0.5 * va::norm(va::cross(dc10, dc20n));
    T ap = 0.5 * va::norm(va::cross(dc10, dc20p));

    return an + ap;
  }

  T baryArea(vertex_ptr v) {
    face_vertex_ptr itb = v->fbegin();
    face_vertex_ptr ite = v->fend();
    bool at_head = false;
    T aTotal = 0;
    int i = 0;
    if (itb && ite) {
      while (!at_head && i < 40) {
        at_head = itb == ite;
        T aij = baryArea(itb);
        aTotal += aij;
        i++;
      }
    }
    return aTotal;
  }

  T getEdgeWeight(edge_ptr ei) {
    face_vertex_ptr fv1 = ei->v1();
    face_vertex_ptr fv2 = ei->v2();
    return fv1->cotan() + fv1->cotan();
  }

  T willmore(face_vertex_ptr fv) {

    coordinate_type ci = fv->coordinate();
    coordinate_type cj = fv->next()->coordinate();
    coordinate_type ck = fv->prev()->coordinate();
    coordinate_type cl = fv->vnext()->next()->coordinate();
    coordinate_type A = cj - ck;
    A.normalize();
    coordinate_type B = cl - cj;
    B.normalize();
    coordinate_type C = cl - ci;
    B.normalize();
    coordinate_type D = ci - ck;
    B.normalize();
    return va::dot(A, C) * va::dot(B, D) - va::dot(A, B) * va::dot(C, D) -
           va::dot(B, C) * va::dot(D, A);
  }

  template <typename TYPE>
  void calcDiffuseQuantity(m2::surf<SPACE> &in, vector<TYPE> &vertexWeights,
                           T amt) {

    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<edge_ptr> &tedges = in.get_edges();
    for (int i = 0; i < 4; i++) {
      vector<TYPE> tempWeights = vertexWeights;
      for (long i = 0; i < tverts.size(); i++) {
        if (tverts[i] && tverts[i]->size() > 0) {
          vertex_ptr v = tverts[i];
          face_vertex_ptr itb = v->fbegin();
          face_vertex_ptr ite = v->fend();
          bool at_head = false;
          T wTotal = 0;
          TYPE qTotal = 0;

          int i = 0;
          TYPE qi = tempWeights[v->position_in_set()];
          if (itb && ite) {
            while (!at_head && i < 40) {
              at_head = itb == ite;
              TYPE qj = tempWeights[itb->next()->vertex()->position_in_set()];
              // T wij = getEdgeWeight(itb->edge());;
              T wij = 1.0;
              TYPE qij = wij * (qj - qi);
              qTotal += qij;
              wTotal += wij;
              itb = itb->vnext();
              i++;
            }
            // vertexWeights[v->position_in_set()] = wTotal;
            if (wTotal < 1e-9)
              continue;
            vertexWeights[v->position_in_set()] += amt * qTotal / wTotal;
          }
        }
      }
    }
  }

  void calcDiffuseQuantity(m2::surf<SPACE> &in, vector<T> &vertexWeights,
                           T amt) {

    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<edge_ptr> &tedges = in.get_edges();
    for (int i = 0; i < 4; i++) {
      vector<T> tempWeights = vertexWeights;
      for (long i = 0; i < tverts.size(); i++) {
        if (tverts[i] && tverts[i]->size() > 0) {
          vertex_ptr v = tverts[i];
          face_vertex_ptr itb = v->fbegin();
          face_vertex_ptr ite = v->fend();
          bool at_head = false;
          T wTotal = 0;
          T qTotal = 0;

          int i = 0;
          T qi = tempWeights[v->position_in_set()];
          if (itb && ite) {
            while (!at_head && i < 40) {
              at_head = itb == ite;
              T qj = tempWeights[itb->next()->vertex()->position_in_set()];
              // T wij = getEdgeWeight(itb->edge());;
              T wij = 1.0;
              qTotal += wij * (qj - qi);
              wTotal += wij;
              itb = itb->vnext();
              i++;
            }
            // vertexWeights[v->position_in_set()] = wTotal;
            if (wTotal < 1e-9)
              continue;
            vertexWeights[v->position_in_set()] += amt * qTotal / wTotal;
          }
        }
      }
    }
  }

  void calcCurveFlowNormal(m2::surf<SPACE> &in, vector<T> &vertexWeights,
                           vector<T> &edgeWeights) {
    // TIMER function//TIMER(__FUNCTION__);

    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<edge_ptr> &tedges = in.get_edges();
    vertexWeights.resize(tverts.size());
    edgeWeights.resize(tedges.size());
    for (long i = 0; i < tedges.size(); i++) {
      if (!tedges[i])
        continue;
      T l = tedges[i]->length();
      T A1 = tedges[i]->v1()->face()->area();
      T A2 = tedges[i]->v2()->face()->area();
      if (l > 1e-12 && A1 > 1e-12 && A2 > 1e-12) {
        edgeWeights[tedges[i]->position_in_set()] = getEdgeWeight(tedges[i]);
      } else
        edgeWeights[tedges[i]->position_in_set()] = 0;
    }

    for (long i = 0; i < tverts.size(); i++) {
      if (tverts[i] && tverts[i]->size() > 0) {
        vertex_ptr v = tverts[i];
        face_vertex_ptr itb = v->fbegin();
        face_vertex_ptr ite = v->fend();
        bool at_head = false;
        coordinate_type kA(0, 0, 0);
        T wTotal = 0;
        T aTotal = 0;
        int i = 0;
        if (itb && ite) {
          while (!at_head && i < 40) {
            at_head = itb == ite;

            T wij = edgeWeights[itb->edge()->position_in_set()];
            T aij = baryArea(itb);
            T polyArea = itb->face()->area();
            coordinate_type c0 = itb->coordinate();
            coordinate_type c1 = itb->next()->coordinate();
            wTotal += wij;
            aTotal += aij;
            kA += wij * (c1 - c0);

            itb = itb->vnext();
            i++;
          }

          // vertexWeights[v->position_in_set()] = wTotal;
          coordinate_type kAi = kA / 2.0 * aTotal;
          vertexWeights[v->position_in_set()] = kA.norm();
        }
      }
    }
  }

  void calcWillmoreEnergy(m2::surf<SPACE> &in, vector<T> &vertexWeights) {
    // TIMER function//TIMER(__FUNCTION__);

    vector<vertex_ptr> &tverts = in.get_vertices();
    vector<edge_ptr> &tedges = in.get_edges();
    vertexWeights.resize(tverts.size());

    for (long i = 0; i < tverts.size(); i++) {
      if (tverts[i] && tverts[i]->size() > 0) {
        vertex_ptr v = tverts[i];
        face_vertex_ptr itb = v->fbegin();
        face_vertex_ptr ite = v->fend();
        bool at_head = false;
        coordinate_type kA(0, 0, 0);
        T wTotal = 0;
        int i = 0;
        if (itb && ite) {
          while (!at_head && i < 40) {
            at_head = itb == ite;

            T wij = willmore(itb);

            T polyArea = itb->face()->area();
            coordinate_type c0 = itb->coordinate();
            coordinate_type c1 = itb->next()->coordinate();
            wTotal += wij;

            itb = itb->vnext();
            i++;
          }

          // vertexWeights[v->position_in_set()] = wTotal;
          vertexWeights[v->position_in_set()] = wTotal;
        }
      }
    }
  }

  void calcDiffuseCurveFlowNormal(m2::surf<SPACE> &in,
                                  vector<T> &vertexWeights) {
    // TIMER function//TIMER(__FUNCTION__);
    // vector<T> vertexWeights;
    vector<T> edgeWeights;
    calcCurveFlowNormal(in, vertexWeights, edgeWeights);
    calcDiffuseQuantity(in, vertexWeights, 0.1);
  }

  std::vector<vertex_ptr> getLocalVertices(m2::surf<SPACE> &in,
                                           vertex_ptr seed,
                                           int maxRecDepth = 3) {
    // TIMER function//TIMER(__FUNCTION__);
    std::vector<vertex_ptr> out;
    typedef std::pair<vertex_ptr, int> recPair;
    std::stack<recPair> stack;
    coordinate_type c0 = seed->coordinate();

    stack.push(recPair(seed, 0));
    while (stack.size() > 0) {
      recPair rp = stack.top();
      stack.pop();
      vertex_ptr vi = rp.first;
      int recDepth = rp.second;
      face_vertex_ptr fvb = vi->fbegin();
      face_vertex_ptr fve = vi->fend();
      coordinate_type ci = vi->coordinate();
      bool iterating = true;
      vi->flag = 1;
      out.push_back(vi);
      while (iterating) {
        iterating = fvb != fve;
        vertex_ptr vj = fvb->next()->vertex();
        coordinate_type cj = vj->coordinate();
        if (vj->flag == 0 && recDepth < maxRecDepth) {
          stack.push(recPair(vj, recDepth + 1));
        }
        fvb = fvb->vnext();
      }
    }

    for (int i = 0; i < out.size(); i++) {
      out[i]->flag = 0;
    }
    return out;
  }

  void calcCovariance(m2::surf<SPACE> &in, vertex_ptr v,
                      coordinate_type &covVals, mat3 &covTens, T dx) {
    if (v->size() == 0)
      return;
    vector<vertex_ptr> lverts = getLocalVertices(in, v, 2);
    for (int i = 0; i < lverts.size(); i++) {

      coordinate_type c0 = v->coordinate();
      coordinate_type c1 = lverts[i]->coordinate();
      coordinate_type dc = c1 - c0;
      T dist = norm(dc);
      // T wij = 1.0;
      T exponent = -dist * dist / (dx * dx);
      T wij = exp(exponent);
      for (int l = 0; l < 3; l++)
        for (int m = 0; m < 3; m++) {
          covTens(l, m) += wij * dc[l] * dc[m];
        }
    }

    calcSVD<SPACE>(covTens, covVals);
    // std::cout << lverts.size() << " " << covVals << std::endl;
  }

  std::vector<edge_ptr> getEdgesNearPoint(m2::surf<SPACE> &in,
                                          vertex_ptr seed, T eps,
                                          int maxRecDepth = 3) {
    // TIMER function//TIMER(__FUNCTION__);
    std::vector<edge_ptr> out;
    std::vector<vertex_ptr> cleanup;
    typedef std::pair<vertex_ptr, int> recPair;
    std::stack<recPair> stack;
    coordinate_type c0 = seed->coordinate();

    stack.push(recPair(seed, 0));
    while (stack.size() > 0) {
      recPair rp = stack.top();
      stack.pop();
      vertex_ptr vi = rp.first;
      int recDepth = rp.second;
      face_vertex_ptr fvb = vi->fbegin();
      face_vertex_ptr fve = vi->fend();
      coordinate_type ci = vi->coordinate();
      bool iterating = true;

      while (iterating) {
        iterating = fvb != fve;
        vertex_ptr vj = fvb->next()->vertex();
        edge_ptr ej = fvb->edge();
        coordinate_type cj1 = ej->v1()->vertex()->coordinate();
        coordinate_type cj2 = ej->v2()->vertex()->coordinate();
        distance_calculator<SPACE> calc;
        T s;
        T d = calc.distanceFromLine(cj1, cj2, c0, s);

        if (vj->flag != 1 && recDepth < maxRecDepth) {
          vj->flag = 1;
          stack.push(recPair(vj, recDepth + 1));
          cleanup.push_back(vj);
        }

        if (ej->flag != 1 && d < eps) {
          if (s > 0.999 || s < 1e-4) {
            if (ej->v1()->vertex()->flag != 1) {
              out.push_back(ej);
              ej->v1()->vertex()->flag = 1;
            }
            if (ej->v2()->vertex()->flag != 1) {
              out.push_back(ej);
              ej->v2()->vertex()->flag = 1;
            }
          } else
            out.push_back(ej);

          ej->flag = 1;
        }
        fvb = fvb->vnext();
      }
    }

    for (int i = 0; i < out.size(); i++) {
      out[i]->flag = 0;
    }
    for (int i = 0; i < cleanup.size(); i++) {
      cleanup[i]->flag = 0;
    }
    return out;
  }

  std::vector<edge_ptr> getLocalEdges(m2::surf<SPACE> &in, vertex_ptr seed,
                                      T eps) {
    // TIMER function//TIMER(__FUNCTION__);
    std::vector<edge_ptr> out;
    int maxRecDepth = 2;
    typedef std::pair<vertex_ptr, int> recPair;
    std::stack<recPair> stack;
    coordinate_type c0 = seed->coordinate();
    stack.push(recPair(seed, 0));
    while (stack.size() > 0) {
      recPair rp = stack.top();
      stack.pop();
      vertex_ptr vi = rp.first;
      int recDepth = rp.second;
      vi->flag = 1;
      face_vertex_ptr fvb = vi->fbegin();
      face_vertex_ptr fve = vi->fend();
      coordinate_type ci = vi->coordinate();
      bool iterating = true;
      while (iterating) {
        iterating = fvb != fve;
        vertex_ptr vj = fvb->next()->vertex();
        coordinate_type cj = vj->coordinate();
        if (vj->flag == 0 && (cj - c0).norm() < eps && recDepth < maxRecDepth) {
          stack.push(recPair(vj, recDepth + 1));
        }
        if (fvb->edge()->flag == 0) {
          fvb->edge()->flag = 1;
          out.push_back(fvb->edge());
        }
        fvb = fvb->vnext();
      }
    }

    for (int i = 0; i < out.size(); i++) {
      out[i]->v1()->vertex()->flag = 0;
      out[i]->v2()->vertex()->flag = 0;
      out[i]->flag = 0;
    }
    return out;
  }

  void calculateBiDirection(m2::surf<SPACE> &in, vertex_ptr v,
                            coordinate_type &w, mat3 &cov, T dx) {

    // TIMER function//TIMER(__FUNCTION__);

    vector<edge_ptr> edges = getLocalEdges(in, v, 2.0 * dx);
    T cumArea = 0;

    for (int j = 0; j < edges.size(); j++) {
      edge_ptr ej = edges[j];
      cumArea += ej->v1()->face()->area();
      cumArea += ej->v2()->face()->area();
    }
    cumArea *= 0.5;
    for (int j = 0; j < edges.size(); j++) {
      edge_ptr ej = edges[j];
      coordinate_type c1 = ej->v1()->coordinate();
      coordinate_type c2 = ej->v2()->coordinate();
      coordinate_type n1 = ej->v1()->face()->normal();
      coordinate_type n2 = ej->v2()->face()->normal();
      T B = dot(n1, n2) / cumArea;
      coordinate_type dc = c1 - c2;
      for (int l = 0; l < 3; l++)
        for (int m = 0; m < 3; m++) {
          cov(l, m) += B * dc[l] * dc[m];
        }
    }
    calcSVD<SPACE>(cov, w);
  }

  void calculateBiDirectionField(surf_ref in, vector<mat3> &directionField,
                                 vector<coordinate_type> &directionWeights,
                                 T dx) {

    // TIMER function//TIMER(__FUNCTION__);
    vector<vertex_ptr> &tverts = in.get_vertices();
    for (long i = 0; i < tverts.size(); i++) {
      if (tverts[i] && tverts[i]->fbegin() && tverts[i]->fend()) {
        vertex_ptr v = tverts[i];
        mat3 cov;
#if 1
        vector<edge_ptr> edges = getLocalEdges(in, v, 2.0 * dx);
        T cumArea = 0;

        for (int j = 0; j < edges.size(); j++) {
          edge_ptr ej = edges[j];
          cumArea += ej->v1()->face()->area();
          cumArea += ej->v2()->face()->area();
        }
        cumArea *= 0.5;
        for (int j = 0; j < edges.size(); j++) {
          edge_ptr ej = edges[j];
          coordinate_type c1 = ej->v1()->coordinate();
          coordinate_type c2 = ej->v2()->coordinate();
          coordinate_type n1 = ej->v1()->face()->normal();
          coordinate_type n2 = ej->v2()->face()->normal();
          T B = dot(n1, n2) / cumArea;
          coordinate_type dc = c1 - c2;
          for (int l = 0; l < 3; l++)
            for (int m = 0; m < 3; m++) {
              directionField[v->position_in_set()](l, m) += B * dc[l] * dc[m];
            }
        }

#endif
      }
    }
    // calcDiffuseQuantity<mat3>(in,directionField,0.1);
    for (int i = 0; i < directionField.size(); i++) {
      calcSVD<SPACE>(directionField[i], directionWeights[i]);
    }
  }

  void shadeVertices(surf_ref mesh) {
    vector<vertex_ptr> &verts = mesh.get_vertices();
    vector<T> edgeWeights;
    vector<T> vertexWeights;
    // calcCurveFlowNormal(mesh, vertexWeights, edgeWeights);

    // vector<T> vertexWeights;
    calcDiffuseCurveFlowNormal(mesh, vertexWeights);
    T maxW = 0;
    T minW = 0;

    for (int i = 0; i < vertexWeights.size(); i++) {
      T wi = vertexWeights[i];
      maxW = wi > maxW ? wi : maxW;
      minW = wi < minW ? wi : minW;
    }

    for (int i = 0; i < verts.size(); i++) {
      T denom = maxW - minW;
      denom = denom < 40.00 ? denom : 40.0;
      // denom = denom > 10.0 ? 10.0:denom;
      T r = vertexWeights[i] / denom;
      verts[i]->color.r = 0.75 - 0.1 * r;
      verts[i]->color.g = 0.75 - 0.1 * r;
      verts[i]->color.b = 0.75 - 0.1 * r;
    }
  }

  void shadeVerticesWillmore(surf_ref mesh) {
    vector<vertex_ptr> &verts = mesh.get_vertices();
    vector<T> vertexWeights;
    calcWillmoreEnergy(mesh, vertexWeights);

    // vector<T> vertexWeights;
    // calcDiffuseCurveFlowNormal(mesh, vertexWeights);
    T maxW = 0;
    T minW = 0;

    for (int i = 0; i < vertexWeights.size(); i++) {
      T wi = vertexWeights[i];
      maxW = wi > maxW ? wi : maxW;
      minW = wi < minW ? wi : minW;
    }

    for (int i = 0; i < verts.size(); i++) {
      // minW = -2.0e-05;
      // maxW =  2.0e-05;
      T denom = maxW - minW;
      // denom = denom < 1000.00 ? denom:1000.0;
      // denom = denom > 10.0 ? 10.0:denom;
      T r = (vertexWeights[i] - minW) / denom;
      // std::cout << r << " " << vertexWeights[i] << " " << minW << " "
      // 	  << maxW << " " << denom << std::endl;
      verts[i]->color.r = 0.2;
      verts[i]->color.g = 0.75 - 0.05 * r;
      verts[i]->color.b = 1.0 - 0.1 * r;
    }
  }

  void shadeVerticesWinding(surf_ref mesh) {
    vector<vertex_ptr> &verts = mesh.get_vertices();

    T maxW = 0;
    T minW = 0;

    for (int i = 0; i < verts.size(); i++) {
      T wi = verts[i]->winding;
      maxW = wi > maxW ? wi : maxW;
      minW = wi < minW ? wi : minW;
    }

    for (int i = 0; i < verts.size(); i++) {
      T denom = maxW - minW;
      denom = denom < 100.00 ? denom : 100.0;
      // denom = denom > 10.0 ? 10.0:denom;
      T r = verts[i]->winding / denom;
      verts[i]->color.r = r;
    }
  }
};
} // namespace m2
#endif
