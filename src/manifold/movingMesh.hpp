#ifndef __M2MOVING__
#define __M2MOVING__

#include <stack> 

#ifdef _OPENMP
# include <omp.h>
#endif

#include "m2Includes.h"
#include "geometry_types.hpp"
#include "conj_grad.hpp"
#include "modify.hpp"
#include "remesh.hpp"
#include "octree.hpp"
#include "bins.hpp"
#include "debugger.h"

#include <cmath>

namespace m2 {

  template <typename SPACE>
  inline void calcSVD(typename SPACE::coordinate_type * vec,
		      typename SPACE::coordinate_type & val){
    M2_TYPEDEFS;
    Eigen::Matrix3f m3;
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
	m3(i,j) = vec[i][j];
    Eigen::JacobiSVD<Eigen::Matrix3f>  svd(m3,Eigen::ComputeFullU);

    const Eigen::Matrix3f U = svd.matrixU();
    const Eigen::VectorXf S = svd.singularValues();
    val = coordinate_type(S[0],S[1],S[2]);
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
	vec[j][i] = U(i,j);
  }

  template <typename SPACE>
  inline void calcSVD(typename SPACE::matrix3 & mi,
		      typename SPACE::coordinate_type & val){
    M2_TYPEDEFS;
    Eigen::Matrix3f m3;
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
	m3(i,j) = mi(i,j);
    Eigen::JacobiSVD<Eigen::Matrix3f>  svd(m3,Eigen::ComputeFullU);

    const Eigen::Matrix3f U = svd.matrixU();
    const Eigen::VectorXf S = svd.singularValues();
    val = coordinate_type(S[0],S[1],S[2]);

    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
	mi(i,j) = U(i,j);
  }

  template <typename SPACE>
  class mesh_calculator{ M2_TYPEDEFS;
  public:
    T baryArea(face_vertex_ptr fv){
      //assumes triangels
      coordinate_type c0 = fv->coordinate();
      coordinate_type c1 = fv->next()->coordinate();
      coordinate_type c2n = fv->vnext()->next()->coordinate();
      coordinate_type c2p = fv->vprev()->next()->coordinate();

      coordinate_type c1h = 0.5*(c0+c1);
      coordinate_type c2nh = (c0+c1+c2n)/3.0;
      coordinate_type c2ph = (c0+c1+c2p)/3.0;

      coordinate_type dc10 = c1h - c0;
      coordinate_type dc20n = c2nh - c0;
      coordinate_type dc20p = c2ph - c0;
      T an = 0.5*cross(dc10,dc20n).norm();
      T ap = 0.5*cross(dc10,dc20p).norm();
      

      return an + ap;
    }

    T baryArea(vertex_ptr v){      
      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;
      T aTotal = 0;
      int i = 0;
      if(itb && ite){
	while (!at_head && i < 40) {
	  at_head = itb==ite;	  
	  T aij = baryArea(itb);
	  aTotal += aij;
	  i++;
	}
	
      }
      return aTotal;
    }

    T cotan(coordinate_type c0,
	    coordinate_type c1,
	    coordinate_type c2){;

      coordinate_type dc10 = c1 - c0;
      coordinate_type dc20 = c2 - c0;
      //T denom = abs(dc10)*abs(dc20);
      T cosP = dot(dc10,dc20);
      T sinP = cross(dc10,dc20).norm();
      T cotP = cosP/sinP;
      if(sinP > 1e-12)
	return cotP;  
      else
	return 1.0;
    }

    T cotan(face_vertex_ptr fv){
      //assumes triangls
      coordinate_type c0 = fv->prev()->coordinate();
      coordinate_type c1 = fv->coordinate();
      coordinate_type c2 = fv->next()->coordinate();
      return cotan(c0,c1,c2);
    }


    T getEdgeWeight(edge_ptr ei){
      face_vertex_ptr fv1 = ei->v1();
      face_vertex_ptr fv2 = ei->v2();      
      return cotan(fv1) + cotan(fv2);      
    }


    T willmore(face_vertex_ptr fv){

      coordinate_type ci = fv->coordinate();
      coordinate_type cj = fv->next()->coordinate();
      coordinate_type ck = fv->prev()->coordinate();
      coordinate_type cl = fv->vnext()->next()->coordinate();
      coordinate_type A = cj - ck; A.normalize();
      coordinate_type B = cl - cj; B.normalize();
      coordinate_type C = cl - ci; B.normalize();
      coordinate_type D = ci - ck; B.normalize();
      return dot(A,C)*dot(B,D) - dot(A,B)*dot(C,D) - dot(B,C)*dot(D,A);
    }

    template<typename TYPE>
    void calcDiffuseQuantity(m2::control<SPACE>& in,
			     vector<TYPE> & vertexWeights,
			     T amt){

      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<edge_ptr>& tedges = in.get_edges();
      for(int i = 0; i < 4; i++){
	vector<TYPE> tempWeights = vertexWeights;
	for (long i = 0; i < tverts.size(); i++) {
	  if(tverts[i] && tverts[i]->size() > 0){
	    vertex_ptr v = tverts[i];
	    face_vertex_ptr itb = v->fbegin();
	    face_vertex_ptr ite = v->fend();
	    bool at_head = false;
	    T wTotal = 0;
	    TYPE qTotal = 0;

	    int i = 0;
	    TYPE qi = tempWeights[v->position_in_set()];
	    if(itb && ite){
	      while (!at_head && i < 40) {
		at_head = itb==ite;	      
		TYPE qj = tempWeights[itb->next()->vertex()->position_in_set()];
		//T wij = getEdgeWeight(itb->edge());;
		T wij = 1.0;
		TYPE qij = wij*(qj-qi);
		qTotal += qij;
		wTotal += wij;
		itb = itb->vnext();
		i++;
	      }	
	      //vertexWeights[v->position_in_set()] = wTotal;
	      if(wTotal < 1e-9) continue;
	      vertexWeights[v->position_in_set()] += amt*qTotal/wTotal;
	    }
	  }
	}
      }
    }

    void calcDiffuseQuantity(m2::control<SPACE>& in,
			     vector<T> & vertexWeights,
			     T amt){

      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<edge_ptr>& tedges = in.get_edges();
      for(int i = 0; i < 4; i++){
	vector<T> tempWeights = vertexWeights;
	for (long i = 0; i < tverts.size(); i++) {
	  if(tverts[i] && tverts[i]->size() > 0){
	    vertex_ptr v = tverts[i];
	    face_vertex_ptr itb = v->fbegin();
	    face_vertex_ptr ite = v->fend();
	    bool at_head = false;
	    T wTotal = 0;
	    T qTotal = 0;

	    int i = 0;
	    T qi = tempWeights[v->position_in_set()];
	    if(itb && ite){
	      while (!at_head && i < 40) {
		at_head = itb==ite;	      
		T qj = tempWeights[itb->next()->vertex()->position_in_set()];
		//T wij = getEdgeWeight(itb->edge());;
		T wij = 1.0;
		qTotal += wij*(qj-qi);
		wTotal += wij;
		itb = itb->vnext();
		i++;
	      }	
	      //vertexWeights[v->position_in_set()] = wTotal;
	      if(wTotal < 1e-9) continue;
	      vertexWeights[v->position_in_set()] += amt*qTotal/wTotal;
	    }
	  }
	}
      }
    }

    void calcCurveFlowNormal(m2::control<SPACE>& in,
			     vector<T> & vertexWeights,
			     vector<T> & edgeWeights){
      TIMER functionTimer(__FUNCTION__);

      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<edge_ptr>& tedges = in.get_edges();
      vertexWeights.resize(tverts.size());
      edgeWeights.resize(tedges.size());
      for (long i = 0; i < tedges.size(); i++) {
	if(!tedges[i]) continue;
	T l = tedges[i]->length();
	T A1 = tedges[i]->v1()->face()->area();
	T A2 = tedges[i]->v2()->face()->area();
	if(l > 1e-12 && A1 > 1e-12 && A2 > 1e-12 ){
	  edgeWeights[tedges[i]->position_in_set()] 
	    = getEdgeWeight(tedges[i]);
	}
	else edgeWeights[tedges[i]->position_in_set()] = 0;
      }

      for (long i = 0; i < tverts.size(); i++) {
	if(tverts[i] && tverts[i]->size() > 0){
	  vertex_ptr v = tverts[i];
	  face_vertex_ptr itb = v->fbegin();
	  face_vertex_ptr ite = v->fend();
	  bool at_head = false;
	  coordinate_type kA(0,0,0);
	  T wTotal = 0;
	  T aTotal = 0;
	  int i = 0;
	  if(itb && ite){
	    while (!at_head && i < 40) {
	      at_head = itb==ite;

	      T wij = edgeWeights[itb->edge()->position_in_set()];
	      T aij = baryArea(itb);
	      T polyArea = itb->face()->area();
	      coordinate_type c0 = itb->coordinate();
	      coordinate_type c1 = itb->next()->coordinate();
	      wTotal += wij;
	      aTotal += aij;
	      kA += wij*(c1-c0);
	      
	      itb = itb->vnext();
	      i++;
	    }
	
	    //vertexWeights[v->position_in_set()] = wTotal;
	    vertexWeights[v->position_in_set()] = norm(kA/(2.0*aTotal));
	  }
	}
      }
    }

    void calcWillmoreEnergy(m2::control<SPACE>& in,
			    vector<T> & vertexWeights){
      TIMER functionTimer(__FUNCTION__);

      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<edge_ptr>& tedges = in.get_edges();
      vertexWeights.resize(tverts.size());

      for (long i = 0; i < tverts.size(); i++) {
	if(tverts[i] && tverts[i]->size() > 0){
	  vertex_ptr v = tverts[i];
	  face_vertex_ptr itb = v->fbegin();
	  face_vertex_ptr ite = v->fend();
	  bool at_head = false;
	  coordinate_type kA(0,0,0);
	  T wTotal = 0;
	  int i = 0;
	  if(itb && ite){
	    while (!at_head && i < 40) {
	      at_head = itb==ite;

	      T wij = willmore(itb);

	      T polyArea = itb->face()->area();
	      coordinate_type c0 = itb->coordinate();
	      coordinate_type c1 = itb->next()->coordinate();
	      wTotal += wij;
	      
	      itb = itb->vnext();
	      i++;
	    }
	
	    //vertexWeights[v->position_in_set()] = wTotal;
	    vertexWeights[v->position_in_set()] = wTotal;

	  }
	}
      }
    }

    void calcDiffuseCurveFlowNormal(m2::control<SPACE>& in,
				    vector<T> & vertexWeights){
      TIMER functionTimer(__FUNCTION__);
      //vector<T> vertexWeights;
      vector<T> edgeWeights;
      calcCurveFlowNormal(in, vertexWeights, edgeWeights);
      calcDiffuseQuantity(in,vertexWeights,0.1);
    }

    std::vector<vertex_ptr> getLocalVertices(m2::control<SPACE>& in,
					     vertex_ptr seed,
					     int maxRecDepth = 3){
      TIMER functionTimer(__FUNCTION__);   
      std::vector<vertex_ptr> out;
      typedef std::pair<vertex_ptr,int> recPair;
      std::stack<recPair> stack;
      coordinate_type c0 = seed->coordinate();

      stack.push(recPair(seed,0));
      while(stack.size() > 0){
	recPair rp = stack.top(); stack.pop();
	vertex_ptr vi = rp.first;
	int recDepth = rp.second;
	face_vertex_ptr fvb = vi->fbegin();
	face_vertex_ptr fve = vi->fend();
	coordinate_type ci = vi->coordinate();
	bool iterating = true;
	vi->flag = 1;
	out.push_back(vi);
	while(iterating){
	  iterating = fvb != fve;
	  vertex_ptr vj = fvb->next()->vertex();
	  coordinate_type cj = vj->coordinate();
	  if(vj->flag == 0 &&
	     recDepth <  maxRecDepth){	    
	    stack.push(recPair(vj,recDepth + 1));
	  }
	  fvb = fvb->vnext();
	}
      }

      for(int i = 0; i < out.size(); i++){
	out[i]->flag = 0;
      }
      return out;
    }

    void calcCovariance(m2::control<SPACE>& in, 
			vertex_ptr v, 
			coordinate_type& covVals,
			mat3 & covTens, T dx){
      if(v->size() ==  0) return;
      vector<vertex_ptr> lverts = getLocalVertices(in, v, 2);
      for(int i = 0; i < lverts.size(); i++){
	
	coordinate_type c0 = v->coordinate();
	coordinate_type c1 = lverts[i]->coordinate();
	coordinate_type dc = c1-c0;
	T dist = norm(dc);
	//T wij = 1.0;	
	T exponent = -dist*dist/(dx*dx);
	T  wij = exp(exponent);
	for(int l = 0; l < 3; l++)
	  for(int m = 0; m < 3; m++){
	    covTens(l,m) +=wij*dc[l]*dc[m];
	  }
      }

      calcSVD<SPACE>(covTens, covVals);
      //std::cout << lverts.size() << " " << covVals << std::endl;
    }

    std::vector<edge_ptr> getEdgesNearPoint(m2::control<SPACE>& in,
					    vertex_ptr seed,
					    T eps, int maxRecDepth = 3){
      TIMER functionTimer(__FUNCTION__);   
      std::vector<edge_ptr> out;
      std::vector<vertex_ptr> cleanup;
      typedef std::pair<vertex_ptr,int> recPair;
      std::stack<recPair> stack;
      coordinate_type c0 = seed->coordinate();

      stack.push(recPair(seed,0));
      while(stack.size() > 0){
	recPair rp = stack.top(); stack.pop();
	vertex_ptr vi = rp.first;
	int recDepth = rp.second;
	face_vertex_ptr fvb = vi->fbegin();
	face_vertex_ptr fve = vi->fend();
	coordinate_type ci = vi->coordinate();
	bool iterating = true;

	while(iterating){
	  iterating = fvb != fve;
	  vertex_ptr vj = fvb->next()->vertex();
	  edge_ptr ej = fvb->edge();
	  coordinate_type cj1 = ej->v1()->vertex()->coordinate();
	  coordinate_type cj2 = ej->v2()->vertex()->coordinate();
	  distance_calculator<SPACE> calc;
	  T s;
	  T d = calc.distanceFromLine(cj1,cj2,c0,s);

	  if(vj->flag != 1 &&
	     recDepth <  maxRecDepth){	    
	    vj->flag = 1;
	    stack.push(recPair(vj,recDepth + 1));
	    cleanup.push_back(vj);
	  }

	  if(ej->flag != 1 &&
	     d < eps){
	    if(s > 0.999 || s < 1e-4){
	      if(ej->v1()->vertex()->flag != 1){
		out.push_back(ej);	
		ej->v1()->vertex()->flag = 1;
	      }
	      if(ej->v2()->vertex()->flag != 1){
		out.push_back(ej);
		ej->v2()->vertex()->flag = 1;	
	      }
	    }
	    else
	      out.push_back(ej);

	    ej->flag = 1;
	  }
	  fvb = fvb->vnext();
	}
      }

      for(int i = 0; i < out.size(); i++){
	out[i]->flag = 0;
      }
      for(int i = 0; i < cleanup.size(); i++){
	cleanup[i]->flag = 0;
      }
      return out;
    }
    
    std::vector<edge_ptr> getLocalEdges(m2::control<SPACE>& in,
					vertex_ptr seed,
					T eps){
      TIMER functionTimer(__FUNCTION__);   
      std::vector<edge_ptr> out;
      int maxRecDepth = 2;
      typedef std::pair<vertex_ptr,int> recPair;
      std::stack<recPair> stack;
      coordinate_type c0 = seed->coordinate();
      stack.push(recPair(seed,0));
      while(stack.size() > 0){
	recPair rp = stack.top(); stack.pop();
	vertex_ptr vi = rp.first;
	int recDepth = rp.second;
	vi->flag = 1;
	face_vertex_ptr fvb = vi->fbegin();
	face_vertex_ptr fve = vi->fend();
	coordinate_type ci = vi->coordinate();
	bool iterating = true;
	while(iterating){
	  iterating = fvb != fve;
	  vertex_ptr vj = fvb->next()->vertex();
	  coordinate_type cj = vj->coordinate();
	  if(vj->flag == 0 &&
	     (cj - c0).mag() < eps &&
	     recDepth <  maxRecDepth){	    
	    stack.push(recPair(vj,recDepth + 1));
	  }
	  if(fvb->edge()->flag == 0){
	    fvb->edge()->flag = 1;
	    out.push_back(fvb->edge());	    
	  }
	  fvb = fvb->vnext();
	}
      }

      for(int i = 0; i < out.size(); i++){
	out[i]->v1()->vertex()->flag = 0;
	out[i]->v2()->vertex()->flag = 0;
	out[i]->flag = 0;
      }
      return out;
    }

    void calculateBiDirection(m2::control<SPACE>& in, 
			      vertex_ptr v, 
			      coordinate_type & w,
			      mat3 & cov,
			      T dx){

      TIMER functionTimer(__FUNCTION__);

      vector<edge_ptr> edges = getLocalEdges(in, v, 2.0*dx);
      T cumArea = 0;

      for(int j = 0; j < edges.size(); j++){
	edge_ptr ej = edges[j];
	cumArea += ej->v1()->face()->area();
	cumArea += ej->v2()->face()->area();
      }
      cumArea *= 0.5;
      for(int j = 0; j < edges.size(); j++){
	edge_ptr ej = edges[j];
	coordinate_type c1 = ej->v1()->coordinate();
	coordinate_type c2 = ej->v2()->coordinate();
	coordinate_type n1 = ej->v1()->face()->normal();
	coordinate_type n2 = ej->v2()->face()->normal();
	T B = dot(n1,n2)/cumArea;
	coordinate_type dc = c1-c2;
	for(int l = 0; l < 3; l++)
	  for(int m = 0; m < 3; m++){
	    cov(l,m) += B*dc[l]*dc[m];
	  }
      }
      calcSVD<SPACE>(cov, w);
    }
    
    void calculateBiDirectionField(control_ref in,
				   vector<mat3> & directionField,
				   vector<coordinate_type>  & directionWeights,
				   T dx){

      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr>& tverts = in.get_vertices();
      for (long i = 0; i < tverts.size(); i++) {
	if(tverts[i] && 
	   tverts[i]->fbegin() &&
	   tverts[i]->fend()){
	  vertex_ptr v = tverts[i];
	  mat3 cov;
#if 1
	  vector<edge_ptr> edges = getLocalEdges(in, v, 2.0*dx);
	  T cumArea = 0;

	  for(int j = 0; j < edges.size(); j++){
	    edge_ptr ej = edges[j];
	    cumArea += ej->v1()->face()->area();
	    cumArea += ej->v2()->face()->area();
	  }
	  cumArea *= 0.5;
	  for(int j = 0; j < edges.size(); j++){
	    edge_ptr ej = edges[j];
	    coordinate_type c1 = ej->v1()->coordinate();
	    coordinate_type c2 = ej->v2()->coordinate();
	    coordinate_type n1 = ej->v1()->face()->normal();
	    coordinate_type n2 = ej->v2()->face()->normal();
	    T B = dot(n1,n2)/cumArea;
	    coordinate_type dc = c1-c2;
	    for(int l = 0; l < 3; l++)
	      for(int m = 0; m < 3; m++){
		directionField[v->position_in_set()](l,m) += B*dc[l]*dc[m];
	      }
	  }

#endif

	}
      }
      //calcDiffuseQuantity<mat3>(in,directionField,0.1);
      for(int i = 0; i < directionField.size(); i++){
	calcSVD<SPACE>(directionField[i], directionWeights[i]);
      }

    }

    void shadeVertices(control_ref mesh){
      vector<vertex_ptr>& verts = mesh.get_vertices();
      vector<T> edgeWeights;
      vector<T> vertexWeights;
      //calcCurveFlowNormal(mesh, vertexWeights, edgeWeights);

      //vector<T> vertexWeights;
      calcDiffuseCurveFlowNormal(mesh, vertexWeights);
      T maxW = 0;
      T minW = 0;

      for(int i = 0; i < vertexWeights.size(); i++){
	T wi = vertexWeights[i];
	maxW = wi > maxW ? wi:maxW;
	minW = wi < minW ? wi:minW;
      }

      for (int i = 0; i < verts.size(); i++){
	T denom = maxW-minW;
	denom = denom < 40.00 ? denom:40.0;
	//denom = denom > 10.0 ? 10.0:denom;
	T r = vertexWeights[i]/denom;
	verts[i]->color.r = 0.75-0.1*r;
	verts[i]->color.g = 0.75-0.1*r;
	verts[i]->color.b = 0.75-0.1*r;
      }
    }

    void shadeVerticesWillmore(control_ref mesh){
      vector<vertex_ptr>& verts = mesh.get_vertices();
      vector<T> vertexWeights;
      calcWillmoreEnergy(mesh, vertexWeights);

      //vector<T> vertexWeights;
      //calcDiffuseCurveFlowNormal(mesh, vertexWeights);
      T maxW = 0;
      T minW = 0;

      for(int i = 0; i < vertexWeights.size(); i++){
	T wi = vertexWeights[i];
	maxW = wi > maxW ? wi:maxW;
	minW = wi < minW ? wi:minW;
      }

      for (int i = 0; i < verts.size(); i++){
	// minW = -2.0e-05;
	// maxW =  2.0e-05;
	T denom = maxW-minW;
	// denom = denom < 1000.00 ? denom:1000.0;
	// denom = denom > 10.0 ? 10.0:denom;
	T r = (vertexWeights[i] - minW)/denom;
	// std::cout << r << " " << vertexWeights[i] << " " << minW << " " 
	// 	  << maxW << " " << denom << std::endl;
	verts[i]->color.r = 0.2; 
	verts[i]->color.g = 0.75 - 0.05*r;
	verts[i]->color.b = 1.0-0.1*r;
      }
    }

    void shadeVerticesWinding(control_ref mesh){
      vector<vertex_ptr>& verts = mesh.get_vertices();

      T maxW = 0;
      T minW = 0;

      for(int i = 0; i < verts.size(); i++){
	T wi = verts[i]->winding;
	maxW = wi > maxW ? wi:maxW;
	minW = wi < minW ? wi:minW;
      }

      for (int i = 0; i < verts.size(); i++){
	T denom = maxW-minW;
	denom = denom < 100.00 ? denom:100.0;
	//denom = denom > 10.0 ? 10.0:denom;
	T r = verts[i]->winding/denom;
	verts[i]->color.r = r;
      }
    }
  };

  template <typename SPACE>
  class mesh_filter{ M2_TYPEDEFS;
  public:


    void filterCutoff(m2::control<SPACE>& in,
		      T cutoff,
		      vector<T> & vertexWeights,
		      T strength){
      TIMER functionTimer(__FUNCTION__);   
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      m2::mesh_calculator<SPACE> calc;

      for (long i = 0; i < tverts.size(); i++) {

	if(!in.has_vertex(i)) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;

	face_vertex_ptr itb = v->fbegin();
	face_vertex_ptr ite = v->fend();
	bool at_head = false;
	coordinate_type kA(0,0,0);
	T wTotal = 0, aTotal = 0;
	coordinate_type c0 = tverts[i]->coordinate();
	int k = 0;
	while (!at_head) {
	  at_head = itb==ite;
	  //mesh_calculator<SPACE> calc;
	  T wij = calc.getEdgeWeight(itb->edge());
	  //T wij = 1.0;
	  coordinate_type c1 = itb->next()->coordinate();	  
	  T aij = calc.baryArea(itb);
	  
	  aTotal += aij;
	  wTotal += wij;

	  kA += wij*(c1-c0);
	  itb = itb->vnext();
	  k++;
	}
	T kT = vertexWeights[v->position_in_set()];

	if(wTotal > 1e-6 &&  kT < cutoff)
	  filteredCoordinates[i] = kA/wTotal;
	else
	  filteredCoordinates[i] = coordinate_type(0,0,0);
      }

      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	coordinate_type ci =  tverts[i]->coordinate();
	tverts[i]->coordinate() = ci + strength*filteredCoordinates[i];
      }
    }

    void filter(m2::control<SPACE>& in, T strength){
      TIMER functionTimer(__FUNCTION__);   
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      for (long i = 0; i < tverts.size(); i++) {

	if(!in.has_vertex(i)) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;

	face_vertex_ptr itb = v->fbegin();
	face_vertex_ptr ite = v->fend();
	bool at_head = false;
	coordinate_type kA(0,0,0);
	T wTotal = 0, aTotal = 0;
	coordinate_type c0 = tverts[i]->coordinate();
	int k = 0;
	while (!at_head) {
	  at_head = itb==ite;
	  mesh_calculator<SPACE> calc;
	  T wij = calc.getEdgeWeight(itb->edge());
	  //T wij = 1.0;
	  coordinate_type c1 = itb->next()->coordinate();	  
	  T aij = calc.baryArea(itb);
	  
	  aTotal += aij;
	  wTotal += wij;

	  kA += wij*(c1-c0);
	  itb = itb->vnext();
	  k++;
	}
	T kT = norm(kA/(2.0*aTotal));

	if(wTotal > 1e-10)
	  filteredCoordinates[i] = kA/wTotal;
	else
	  filteredCoordinates[i] = coordinate_type(0,0,0);
      }

      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	coordinate_type ci =  tverts[i]->coordinate();
	tverts[i]->coordinate() = ci + strength*filteredCoordinates[i];
      }
    }

    void cuspFilter(m2::control<SPACE>& in, T strength){
      TIMER functionTimer(__FUNCTION__);   
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      for (long i = 0; i < tverts.size(); i++) {

	if(!in.has_vertex(i)) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;

	face_vertex_ptr itb = v->fbegin();
	face_vertex_ptr ite = v->fend();
	bool at_head = false;
	coordinate_type kA(0,0,0);
	T wTotal = 0, aTotal = 0;
	coordinate_type c0 = tverts[i]->coordinate();
	int k = 0;
	while (!at_head) {
	  at_head = itb==ite;
	  mesh_calculator<SPACE> calc;
	  T wij = calc.getEdgeWeight(itb->edge());
	  edge_ptr e = itb->edge();
	  coordinate_type n1 = e->v1()->face()->normal();
	  coordinate_type n2 = e->v2()->face()->normal();
	  wij = 1.0 + dot(n1,n2);
	  //wij = dot(n1,n2);
	  coordinate_type c1 = itb->next()->coordinate();	  
	  
	  wTotal += wij;

	  kA += wij*(c1-c0);
	  itb = itb->vnext();
	  k++;
	}
	T kT = norm(kA/(2.0*aTotal));

	if(wTotal > 1e-10)
	  filteredCoordinates[i] = kA/wTotal;
	else
	  filteredCoordinates[i] = coordinate_type(0,0,0);
      }

      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	coordinate_type ci =  tverts[i]->coordinate();
	tverts[i]->coordinate() = ci + strength*filteredCoordinates[i];
      }
    }

    coordinate_type laplacianFilterVertex(vertex_ptr v){
      if(v->size() ==  0) return coordinate_type(0,0,0);
      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;
      coordinate_type kA(0,0,0);
      T wTotal = 0;
      coordinate_type c0 = v->coordinate();
      int k = 0;
      while (!at_head) {
	at_head = itb==ite;
	T wij = 1.0;
	coordinate_type c1 = itb->next()->coordinate();
	if(itb->next()->vertex()->flag == 0){
	  itb->next()->vertex()->flag += 1;
	  wTotal +=  wij;
	  kA += wij*(c1-c0);
	  k++;
	}
	itb = itb->vnext();
      }
      
      itb = v->fbegin();
      ite = v->fend();
      at_head = false;

      while (!at_head) {
	at_head = itb==ite;
	itb->next()->vertex()->flag = 0;
	itb = itb->vnext();	
      }

      if(wTotal > 1e-10){
	v->update_normal();
	coordinate_type w;
	coordinate_type dv = kA/wTotal;
	return dv;
      }
      else return coordinate_type(0,0,0);
    }
    
    coordinate_type projectOntoLine(const coordinate_type& v0, 
				    const coordinate_type& v1, 
				    const coordinate_type& pt){
      coordinate_type s = v1 - v0;
      coordinate_type v = pt - v0;
      coordinate_type vp = dot(v,s)/dot(s,s)*s;
      coordinate_type ptl = v0 + vp;
      return ptl;
    }
    
    void mlsFilter(m2::control<SPACE>& in, T strength, T dx){
      TIMER functionTimer(__FUNCTION__);   
      vector<T> vertexWeights;
      vector<T> edgeWeights;
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;

	vertex_ptr v = tverts[i];
	coordinate_type c0 = tverts[i]->coordinate();
	if(v->size() ==  0) continue;
	coordinate_type ca(0,0,0); T wTot = 0;
	m2::mesh_calculator<SPACE> calc;
	vector<vertex_ptr> lverts = calc.getLocalVertices(in, v, 2);

	for(int j = 0; j < lverts.size(); j++){
	  coordinate_type c1 = lverts[j]->coordinate();
	  coordinate_type dc = c1-c0;
	  T dist = norm(dc);
	  //T wij = 1.0;	
	  T exponent = -dist*dist/(dx*dx);
	  T  wij = exp(exponent);
	  wTot += wij;
	  ca += wij*dc;
	}
	ca /= wTot;
	ca += c0;

	mat3 m; coordinate_type w;
	calc.calcCovariance(in, v,w,m, dx);
	coordinate_type d(m(0,0), m(1,0), m(2,0));
	coordinate_type d2(m(0,2), m(1,2), m(2,2));
	

	coordinate_type cp = projectOntoLine(ca - 0.1*d, ca + 0.1*d, c0);
	filteredCoordinates[i] = cp - c0;
      }

      for (long i = 0; i < filteredCoordinates.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	tverts[i]->coordinate() += strength*filteredCoordinates[i];
      }
    }

    void laplacianFilter(m2::control<SPACE>& in, T strength){
      TIMER functionTimer(__FUNCTION__);   
      vector<T> vertexWeights;
      vector<T> edgeWeights;
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;

	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;
	filteredCoordinates[i] = laplacianFilterVertex(v);
      }

      for (long i = 0; i < filteredCoordinates.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	tverts[i]->coordinate() += strength*filteredCoordinates[i];
      }
    }

    void taubinFilter(m2::control<SPACE>& in, T a,T b){
      TIMER functionTimer(__FUNCTION__);   
      vector<T> vertexWeights;
      vector<T> edgeWeights;
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;

	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;
	filteredCoordinates[i] = laplacianFilterVertex(v);
      }

      for (long i = 0; i < filteredCoordinates.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	tverts[i]->coordinate() += a*filteredCoordinates[i];
      }

      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;
	filteredCoordinates[i] = laplacianFilterVertex(v);
      }

      for (long i = 0; i < filteredCoordinates.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	tverts[i]->coordinate() -= b*filteredCoordinates[i];
      }
    }

    void cacheTensor(m2::control<SPACE>& in, 
		     vector<coordinate_type*> & tensorArray,
		     vector<coordinate_type>  & singularArray){
      vector<vertex_ptr>& tverts = in.get_vertices();
      for (long i = 0; i < tverts.size(); i++) {
	
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;
	face_vertex_ptr itb = v->fbegin();
	face_vertex_ptr ite = v->fend();
	bool at_head = false;
	coordinate_type * cov = new coordinate_type[3];
	while (!at_head) {
	  at_head = itb==ite;
	  T wij = itb->face()->area();
	  coordinate_type N = itb->face()->normal();
	  for(int l = 0; l < 3; l++)
	    for(int m = 0; m < 3; m++){
	      cov[l][m] +=wij*N[l]*N[m];
	    }
	  itb = itb->vnext();
	
	}
	coordinate_type w;
	calcSVD<SPACE>(cov, w);
	tensorArray[i] = cov;
	singularArray[i] = w;
      }
    }

    coordinate_type nullLaplacianFilterVertex(vertex_ptr v,
					      coordinate_type cov[3], 
					      coordinate_type w){
      if(v->size() ==  0) return coordinate_type(0,0,0);
      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;
      coordinate_type kA(0,0,0);
      T wTotal = 0;
      coordinate_type c0 = v->coordinate();
      int k = 0;
      //coordinate_type cov[3];
      while (!at_head) {
	at_head = itb==ite;
	T wij = itb->face()->area();
	coordinate_type c1 = itb->next()->coordinate();
	wTotal +=  wij;
	kA += wij*(c1-c0);
	itb = itb->vnext();
	k++;
      }
      coordinate_type dv = kA/wTotal;
      if(wTotal > 1e-10){
	v->update_normal();
	
	coordinate_type u = cov[1]; u.normalize();
	coordinate_type v = cov[2]; v.normalize();
	coordinate_type outu = dot(dv,u)*u;
	coordinate_type outv = dot(dv,v)*v;
	
	T wp = w[1] + w[2];
	if(w[1] < 1e-3&& w[2] < 1e-3)
	  return outu + outv;
	else if(w[2] < 1e-3)
	  return outv;
	else return coordinate_type(0,0,0); 
      }
      else return coordinate_type(0,0,0);
    }

    void nullLaplacianFilter(m2::control<SPACE>& in, T strength){
      TIMER functionTimer(__FUNCTION__);   
      vector<T> vertexWeights;
      vector<T> edgeWeights;
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      vector<coordinate_type*> tensorArray;
      tensorArray.resize(tverts.size());
      vector<coordinate_type> singularArray;
      singularArray.resize(tverts.size());
      cacheTensor(in, tensorArray, singularArray);
      for(int k = 0; k<2; k++){
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (long i = 0; i < tverts.size(); i++) {
	  if(!in.has_vertex(i)) continue;
	  if(tverts[i]->pinned) continue;
	  vertex_ptr v = tverts[i];
	  if(v->size() ==  0) continue;
	  filteredCoordinates[i] = 
	    nullLaplacianFilterVertex(v, tensorArray[i], singularArray[i]);
	}

	for (long i = 0; i < filteredCoordinates.size(); i++) {
	  if(!in.has_vertex(i)) continue;
	  if(tverts[i]->pinned) continue;
	  tverts[i]->coordinate() += strength*filteredCoordinates[i];
	}
      }
      for(int i = 0; i < tensorArray.size(); i++){
	delete tensorArray[i];
      }
    }

    void flaggedFilter(m2::control<SPACE>& in, T strength){
      TIMER functionTimer(__FUNCTION__);   
      vector<T> vertexWeights;
      vector<T> edgeWeights;
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> filteredCoordinates;
      filteredCoordinates.resize(tverts.size());
      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;
	if(v->flag != 1)  continue;
	//filteredCoordinates[i] = laplacianFilterVertex(v);
	filteredCoordinates[i] = nullLaplacianFilterVertex(v);
      }

      for (long i = 0; i < filteredCoordinates.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	if(tverts[i]->flag != 1)  continue;
	tverts[i]->coordinate() += strength*filteredCoordinates[i];
      }
    }
  };

  template <typename SPACE>
  class set_operations{ M2_TYPEDEFS;
  public:
    bool flip_edges(control_ptr obj){
      TIMER functionTimer(__FUNCTION__);
      edge_array& edges = obj->get_edges();
      m2::construct<SPACE> cons;
      bool flipped = false;
      edge_array permEdges;
      for(int i = 0; i < edges.size(); i++){
	if(!obj->has_edge(i)) continue;
	edge_ptr e = edges[i];
	if(e->v1()->face()->size() == 3 && e->v2()->face()->size() == 3)
	  permEdges.push_back(edges[i]);
      }

      for(int i = 0; i < permEdges.size(); i++){
	int card = rand()%permEdges.size();
	edge_ptr et = permEdges[i];
	permEdges[i] = permEdges[card];
	permEdges[card] = et;
      }

      for(int i = 0; i < permEdges.size(); i++){
	edge_ptr e = permEdges[i];
	bool pinned =
	  e->v1()->vertex()->pinned == true &&
	  e->v2()->vertex()->pinned == true;
	bool notPinned =
	  e->v1()->vertex()->pinned != true &&
	  e->v2()->vertex()->pinned != true;
	if(pinned) continue;	
	if(notPinned) continue;


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
	T m01 = 1.0/(c0-c1).mag();
	T m12 = 1.0/(c1-c2).mag();
	T m23 = 1.0/(c2-c3).mag();
	T m30 = 1.0/(c3-c0).mag();

	T cos0 = dot((c1-c0),(c3-c0))*m01*m30;
	T cos1 = dot((c0-c1),(c2-c1))*m01*m12;
	T cos2 = dot((c1-c2),(c3-c2))*m12*m23;
	T cos3 = dot((c0-c3),(c2-c3))*m30*m23;
	//half angle cos^2(2a) = 0.5*(1+cos(a))
	T cSame =  acos(cos1) + acos(cos3); //corresponds to flipped edge
	T cFlip =  acos(cos0) + acos(cos2); //corresponds to flipped edge
	//T cSame =  (M_PI - (acos(cos1) + acos(cos3))); //corresponds to flipped edge

	T eFlip = cFlip*cFlip;
	T eSame = cSame*cSame;
#if 0
	T sin0 = cross(c1-c0,c3-c0).mag();
	T sin1 = cross(c0-c1,c2-c1).mag();
	T sin2 = cross(c1-c2,c3-c2).mag();
	T sin3 = cross(c0-c3,c2-c3).mag();
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

	//if(cSame > M_PI){
	if(cFlip < cSame){
	  m2::construct<SPACE> cons;
	  cons.flip_edge(obj, e);
	  // e->v1()->data *= -1;
	  // e->v2()->data *= -1;
	  flipped = true;
	}

      }
      return flipped;
    }
  };

  template <typename SPACE>
  class join_mesh{ M2_TYPEDEFS;
  public:
    
    //---------------------------------------------------------------------------
    //local structures
    //---------------------------------------------------------------------------    

    struct contact_manifold{
    public:
      
      int p1,p2;
      int p1type, p2type; //means that we can querie the collision type
      control_ptr mMesh;
      T length(){
      	if(p1type == p2type){
	  vector<vertex_ptr> & verts = this->mMesh->get_vertices();
	  T d0 = norm2(verts[p1]->coordinate() - verts[p2]->coordinate());
	  return d0;
	}
      }
    };

    struct contact_manifold_sort {
    public:
      bool operator() (contact_manifold c0, contact_manifold c1) { 
	//point 0
	//edge  1
	//tri   2
	T d0, d1;
	if(c0.p1type == c0.p2type){
	  vector<vertex_ptr> & verts = this->mMesh->get_vertices();
	  d0 = norm2(verts[c0.p1]->coordinate() - verts[c0.p2]->coordinate());
	}
	if(c1.p1type == c1.p2type){
	  vector<vertex_ptr> & verts = this->mMesh->get_vertices();
	  d1 = norm2(verts[c1.p1]->coordinate() - verts[c1.p2]->coordinate());
	}
	
	return (d0 < d1);
      }

      control_ptr mMesh;
    } mContactSorter;
    
    //in order to run unit test
    join_mesh(){
    };

    join_mesh(control_ptr mesh, T tolerance){
      tol = tolerance;
      mMesh = mesh;
    };

    //---------------------------------------------------------------------------
    //constructor and init
    //---------------------------------------------------------------------------   

    void init(){
      //initializes the trees necessary to make this all work
      vector<face_ptr> & faces = mMesh->get_faces();
      for(int i = 0; i < faces.size(); i++){
	triangle_type tempFace =  makeTriangle(faces[i]->fbegin(),
					       faces[i]->fbegin()->next(),
					       faces[i]->fbegin()->prev());
	oldFaces.push_back(tempFace);
      }      
      exteriorTree.build(oldFaces,10);
    }

    inline void dimPolynomial(coordinate_type a, coordinate_type va, 
			      coordinate_type b, coordinate_type vb, 
			      coordinate_type c, coordinate_type vc,
			      int i, int j, int k, 
			      T & o, T & p, T & q, T & r){
      o = va[i]*vb[j]*vc[k];
      p = a[i]*vb[j]*vc[k] + va[i]*b[j]*vc[k] + va[i]*vb[j]*c[k];
      q = va[i]*b[j]*c[k] + a[i]*vb[j]*c[k] + a[i]*b[j]*vc[k];
      r = a[i]*b[j]*c[k];
      //std::cout << o << " " << p << " " << q << " " << r << std::endl;
    }

    void buildPolynomial(coordinate_type a, coordinate_type va, 
			 coordinate_type b, coordinate_type vb, 
			 coordinate_type c, coordinate_type vc, 
			 T & p, T & q, T & r){
      T o0, o1, o2, o3, o4, o5;
      T p0, p1, p2, p3, p4, p5;
      T q0, q1, q2, q3, q4, q5;
      T r0, r1, r2, r3, r4, r5;
      dimPolynomial(a, va, b, vb, c, vc, 1,2,0, o0, p0, q0, r0);
      dimPolynomial(a, va, b, vb, c, vc, 2,1,0, o1, p1, q1, r1);
      dimPolynomial(a, va, b, vb, c, vc, 2,0,1, o2, p2, q2, r2);
      dimPolynomial(a, va, b, vb, c, vc, 0,2,1, o3, p3, q3, r3);
      dimPolynomial(a, va, b, vb, c, vc, 0,1,2, o4, p4, q4, r4);
      dimPolynomial(a, va, b, vb, c, vc, 1,0,2, o5, p5, q5, r5);
      T o = o0 - o1 + o2 - o3 + o4 - o5;
      p = p0 - p1 + p2 - p3 + p4 - p5;
      q = q0 - q1 + q2 - q3 + q4 - q5;
      r = r0 - r1 + r2 - r3 + r4 - r5;
      //std::cout << o << " " << p << " " << q << " " << r << std::endl;
      //o = o < 1e-12 ? 1e-12 : o;
      o += 1e-12;
      p/=o;
      q/=o;
      r/=o;
      //std::cout << o << " " << p << " " << q << " " << r << std::endl;
    }
    

    T fRand(){
      T f = (T)rand()/(T)RAND_MAX;
      return f;
    }

    T qRand(){
      T f = (T)rand()/(T)RAND_MAX;
      return 1.0 - 2.0*f;
    }

    void testCollision(){

      // we'll build a random triangle, and a random point from the interior of the triangle 
      
      coordinate_type T0(qRand(),qRand(),qRand());
      coordinate_type T1(qRand(),qRand(),qRand());
      coordinate_type T2(qRand(),qRand(),qRand());
      T q0 = 0.33*fRand(); 
      T q1 = 0.33*fRand();
      T q2 = 1.0 - q0 - q1;
      
      //coordinate_type T0(-1.0, 0.5, 0.2);
      //coordinate_type T1( 0.5, 1.0, 0.2);
      //coordinate_type T2( 1.0,-1.0,-0.2);
      
      //T q0 = 0.33; 
      //T q1 = 0.33;
      //T q2 = 1.0 - q0 - q1;
      
      coordinate_type x = q0*T0 + q1*T1 + q2*T2;
      coordinate_type tri[3] = {T0, T1, T2}; 
      
      m2::distance_calculator<SPACE> calc; 
      T dist = calc.distanceFromTriangle(tri, x);
      std::cout << " dist from tri: " << dist << std::endl;
      
      coordinate_type v0(qRand(),qRand(),qRand());
      coordinate_type v1(qRand(),qRand(),qRand());
      coordinate_type v2(qRand(),qRand(),qRand());
      coordinate_type vx(qRand(),qRand(),qRand());
      
      //coordinate_type v0(-0.1,-0.22,-0.1);
      //coordinate_type v1(-0.2,-0.1,-0.11);
      //coordinate_type v2(-0.1,-0.2, 0.3);
      //coordinate_type vx( 0.1, 0.1, 0.2);
      
      //choose a random time, displace the triangle and point a random vector away 
      T t = fRand();
      coordinate_type Tp0 = T0 + v0*t;
      coordinate_type Tp1 = T1 + v1*t;
      coordinate_type Tp2 = T2 + v2*t;
      coordinate_type xp = x + vx*t;
      
      //now run the continuous collision detection backwards and see if we get a negative time
      //as one of the roots
      
      T p,q,r, roots[3];
      buildPolynomial(Tp1 - Tp0, -v1 + v0, 
		      Tp2 - Tp0, -v2 + v0, 
		      xp  - Tp0, -vx + v0, p, q, r);
      int rootcount
	= magnet::math::cubicSolve(p,q,r, roots[0], roots[1], roots[2]);
      if(rootcount == 0)
	std::cout << " no roots  "  << std::endl;
      if(rootcount == 1)
	std::cout << " poly   : " << p << " " << q << " " << r << " t: " << t << " roots: " 
		  << roots[0] << std::endl;
      if(rootcount == 2)
	std::cout << " poly   : " << p << " " << q << " " << r << " t: " << t << " roots: " 
		  << roots[0] << " " << roots[1] <<std::endl;
      if(rootcount == 3)
	std::cout << " poly   : " << p << " " << q << " " << r << " t: " << t << " roots: " 
		  << roots[0] << " " << roots[1] << " " << roots[2] << std::endl;      
    }

    void firstImpulse(vector<coordinate_type> &  vels, T eps){
      TIMER functionTimer(__FUNCTION__);
      vector<contact_manifold> collisions;

      vector<vertex_ptr> & vertices = mMesh->get_vertices();
      vector<coordinate_type> points; 
      points.resize(vertices.size());
      for(int i = 0; i < vertices.size(); i++){
	points[i] = vertices[i]->coordinate();
      }
      vector<face_ptr>   & faces    = mMesh->get_faces();
      vector<int>           triIndex; triIndex.resize(3*faces.size());
      vector<triangle_type> tris;      tris.resize(faces.size());

      for(int i = 0; i < faces.size(); i++){
	tris[i][0] = faces[i]->fbegin()->coordinate();         //0
	tris[i][1] = faces[i]->fbegin()->next()->coordinate(); //1
	tris[i][2] = faces[i]->fbegin()->prev()->coordinate(); //2 {next->next == prev}

	triIndex[3*i+0] = faces[i]->fbegin()->vertex()->position_in_set();
	triIndex[3*i+1] = faces[i]->fbegin()->next()->vertex()->position_in_set();
	triIndex[3*i+2] = faces[i]->fbegin()->prev()->vertex()->position_in_set();
      }

      typedef aabb_tree<SPACE,triangle_type> triangle_tree;
      typedef typename triangle_tree::node_type edge_node_type;
      triangle_tree tree(tris);
      for(int i = 0; i < points.size(); i++){
	vector<int> collectedTris;
	vector<int> filteredTris;
	coordinate_type ci = points[i];
	swept_point_type li(ci, vels[i]); li.dt = dt;
	coordinate_type veli = vels[i];
	int indexMin = 0;
	coordinate_type cmin;
	T d = this->minLength;
	getAllNearest
	  <SPACE,triangle_type,coordinate_type>
	  (ci,tree,tris,collectedTris, 2.0*d);
	coordinate_type avgNormal(0.0,0.0,0.0);
	coordinate_type avgVelocity(0.0,0.0,0.0);
	T accumWeight = 0;
	
	for(int j = 0; j < collectedTris.size(); j++){
	  int pj = tree.permutation[j];
	  int j0 = triIndex[3*collectedTris[j] + 0];
	  int j1 = triIndex[3*collectedTris[j] + 1];
	  int j2 = triIndex[3*collectedTris[j] + 2];
	  if(j0 == i) continue;
	  if(j1 == i) continue;
	  if(j2 == i) continue;
	  vertex_ptr v4 = vertices[i];
	  vertex_ptr vj[3] = {vertices[j0],vertices[j1],vertices[j2]};
	  coordinate_type velj[3] = {vels[j0],vels[j1],vels[j2]};
	  
	  for(int k = 0; k < 3; k++){
	    int jk = triIndex[3*collectedTris[j] + k];
	    if(v4->shares_edge_with(vj[k])) continue;
	    if(jk == i) continue;
	    
	    coordinate_type cj = vertices[jk]->coordinate();
	    
	    T dist = norm2(ci-cj);
	    
	    T w = 1.0/powf(dist*dist + d*d,1.0);
	    avgNormal   += w*vj[k]->normal();
	    avgVelocity += w*velj[k];
	    accumWeight += w;
	  }
	}
	
	if(accumWeight > 1e-16){	  
	  avgNormal   /= accumWeight;
	  avgVelocity /= accumWeight;
	  
	  T J = (eps-d)/dt - dot(avgNormal, avgVelocity - veli);
	  std::cout << accumWeight << " " << avgVelocity 
		    << " " << dt << " " << J << std::endl;
	  vels[i] += J*avgNormal;
	}
      }
    }

    vector<contact_manifold> getCollisions(const vector<coordinate_type> &  vels){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & vertices = mMesh->get_vertices();
      vector<coordinate_type> points; points.resize(vertices.size());
      vector<contact_manifold> collisions;	      
      for(int i = 0; i < vertices.size(); i++){
	points[i] = vertices[i]->coordinate();
      }

      vector<face_ptr>   & faces    = mMesh->get_faces();
      vector<int>           triIndex; triIndex.resize(3*faces.size());
      vector<swept_triangle_type> tris;      tris.resize(faces.size());
	      
      for(int i = 0; i < faces.size(); i++){
	tris[i][0] = faces[i]->fbegin()->coordinate();         //0
	tris[i][1] = faces[i]->fbegin()->next()->coordinate(); //1
	tris[i][2] = faces[i]->fbegin()->prev()->coordinate(); //2 {next->next == prev}

	tris[i].v[0] = vels[faces[i]->fbegin()->vertex()->position_in_set()];         //0
	tris[i].v[1] = vels[faces[i]->fbegin()->next()->vertex()->position_in_set()]; //1
	tris[i].v[2] = vels[faces[i]->fbegin()->prev()->vertex()->position_in_set()]; //2 {next->next == prev}
	tris[i].dt = this->dt;

	triIndex[3*i+0] = faces[i]->fbegin()->vertex()->position_in_set();
	triIndex[3*i+1] = faces[i]->fbegin()->next()->vertex()->position_in_set();
	triIndex[3*i+2] = faces[i]->fbegin()->prev()->vertex()->position_in_set();
      }

      typedef aabb_tree<SPACE,swept_triangle_type> triangle_tree;
      typedef typename triangle_tree::node_type edge_node_type;
      triangle_tree tree(tris);
      //tid = omp_get_thread_num();

      for(int i = 0; i < points.size(); i++){
	vector<int> collectedTris;
	vector<int> filteredTris;
	T minDist = 99999;
	coordinate_type ci = points[i];
	swept_point_type li(ci, vels[i]); li.dt = dt;
	coordinate_type veli = velocities[i];
	int indexMin = 0;
	coordinate_type cmin;
	getAllNearest
	  <SPACE,swept_triangle_type,swept_point_type>
	  (li,tree,tris,collectedTris, 0.025*this->minLength);
	
	for(int j = 0; j < collectedTris.size(); j++){
	  int pj = tree.permutation[j];
	  if(triIndex[3*collectedTris[j] + 0] == i) continue;
	  if(triIndex[3*collectedTris[j] + 1] == i) continue;
	  if(triIndex[3*collectedTris[j] + 2] == i) continue;
	  for(int k = 0; k < 3; k++){
	    int index           = triIndex[3*collectedTris[j] + k];
	    coordinate_type cj = vertices[index]->coordinate();
	    T dist = norm2(ci-cj);
	    if(dist < minDist){
	      minDist = dist;
	      indexMin = index;
	      cmin = cj;
	    }
	  }
	}
	
	if(!vertices[i]) continue;
	if(!vertices[indexMin]) continue;
	if(vertices[indexMin]->flag > 0 || vertices[i]->flag > 0) continue;
	vertices[indexMin]->flag += 1;
	vertices[i]->flag += 1;
	cmin = vertices[indexMin]->coordinate();
	vertices[i]->update_normal();
	coordinate_type N = vertices[i]->normal();
	coordinate_type dc = ci-cmin; dc.normalize();
	if(minDist < 999){
	  //if(minDist < 999 && abs(dot(dc,N)) > 0.0){
	  contact_manifold cm; 
	  cm.mMesh = mMesh;
	  cm.p1 = i;     cm.p2 = indexMin;
	  cm.p1type = 0; cm.p2type = 0;
	  collisions.push_back(cm);
	}
      }
	    
      return collisions;
    }

    vector<contact_manifold> getCollisions(){
      TIMER functionTimer(__FUNCTION__);
      vector<contact_manifold> collisions;
      std::vector<contact_manifold* > local; local.resize(8);
      std::vector<int > numCollisions; numCollisions.resize(8);
#if 0
#pragma omp parallel
      {
	int np = omp_get_num_threads();
	int tid  = omp_get_thread_num();
#endif
	std::vector<contact_manifold> localCollisions;
	      
	vector<vertex_ptr> & vertices = mMesh->get_vertices();
	vector<coordinate_type> points; points.resize(vertices.size());
	      
	//#pragma omp parallal for
	for(int i = 0; i < vertices.size(); i++){
	  points[i] = vertices[i]->coordinate();
	}
	vector<face_ptr>   & faces    = mMesh->get_faces();
	vector<int>           triIndex; triIndex.resize(3*faces.size());
	vector<triangle_type> tris;      tris.resize(faces.size());

	//#pragma omp parallal for
	for(int i = 0; i < faces.size(); i++){
	  tris[i][0] = faces[i]->fbegin()->coordinate();         //0
	  tris[i][1] = faces[i]->fbegin()->next()->coordinate(); //1
	  tris[i][2] = faces[i]->fbegin()->prev()->coordinate(); //2 {next->next == prev}

	  triIndex[3*i+0] = faces[i]->fbegin()->vertex()->position_in_set();
	  triIndex[3*i+1] = faces[i]->fbegin()->next()->vertex()->position_in_set();
	  triIndex[3*i+2] = faces[i]->fbegin()->prev()->vertex()->position_in_set();
	}

	typedef aabb_tree<SPACE,triangle_type> triangle_tree;
	typedef typename triangle_tree::node_type edge_node_type;
	triangle_tree tree(tris);

	//#pragma omp parallal for
	for(int i = 0; i < points.size(); i++){
#if 0
	  int tid = omp_get_thread_num();
#endif
	  vector<int> collectedTris;
	  T minDist = 99999;
	  coordinate_type ci = points[i];
	  int indexMin = 0;
	  coordinate_type cmin;
	  getAllNearestTriPoint<SPACE>(ci,tree,tris,collectedTris, tol);
	  for(int j = 0; j < collectedTris.size(); j++){
	    int pj = tree.permutation[j];
	    if(triIndex[3*collectedTris[j] + 0] == i) continue;
	    if(triIndex[3*collectedTris[j] + 1] == i) continue;
	    if(triIndex[3*collectedTris[j] + 2] == i) continue;
	    for(int k = 0; k < 3; k++){
	      int index           = triIndex[3*collectedTris[j] + k];
	      coordinate_type cj = vertices[index]->coordinate();
	      T dist = norm2(ci-cj);
	      if(dist < minDist){
		minDist = dist;
		indexMin = index;
		cmin = cj;
	      }
	    }
	  }
	
	  if(!vertices[i]) continue;
	  if(!vertices[indexMin]) continue;
	  if(vertices[indexMin]->flag > 0 || vertices[i]->flag > 0) continue;
	  vertices[indexMin]->flag += 1;
	  vertices[i]->flag += 1;
	  cmin = vertices[indexMin]->coordinate();
	  vertices[i]->update_normal();
	  coordinate_type N = vertices[i]->normal();
	  coordinate_type dc = ci-cmin; dc.normalize();

	  if(minDist < 999 && abs(dot(dc,N)) > 0.0){
	    contact_manifold cm; 
	    cm.mMesh = mMesh;
	    cm.p1 = i;     cm.p2 = indexMin;
	    cm.p1type = 0; cm.p2type = 0;
	    //localCollisions.push_back(cm);
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
    
    bool checkEdgeInSet(face_vertex_ptr fv0,
			face_vertex_ptr fv1,
			face_vertex_ptr fva[4]){
      vertex_ptr v0 = fv0->vertex();
      vertex_ptr v1 = fv1->vertex();
      for(int i = 0; i < 4; i++){
	vertex_ptr v0i = fva[i]->vertex();
	vertex_ptr v1i = fva[i]->next()->vertex();
	if(v0 == v0i && v1 == v1i) return true;
	if(v1 == v0i && v0 == v1i) return true;
      }
      return false;	
    }

    triangle_type makeTriangle(face_vertex_ptr v0,
			       face_vertex_ptr v1,
			       face_vertex_ptr v2){
      //convenience triangle constructor
      return triangle_type(v0->vertex()->coordinate(),
			   v1->vertex()->coordinate(),
			   v2->vertex()->coordinate());
    };

    face_ptr makeFace(vertex_ptr v0,
		      vertex_ptr v1,
		      vertex_ptr v2){
  
      face_vertex_ptr fv0 = new face_vertex_type();
      face_vertex_ptr fv1 = new face_vertex_type();
      face_vertex_ptr fv2 = new face_vertex_type();
      fv0->next() = fv1;   fv1->prev() = fv0;
      fv1->next() = fv2;   fv2->prev() = fv1;
      fv2->next() = fv0;   fv0->prev() = fv2;
      fv0->vertex() = v0; v0->add_face_vertex(fv0);
      fv1->vertex() = v1; v1->add_face_vertex(fv1);
      fv2->vertex() = v2; v2->add_face_vertex(fv2);
      face_ptr f = new face_type();
      f->setHead(fv0);
      fv0->face() = f;
      fv1->face() = f;
      fv2->face() = f;
      f->update_all();
      return f;
    }

    void stitchFaces(control_ptr control, face_ptr f1, face_ptr f2){
      face_vertex_ptr fv1Minb = f1->fbegin();
      face_vertex_ptr fv2Minb = f2->fbegin();
      bool iteratingOuter = true;
      face_vertex_ptr fv2Compb = f2->fbegin();
      face_vertex_ptr fv2Compe = fv2Compb->next();
      T distSqMax = 999999;
      vector<coordinate_type> minEdgeSet;
      while (iteratingOuter){
	face_vertex_ptr fv1b = f1->fbegin();
	face_vertex_ptr fv1e = f1->fend();
	face_vertex_ptr fv2b = fv2Compb;
	face_vertex_ptr fv2e = fv2Compb->next();
	T distSqSum = 0;
	bool iteratingInner = true;
	int i = 0; int j = 0;

	vector<coordinate_type> curEdgeSet;
	while (iteratingInner){

	  distSqSum += norm2(fv1b->coordinate() - fv2b->coordinate()); 
	  T distSqPrev = norm2(fv1b->coordinate() - fv2b->prev()->coordinate()); 
	  T distSqNext = norm2(fv1b->next()->coordinate() - fv2b->coordinate()); 
	  if(fv1b != fv1e && fv2b != fv2e){
	    if(distSqNext < distSqPrev){
	      curEdgeSet.push_back(fv1b->next()->coordinate());
	      curEdgeSet.push_back(fv2b->coordinate());
	      fv1b = fv1b->next();
	    }
	    else{
	      curEdgeSet.push_back(fv1b->coordinate());
	      curEdgeSet.push_back(fv2b->prev()->coordinate());
	      fv2b = fv2b->prev();
	    }
	  }
	  else if(fv2b == fv2e && fv1b != fv1e){
	    fv1b = fv1b->next();
	  }
	  else if(fv1b == fv1e && fv2b != fv2e){
	    fv2b = fv2b->prev();
	  }
	  else{
	    T distSqPrev = norm2(fv1b->coordinate() - fv2b->prev()->coordinate()); 
	    T distSqNext = norm2(fv1b->next()->coordinate() - fv2b->coordinate()); 
	    if(distSqNext < distSqPrev){
	      curEdgeSet.push_back(fv1b->next()->coordinate());
	      curEdgeSet.push_back(fv2b->coordinate());
	      distSqSum += distSqNext;
	    }
	    else{
	      curEdgeSet.push_back(fv1b->coordinate());
	      curEdgeSet.push_back(fv2b->prev()->coordinate());
	      distSqSum += distSqPrev;
	    }
	    iteratingInner = false;
	  }
	}
	if(distSqSum < distSqMax){
	  minEdgeSet = curEdgeSet;
	  distSqMax = distSqSum;
	  fv2Minb = fv2Compb;
	}
	iteratingOuter = fv2Compb != fv2Compe;
	fv2Compb = fv2Compb->prev();
      }
#if 0
      for(int i = 0; i < minEdgeSet.size(); i+=2){
	m2::Debugger& debug = m2::Debugger::get_instance();
	debug.DebugLines0.push_back(minEdgeSet[i+0]);
	debug.DebugLines0.push_back(minEdgeSet[i+1]);
      }
#endif 
      for(int i = 0; i < minEdgeSet.size(); i+=2){
	coordinate_type ep0 = minEdgeSet[i+0];
	coordinate_type ep1 = minEdgeSet[i+1];
	if(m2::intersectLineTest<SPACE,triangle_type>(ep0, ep1,exteriorTree,oldFaces))  return;
	if(m2::intersectLineTest<SPACE,triangle_type>(ep0, ep1,interiorTree,newFaces)) return;
      }
      face_vertex_ptr fv1b = f1->fbegin();
      face_vertex_ptr fv1e = f1->fend();
      face_vertex_ptr fv2b = fv2Minb;
      face_vertex_ptr fv2e = fv2Minb->next();
      bool iterating = true;
      //m2::Debugger& debug = m2::Debugger::get_instance();
      vector<face_vertex_ptr> newTriangles;
      while (iterating){      
	T distSqPrev = norm2(fv1b->coordinate() - fv2b->prev()->coordinate()); 
	T distSqNext = norm2(fv1b->next()->coordinate() - fv2b->coordinate()); 
	if(fv1b != fv1e && fv2b != fv2e){
	  if(distSqNext < distSqPrev){
	    newTriangles.push_back(fv1b->next());
	    newTriangles.push_back(fv2b);
	    newTriangles.push_back(fv1b);

	    fv1b = fv1b->next();
	  }
	  else{
	    newTriangles.push_back(fv2b->prev());
	    newTriangles.push_back(fv2b);
	    newTriangles.push_back(fv1b);
	    fv2b = fv2b->prev();
	  }
	}
	else if(fv2b == fv2e && fv1b != fv1e){
	  newTriangles.push_back(fv1b->next());
	  newTriangles.push_back(fv2b);
	  newTriangles.push_back(fv1b);

	  fv1b = fv1b->next();
	}
	else if(fv1b == fv1e && fv2b != fv2e){
	  newTriangles.push_back(fv2b->prev());
	  newTriangles.push_back(fv2b);
	  newTriangles.push_back(fv1b);
	  fv2b = fv2b->prev();
	}
	else{
	  T distSqPrev = norm2(fv1b->coordinate() - fv2b->prev()->coordinate()); 
	  T distSqNext = norm2(fv1b->next()->coordinate() - fv2b->coordinate()); 
	  if(distSqNext < distSqPrev){
	    newTriangles.push_back(fv1b->next());
	    newTriangles.push_back(fv2b);
	    newTriangles.push_back(fv1b);
	    newTriangles.push_back(fv2b->prev());
	    newTriangles.push_back(fv2b);	
	    newTriangles.push_back(fv1b->next());

	  }
	  else{
	    newTriangles.push_back(fv2b->prev());
	    newTriangles.push_back(fv2b);
	    newTriangles.push_back(fv1b);
	    newTriangles.push_back(fv1b->next());
	    newTriangles.push_back(fv2b->prev());
	    newTriangles.push_back(fv1b);
	  }
	  iterating = false;
	}
      }
      //this is the point where you would dump out all the new potential edges and test 
      //against the existing set of faces to see if there are any conflicts.
      
      int Nbefore = newFaces.size();
      for(int i = 0; i < newTriangles.size(); i+=3){
	triangle_type tri = makeTriangle(newTriangles[i+0],
					 newTriangles[i+1],
					 newTriangles[i+2]);
	newFaces.push_back(tri);
      }
      int Nafter = newFaces.size();
      for(int i = Nbefore; i < Nafter; i++){
	interiorTree.insert(newFaces,i,10);
      }

      vector<face_ptr> localNewFaces;
      for(int i = 0; i < newTriangles.size(); i+=3){
	face_ptr fi = makeFace(newTriangles[i+0]->vertex(),
			       newTriangles[i+1]->vertex(),
			       newTriangles[i+2]->vertex());
	control->push_face(fi);
	localNewFaces.push_back(fi);
      }
      for(int i = 0; i < newTriangles.size(); i+=3){
	int ip  = ((i - 3) + newTriangles.size())%newTriangles.size();
	coordinate_type ci0 = newTriangles[i+0]->vertex()->coordinate();
	coordinate_type ci1 = newTriangles[i+1]->vertex()->coordinate();
	coordinate_type ci2 = newTriangles[i+2]->vertex()->coordinate();
	coordinate_type cp0 = newTriangles[ip+0]->vertex()->coordinate();
	coordinate_type cp1 = newTriangles[ip+1]->vertex()->coordinate();
	coordinate_type cp2 = newTriangles[ip+2]->vertex()->coordinate();
	face_ptr fi = localNewFaces[i/3];
	face_ptr fpi = localNewFaces[ip/3];
	int ctype = 0;
	if(newTriangles[i+0]->face() == newTriangles[i+1]->face()) ctype = 0;
	else ctype = 1;

	if(ctype == 0){
	  face_vertex_ptr fv0 = fi->fbegin();
	  face_vertex_ptr fvd = newTriangles[i+0]; //discard
	  //discard we're getting the other so the index offset is the same
	  edge_ptr e = fvd->edge(); //retrieving old edge don't have to push
	  face_vertex_ptr fv1 = e->other(fvd);
	  fv0->edge() = e; e->v1() = fv0;
	  fv1->edge() = e; e->v2() = fv1;
	}
	else{
	  face_vertex_ptr fv0 = fi->fbegin()->prev();
	  face_vertex_ptr fvd = newTriangles[i+2]; 
	  //discard we're getting the other so the index offset is the same
	  edge_ptr e = fvd->edge(); //retrieving old edge don't have to push
	  face_vertex_ptr fv1 = e->other(fvd);
	  fv0->edge() = e; e->v1() = fv0;
	  fv1->edge() = e; e->v2() = fv1;
	}

	int ptype = 0;
	if(newTriangles[ip+0]->face() == newTriangles[ip+1]->face()) ptype = 0;
	else ptype = 1;    
    
	if(ctype == 0 && ptype == 0){
	  face_vertex_ptr fv1 = fi->fbegin()->next();
	  face_vertex_ptr fv2 = fpi->fbegin()->prev();
	  edge_ptr e = new edge_type();
	  e->v1() = fv1; fv1->edge() = e;
	  e->v2() = fv2; fv2->edge() = e;
	  control->push_edge(e); //making new edge have to push it
	}
	if(ctype == 0 && ptype == 1){
	  face_vertex_ptr fv1 = fi->fbegin()->next();
	  face_vertex_ptr fv2 = fpi->fbegin();
	  edge_ptr e = new edge_type();
	  e->v1() = fv1; fv1->edge() = e;
	  e->v2() = fv2; fv2->edge() = e;
	  control->push_edge(e);
	}
	if(ctype == 1 && ptype == 0){
	  face_vertex_ptr fv1 = fi->fbegin()->next();
	  face_vertex_ptr fv2 = fpi->fbegin()->prev();
	  edge_ptr e = new edge_type();
	  e->v1() = fv1; fv1->edge() = e;
	  e->v2() = fv2; fv2->edge() = e;
	  control->push_edge(e);
	}
	if(ctype == 1 && ptype == 1){
	  face_vertex_ptr fv1 = fi->fbegin()->next();
	  face_vertex_ptr fv2 = fpi->fbegin();
	  edge_ptr e = new edge_type();
	  e->v1() = fv1; fv1->edge() = e;
	  e->v2() = fv2; fv2->edge() = e;
	  control->push_edge(e);
	}
      }
      vector<face_vertex_ptr>  fva1 = f1->face_vertex_trace();
      vector<face_vertex_ptr>  fva2 = f2->face_vertex_trace();
      for(int i = 0; i < fva1.size(); i++) {
	fva1[i]->vertex()->remove_face_vertex(fva1[i]);
	delete fva1[i];
      }

      for(int i = 0; i < fva2.size(); i++) {
	fva2[i]->vertex()->remove_face_vertex(fva2[i]);
	delete fva2[i];
      }

      control->remove_face(f1->position_in_set());
      control->remove_face(f2->position_in_set());

      //delete f1;
      //delete f2;
    }
    
    void pipeEdge(edge_ptr e0, edge_ptr e1,
		  T d00, T d01,
		  T d10, T d11,
		  dynamic_octree<SPACE,triangle_type> & faceTree, 
		  vector<triangle_type> & corners){
      //if(e0->v1()->vertex() == e1->v1()->vertex()) return;
      //if(e0->v1()->vertex() == e1->v2()->vertex()) return;
      //if(e0->v2()->vertex() == e1->v1()->vertex()) return;
      //if(e0->v2()->vertex() == e1->v2()->vertex()) return;
      coordinate_type c00 = e0->v1()->vertex()->coordinate();
      coordinate_type c01 = e0->v2()->vertex()->coordinate();
      coordinate_type c10 = e1->v1()->vertex()->coordinate();
      coordinate_type c11 = e1->v2()->vertex()->coordinate();
      coordinate_type cAvg = 0.25*(c00+c01+c10+c11);
      m2::Debugger& debug = m2::Debugger::get_instance();
      //std::cout << " current join: " << debug.labelId << std::endl; 
      debug.add_id(cAvg[0],cAvg[1],cAvg[2]);


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
      for(int i = 0; i < 8; i++){
	if(e0->v1() == fva[i]) return;
	if(e0->v2() == fva[i]) return;
	if(e1->v1() == fva[i]) return;
	if(e1->v2() == fva[i]) return;
      }
      //std::cout << e0->v1() << " " << e0->v2() << " " << e1->v1() << " " << e1->v2() << std::endl;
      //std::cout << fv00 << " " << fv11 << " " << fv01 << " " << fv10 << std::endl;
      //std::cout << fv00p << " " << fv11p << " " << fv01p << " " << fv10p << std::endl;

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

      int jb = 0; T minE = 999; //beginning i and j indices
      for(int i = 0; i < 4; i++){
	fva0[i]->vertex()->flag = 1;
	fva1[i]->vertex()->flag = 1;
      }
#if 1
      for(int i = 0; i < 4; i++){
	T E = 0;
	for(int k = 0; k < 4; k++){
	  int km = (k-1+4)%4;
	  int j  = (i-k+4)%4;
	  int jm = (i-k+1+4)%4;
	  coordinate_type cik  = fva0[k]->coordinate();
	  coordinate_type cimk = fva0[km]->coordinate();
	  coordinate_type cjk  = fva1[j]->coordinate();

	  bool checkComp0 = checkEdgeInSet(fva0[k],fva1[j],fva0);
	  bool checkComp1 = checkEdgeInSet(fva0[k],fva1[j],fva1);
	  if(!checkEdgeInSet(fva0[k],fva1[j],fva0) && 
	     !checkEdgeInSet(fva0[k],fva1[j],fva1)){
	    coordinate_type dc  = cik-cjk;
	    E += dot(dc,dc);
	  }
#if 1
	  if(!checkEdgeInSet(fva0[km],fva1[j],fva0) && 
	     !checkEdgeInSet(fva0[km],fva1[j],fva1)){
	    coordinate_type dcm = cimk-cjk;
	    E += dot(dcm,dcm);
	  }
#endif
	}
	if(E < minE){
	  minE = E;    jb = i;
	}
      }

#endif

      vector<int> ppairs;
      for(int i = 0; i < 4; i++){
	int ii = (i)%4;
	int jj = (jb-i+4)%4;
	int im = (i-1+4)%4;
	int jm = (jb-i+1+4)%4;
	if(!checkEdgeInSet(fva0[ii],fva1[jj],fva0) && 
	   !checkEdgeInSet(fva0[ii],fva1[jj],fva1)){
	  ppairs.push_back(ii);
	  ppairs.push_back(jj);
	}
#if 1

	if(!checkEdgeInSet(fva0[im],fva1[jj],fva0) && 
	   !checkEdgeInSet(fva0[im],fva1[jj],fva1)){
	  ppairs.push_back(im);
	  ppairs.push_back(jj);
	}
#endif
      }
      std::cout << std::endl;
      bool hit = false;
      for(int i = 0; i < ppairs.size()/2; i++){
	face_vertex_ptr fv0 = fva0[ppairs[2*i+0]];
	face_vertex_ptr fv1 = fva1[ppairs[2*i+1]];
	
	coordinate_type ep0 = fv0->coordinate();
	coordinate_type ep1 = fv1->coordinate();
	if(m2::intersectLineTest<SPACE,triangle_type>(ep0, ep1,faceTree,corners)) hit = true;
	if (hit) return;
      }
      //the runs that made it through the gauntlet
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

      face_ptr f0 = cons.delete_edge(mMesh,e0);
      face_ptr f1 = cons.delete_edge(mMesh,e1);
      vector<edge_ptr> collectedEdges;
      for(int i = 0; i < ppairs.size()/2; i++){
	face_vertex_ptr fv0 = fva0[ppairs[2*i+0]];
	face_vertex_ptr fv1 = fva1[ppairs[2*i+1]];
	m2::Debugger& debug = m2::Debugger::get_instance();
	debug.DebugLines.push_back(fv0->coordinate());
	debug.DebugLines.push_back(fv1->coordinate());
	edge_ptr ei = cons.insert_edge(mMesh, fv0, fv1);
	collectedEdges.push_back(ei);
      }
      
      for(int i = 0; i < collectedEdges.size(); i++){
	edge_ptr ei = collectedEdges[i];
	face_vertex_ptr fv0 = ei->v1();
	face_vertex_ptr fv1 = ei->v2();
	// m2::Debugger& debug = m2::Debugger::get_instance();      
	//debug.DebugTriangles.push_back(fv0->coordinate());
	//debug.DebugTriangles.push_back(fv0->next()->coordinate());
	//debug.DebugTriangles.push_back(fv0->next()->next()->coordinate());
	//debug.DebugTriangles.push_back(fv1->coordinate());
	//debug.DebugTriangles.push_back(fv1->next()->coordinate());
	//debug.DebugTriangles.push_back(fv1->next()->next()->coordinate());
	corners.push_back(makeTriangle(fv0, fv0->next(),fv0->next()->next()));
	corners.push_back(makeTriangle(fv1, fv1->next(),fv1->next()->next()));
	faceTree.insert(corners, corners.size()-1, 8);
	faceTree.insert(corners, corners.size()-2, 8);
	//faceTree.insert(makeTriangle(fv1, fv1->next(),fv1->next()->next()));
	// if(fv0->face()->size() > 3){
	//   face_vertex_ptr fvn = fv0->next()->next();
	//   debug.DebugTriangles.push_back(fvn->coordinate());
	//   debug.DebugTriangles.push_back(fvn->next()->coordinate());
	//   debug.DebugTriangles.push_back(fvn->next()->next()->coordinate());
	//   insertedTriangles.push_back(triangle(fvn, fvn->next(),fvn->next()->next()));
	// }
      }
      
    }

    bool join(){
      bool topology_change = false;
      vector<contact_manifold> collisions;
      if(this->velocities.size() > 0)
	collisions = this->getCollisions(this->velocities);
      else
	collisions = this->getCollisions();
      vector<vertex_ptr> verts = mMesh->get_vertices();

      for(int i = 0; i < collisions.size(); i++){
	coordinate_type c1 = verts[collisions[i].p1]->coordinate();
	coordinate_type c2 = verts[collisions[i].p2]->coordinate();
	Debugger& debug = Debugger::get_instance();
	debug.DebugPoints.push_back(c1);
	debug.DebugPoints.push_back(c2);
	debug.DebugLines.push_back(c1);
	debug.DebugLines.push_back(c2);
      }
      if(collisions.size() > 0) topology_change = true;
      
      std::cout << " - joining: " << collisions.size() << " pairs" << std::endl;
      mContactSorter.mMesh = mMesh;
      std::sort(collisions.begin(), collisions.end(), mContactSorter);
#if 1
      if(collisions.size() == 0) return false;
      int sharedEdges = 0;
      int sepEdges = 0;
      for(int i = 0; i < collisions.size(); i++){
	vertex_ptr v1 = verts[collisions[i].p1];
	vertex_ptr v2 = verts[collisions[i].p2];
	construct<SPACE> cons;
	//std::cout << v1->position_in_set() << " " << v1->size() << " | "
	//	  << v2->position_in_set() << " " << v2->size() << " " << std::endl;
	if(v1 == v2)  continue;
	if(v1->size() < 3) continue;
	if(v2->size() < 3) continue;

	if(v1->shares_edge_with(v2)){	  
	  edge_ptr e = v1->get_shared_edge(v2);
	  cons.collapse_edge(mMesh,e);
	  sharedEdges++;
	}
	else{
	  //we have to deal with two odd topologies I can't account for, yet
	  // f1 and f2 somehow are equal, yet they didn't start from the same valence set
	  // v2 is somehow in f1 after the vertex deletion, so we abort and let triangulation reclaim the 
	  // the new face
	  face_ptr f1 = cons.delete_vertex_primitive(mMesh, v1);
	  if(!f1) continue;
	  if(f1->has_vertex(v2)) continue;
	  face_ptr f2 = cons.delete_vertex_primitive(mMesh, v2);
	  if(!f2) continue;
	  if(f1 != f2)
	    stitchFaces(mMesh,f1,f2);
	  sepEdges++;
	}
      }
      
      std::cout << "  -collapse: " << sharedEdges << " edges, and piped: " << sepEdges << " edges" << std::endl;
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
    control_ptr		 mMesh;
    dynamic_octree<SPACE,triangle_type> exteriorTree;
    dynamic_octree<SPACE,triangle_type> interiorTree;
    vector<triangle_type> oldFaces;
    vector<triangle_type> newFaces;
    vector<coordinate_type> velocities;
    T dt;
    T minLength;
  };

  template <typename SPACE>
  class calcSdf{ M2_TYPEDEFS
  public:
    typedef aabb_tree<SPACE,triangle_type> triangle_tree;
    typedef typename triangle_tree::node_type edge_node_type;
    control_ptr mMesh;

    vector<coordinate_type> normals;
    vector<coordinate_type> vertices;
    vector<triangle_type>  tris;
    vector<int>            triIndices;      
    triangle_tree tree =   triangle_tree(this->tris);

    calcSdf(control_ptr targetMesh){
      this->mMesh = targetMesh;
    };

    void initialize(){
      std::cout << " in init sequence" << std::endl;
      this->normals = mMesh->get_normals();
      this->vertices  = mMesh->get_coordinates();

      this->tris       = this->getTris(mMesh);
      this->triIndices = this->getTriIndices(mMesh);      
      this->tree = triangle_tree(this->tris);
      
    }

    vector<triangle_type> getTris(control_ptr targetMesh){
      vector<face_ptr>   & faces    = targetMesh->get_faces();
      vector<triangle_type> ltris;      
      ltris.resize(faces.size());
	for(int i = 0; i < faces.size(); i++){
	ltris[i][0] = faces[i]->fbegin()->coordinate();         //0
	ltris[i][1] = faces[i]->fbegin()->next()->coordinate(); //1
	ltris[i][2] = faces[i]->fbegin()->prev()->coordinate(); //2 {next->next == prev}
      }
      return ltris;
    }
    
    vector<int> getTriIndices(control_ptr targetMesh){
      vector<face_ptr>   & faces    = targetMesh->get_faces();
      vector<int>           ltriIndex; 
      ltriIndex.resize(3*faces.size());
      for(int i = 0; i < faces.size(); i++){
	ltriIndex[3*i+0] = faces[i]->fbegin()->vertex()->position_in_set();
	ltriIndex[3*i+1] = faces[i]->fbegin()->next()->vertex()->position_in_set();
	ltriIndex[3*i+2] = faces[i]->fbegin()->prev()->vertex()->position_in_set();
      }
      return ltriIndex;
    }

    inline T clamp(T x, T a, T b){
        return x < a ? a : (x > b ? b : x);
    };

    inline void  distanceFromTriangle(coordinate_type p, int tri,
				  T & d, coordinate_type & Ni)
    {
      int i0 = this->triIndices[3*tri+0];
      int i1 = this->triIndices[3*tri+1];
      int i2 = this->triIndices[3*tri+2];
      const coordinate_type& v0 = this->vertices[i0];
      const coordinate_type& v1 = this->vertices[i1];
      const coordinate_type& v2 = this->vertices[i2];

      coordinate_type u = v1-v0;
      coordinate_type v = v2-v0;
      coordinate_type N = cross(u,v);
      T iN2 = 1.0/(dot(N,N));
      coordinate_type w = p - v0;
      T b10 = dot(cross(u,w),N)*iN2;
      T b20 = dot(cross(w,v),N)*iN2;
      T b12 = 1.0 - b10 - b20;
      
      if(b10 <= 0){
	//line - v0, v1
	b20 = dot((p-v0),u)/(dot(u,u));
	b20 = clamp(b20, 0.0, 1.0);
	b12 = 1.0 - b20;
	b10 = 0;
      }    
      else if(b20 <= 0){
	//line - v0, v2
	b10 = dot((p-v0),v)/(dot(v,v));
	b10 = clamp(b10, 0.0, 1.0);
	b12 = 1.0 - b10;
	b20 = 0;
      }
      else if (b12 <= 0)
	{
	//line - v1, v2
	coordinate_type x = v2 - v1;
	b10 = dot((p-v1),x)/(dot(x,x));
	b10 = clamp(b10, 0.0, 1.0);
	b20 = 1.0 - b10;
	b12 = 0;
      }
      
      coordinate_type c  = b10*v2 + b20*v1 + b12*v0;
      //std::cout << tri << " | " << c << " | " << b10 << " " << b20 << " " << b12 << " " << b10 + b20 + b12 << std::endl;
      
      m2::Debugger& debug = m2::Debugger::get_instance();
      //debug.DebugLines.push_back(p);
      //debug.DebugLines.push_back(c);
      Ni = b12*this->normals[i0] + b20*this->normals[i1] + b10*this->normals[i2];
      Ni.normalize();
      d = dot(p-c, Ni);
      
    }
    
    vector<pair<T,coordinate_type> > getSdf(vector<coordinate_type> & coordinates){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & vertices = mMesh->get_vertices();
      vector<coordinate_type> points; points.resize(coordinates.size());
      vector<pair<T,coordinate_type> > out;
      out.resize(coordinates.size());
	      
      //#pragma omp parallal for
      for(int i = 0; i < coordinates.size(); i++){
	int cId = 
	  getNearest<SPACE, triangle_type>
	  (coordinates[i], this->tree, this->tris);
	coordinate_type N;
	T d;
	this->distanceFromTriangle(coordinates[i], cId, d, N);
	pair<T,coordinate_type> dN; 

	dN.first = d;
	dN.second = N;
	out[i] = dN;
      }
      return out;
    }
  };

  template <typename SPACE>
  class moving_mesh{ M2_TYPEDEFS;
    typedef coordinate_type	point_type;
    typedef typename SPACE::triangle_type	        triangle;
    typedef typename list<coordinate_type >::iterator	pl_iterator;
    typedef typename m2::Debugger Debugger;
  public:
    moving_mesh(control_ptr obj_in){
      mMesh = obj_in;
      maxCurvature = 3.0;
      minCurvature = 0.001;
      minLength = 0.01;
      minCollapseLength = 0.00025;
      maxLength = 0.06;

      coordinate_type cen = mMesh->calc_center();
      coordinate_type min = mMesh->calc_min();
      coordinate_type max = mMesh->calc_max();
      coordinate_type lengths = max-min;
      int minRes = 64;
      T minLength = lengths[0];
      minLength = minLength < lengths[1] ? minLength : lengths[1];
      minLength = minLength < lengths[2] ? minLength : lengths[2];

      T dx = minLength/(T)(minRes-1);
      remesh<SPACE> rem;
      rem.triangulate(mMesh);
      //construct<SPACE> cons;
      //face_ptr nf = cons.delete_vertex(mMesh,mMesh->get_vertices()[0]);
      //vector<T> edgeWeights;
      //vector<T> vertexWeights;
      //calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);
    }
    face_bin<SPACE> & getHash(){return mFaceHash;}

    void pin_half(){
      mMesh->pack();
      vertex_array& vl = mMesh->get_vertices();
      coordinate_type cen(0,0,0);
      for (int j = 0; j < vl.size(); j++) {
	cen += vl[j]->coordinate();
	vl[j]->pinned = false;
      }

      cen /= (T)vl.size();
      std::cout << " cen: " << cen << std::endl;
      for (int j = 0; j < vl.size(); j++) {
	if(vl[j]->coordinate()[1] < cen[1]){
	  vl[j]->pinned = true;
	}
      }
    }
    

    void drawField(){
      vector<vertex_ptr>& tverts = mMesh->get_vertices();
      if(tverts.size() > 0  && mDirectionWeights.size() > 0 && mDirectionField.size() > 0){
	for (long i = 0; i < tverts.size(); i++) {
	
	  coordinate_type p0 = tverts[i]->coordinate();
	  coordinate_type c = 0.0035*mDirectionWeights[i];
	  coordinate_type dp0 = p0 - c[0]*mDirectionField[i][0];
	  coordinate_type dp1 = p0 + c[0]*mDirectionField[i][0];
	  glBegin(GL_LINES);
	
	  glColor4f(0.25f, 0.25f, 0.8f, 0.5f);
	  glVertex3f(dp0[0], dp0[1], dp0[2]);
	  glVertex3f(dp1[0], dp1[1], dp1[2]);
	  glEnd();

	  dp0 = p0 - c[1]*mDirectionField[i][1];
	  dp1 = p0 + c[1]*mDirectionField[i][1];
	  glBegin(GL_LINES);	    
	  glColor4f(0.25f, 0.8f, 0.25, 0.5f);
	  glVertex3f(dp0[0], dp0[1], dp0[2]);
	  glVertex3f(dp1[0], dp1[1], dp1[2]);
	  glEnd();

	  dp0 = p0 - c[2]*mDirectionField[i][2];
	  dp1 = p0 + c[2]*mDirectionField[i][2];
	  glBegin(GL_LINES);	    
	  glColor4f(0.8f, 0.25f, 0.25, 0.5f);
	  glVertex3f(dp0[0], dp0[1], dp0[2]);
	  glVertex3f(dp1[0], dp1[1], dp1[2]);
	  glEnd();
	}
      }
    }

    coordinate_type getSvdNormal(vertex_ptr v){
      if(v->size() ==  0) return coordinate_type(0,0,0);
      face_vertex_ptr itb = v->fbegin();
      face_vertex_ptr ite = v->fend();
      bool at_head = false;
      coordinate_type kA(0,0,0);
      T wTotal = 0;
      coordinate_type c0 = v->coordinate();
      int k = 0;
      coordinate_type cov[3];
      while (!at_head) {
	at_head = itb==ite;
	T wij = itb->face()->area();
	coordinate_type c1 = itb->next()->coordinate();
	wTotal +=  wij;
	kA += wij*(c1-c0);
	
	coordinate_type dc = itb->face()->normal();
	for(int l = 0; l < 3; l++)
	  for(int m = 0; m < 3; m++){
	    cov[l][m] +=dc[l]*dc[m];
	  }
	itb = itb->vnext();
	k++;
      }
      
      if(wTotal > 1e-10){
	v->update_normal();
	coordinate_type w;
	calcSVD<SPACE>(cov, w);
	coordinate_type Nv = v->normal();
	coordinate_type N = cov[0];
	if(dot(Nv, N) < 0.0) return -N;
	else return N;
      }
      else return  v->normal();
    }

    vector<coordinate_type> getFaceOffsets(m2::control<SPACE>& in){
      TIMER functionTimer(__FUNCTION__);   
      vector<T> vertexWeights;
      vector<T> edgeWeights;
      vector<vertex_ptr>& tverts = in.get_vertices();
      vector<coordinate_type> normals;
      normals.resize(tverts.size());
      for (long i = 0; i < tverts.size(); i++) {
	if(!in.has_vertex(i)) continue;
	if(tverts[i]->pinned) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;
	normals[i] = getSvdNormal(v);
      }
      return normals;
    }


    coordinate_type getButterflyVertexWeight(face_vertex_ptr v){
      // \/\./
      //vn = next 
      //vnn = vertext_next next
      //pvpp = prev vertex_prev next
      //coordinate_type c = v->coordinate();
      // coordinate_type cn = v->next()->coordinate();
      // coordinate_type cvnn = v->vnext()->next()->coordinate();
      // coordinate_type cpvpp = v->prev()->vprev()->prev()->coordinate();
      // return (8.0*c + 2.0*cn - 1.0*(cvnn + cpvpp));
      face_vertex_ptr itb = v;
      face_vertex_ptr ite = v->vprev();
      bool iterating = true;
      int K = 0, j=0;
      while(iterating){
	iterating = itb != ite;
	itb = itb->vnext();
	K++;
      }
      itb = v;
      iterating = true;
      T s = 0;
      coordinate_type cs(0,0,0);
      while(iterating){
	iterating = itb != ite;
	itb = itb->vnext();
	T frac = (T)j/(T)K;
	T sj = (0.25 + cos(2.0*M_PI*frac) + 0.5*cos(4.0*M_PI*frac))/(T)K;
	coordinate_type ci = itb->next()->coordinate();
	cs += sj*ci;
	s += sj;
	j++;
      }      
      coordinate_type c = v->coordinate();
      return 0.75*c + cs;
    } 

    coordinate_type getButterflyWeight(edge_ptr e){
      face_vertex_ptr v1 = e->v1();
      face_vertex_ptr v2 = e->v2();
      coordinate_type c1 = getButterflyVertexWeight(v1);
      coordinate_type c2 = getButterflyVertexWeight(v2);
      return 0.5*(c1 + c2);
    }

    coordinate_type getAverageWeight(edge_ptr e){
      face_vertex_ptr v1 = e->v1();
      face_vertex_ptr v2 = e->v2();
      coordinate_type c1 = v1->coordinate();
      coordinate_type c2 = v2->coordinate();
      return 0.5*(c1 + c2);
    }

    void checkEdge(int edgeNum){
      edge_array & edges = mMesh->get_edges();
      if(edgeNum < edges.size()){
	edge_ptr ei  = mMesh->get_edges()[edgeNum];
	std::cout << "edge check: " 
		  << ei->v1()->vertex() << " "
		  << ei->v2()->vertex() << " "
		  << ei->v1()->vertex()->size() << " "
		  << ei->v2()->vertex()->size() << " " 
		  << ei->length() << std::endl; 
      }
    }

    void relax_mesh(){
      bool relaxing = true;
      int k = 0;
      TIMER functionTimer(__FUNCTION__);
      std::cout << " mesh size: " << mMesh->get_vertices().size() << std::endl;
      // relaxing = this->delete_vertices_curve();
      
      //relaxing = this->insert_edges();

      mMesh->colorVertices();
      std::cout << " max graph color: " << mMesh->maxGraphColor << std::endl;
      while(relaxing){
	relaxing = this->collapse_edges();
	
	delete_degenerate(); // 
      }

      relaxing = true;
      while(relaxing)
	relaxing = this->insert_edges();

      k = 0;
      relaxing = true;
      while(relaxing && k < 2){
	set_operations<SPACE> setOps;
	relaxing = setOps.flip_edges(mMesh);
	k++;
      }
      mesh_filter<SPACE> filter;
      for(int i = 0; i < 5; i++)
	filter.filter(*mMesh,0.01);
      mesh_calculator<SPACE> calc;
      //calc.shadeVertices(*mMesh);
      //23765
    }

    void cache_mesh(){
      face_array& faces = mMesh->get_faces();
      Debugger& debug = Debugger::get_instance();
      for(int i = 0; i < faces.size(); i++){
	debug.CachedTriangles.push_back(faces[i]->fbegin()->vertex()->coordinate());
	debug.CachedTriangles.push_back(faces[i]->fbegin()->next()->vertex()->coordinate());
	debug.CachedTriangles.push_back(faces[i]->fbegin()->next()->next()->vertex()->coordinate());
      }
    }

    void relax_mesh_no_collapse(){
      bool relaxing = true;
      int k = 0;
      TIMER functionTimer(__FUNCTION__);
      std::cout << " mesh size: " << mMesh->get_vertices().size() << std::endl;
      // for(int i = 0; i < 10; i++){
      // 	;
      // }

      Debugger& debug = Debugger::get_instance();
      debug.DebugPoints.clear();
      debug.DebugBoxes.clear();
      debug.DebugLines.clear();
      debug.DebugLines0.clear();
      debug.DebugLines1.clear();
      debug.DebugTriangles.clear();
      debug.DebugTriangles0.clear();
      debug.CachedTriangles.clear();

      int i = 0;
      relaxing = this->collapse_edges();

      // relaxing = this->collapse_low_curve();
      //std::cout << " deleting degenerates" << std::endl;
      i = 0;
      while(relaxing && i < 20){
	relaxing = delete_degenerate();
      }
      relaxing = true;
      while(relaxing)
      	relaxing = this->insert_edges();
      i = 0;

      while(relaxing && i < 20){
	relaxing = delete_degenerate();
      }

      k = 0;
      relaxing = true;

      remesh<SPACE> rem;
      //
      k = 0;
      relaxing = true;
      while(relaxing && k < 0){
      	set_operations<SPACE> setOps;
      	relaxing = setOps.flip_edges(mMesh);
      	k++;
      }
      this->cache_mesh();

      this->hashFaces();
      mMesh->pack();

      mesh_filter<SPACE> filter;
      mMesh->update_all();
      for(int i = 0; i < 2; i++){
      	//filter.taubinFilter(*mMesh,0.1,0.11);
      	filter.nullLaplacianFilter(*mMesh,0.5);
      }

      vector<edge_ptr>   & edges = mMesh->get_edges();
      vector<coordinate_type> endPoints; endPoints.resize(2*edges.size());
      for(int i = 0; i < edges.size(); i++){
	endPoints[2*i+0] = edges[i]->v1()->coordinate();
	endPoints[2*i+1] = edges[i]->v2()->coordinate();
      }

      //typedef aabb_tree<SPACE,2> edgeTree;
      //typedef typename edgeTree::node_type node_type;
      //edgeTree newTree(endPoints);
      //mAaBb =  newTree;

      //this->joinEdges();
      mesh_calculator<SPACE> calc;
      //calc.shadeVertices(*mMesh);
      //calc.shadeVertices(*mMesh);
      //23765
    }

    bool join(vector<coordinate_type> velocities, T dt){
      join_mesh<SPACE> joiner(mMesh, 0.01*minLength);
      joiner.init();
      joiner.testCollision();
      joiner.velocities = velocities;
      joiner.minLength = this->minLength;
      joiner.dt = dt;
      bool hit = joiner.join();
      int i = 0;
      bool relaxing = true;
      while(relaxing && i < 20){
	relaxing = delete_degenerate();
	i++;
      }
      return hit;
    }

    bool join(){
      join_mesh<SPACE> joiner(mMesh, 0.02*minLength);
      joiner.init();
      bool hit = joiner.join();
      int i = 0;
      bool relaxing = true;
      while(relaxing && i < 20){
	relaxing = delete_degenerate();
	i++;
      }
      return hit;
    }

    edge_ptr subdivideFace(face_vertex_ptr fv){
      TIMER functionTimer(__FUNCTION__);
      //assuming that the face vertex is the newly inserted one.
      face_vertex_ptr fv1 = fv; //this is the 
      face_vertex_ptr fv2 = fv->next();
      face_vertex_ptr fv3 = fv2->next();
      face_vertex_ptr fv4 = fv3->next();
      m2::construct<SPACE> cons;
      edge_ptr enew = cons.insert_edge(mMesh, fv1, fv3);

      return enew;
    }

    void split_edge(edge_ptr e){
      TIMER functionTimer(__FUNCTION__);
      face_vertex_ptr fv1 = e->v1();
      face_vertex_ptr fv2 = e->v2();
      coordinate_type dp = 0.5*(e->v1()->vertex()->data +
				e->v2()->vertex()->data);
      T dp2 = 0.5*(e->v1()->vertex()->data2 +
		   e->v2()->vertex()->data2);
      coordinate_type c = getButterflyWeight(e);
      subdivide<SPACE> subd;
      vertex_ptr nv = subd.subdivide_edge(mMesh,e);
      nv->coordinate() = c;
      nv->data = dp;
      nv->data2 = dp2;
      subdivideFace(fv1->next());
      subdivideFace(fv2->next());
    }

    bool insert_edges(){
      TIMER functionTimer(__FUNCTION__);
      edge_array & edges = mMesh->get_edges();
      vector<vertex_ptr>& verts = mMesh->get_vertices();

      face_array& faces = mMesh->get_faces();
      int N = faces.size();

      bool topology_change = false;
      vector<T> edgeWeights;
      vector<T> vertexWeights;
      mMesh->update_all();
      mMesh->reset_flags();
      mMesh->pack();
      mesh_calculator<SPACE> calc;

      vector<edge_ptr> edgesToSplit;
#if 1
      for (int i = 0; i < edges.size(); i++){
      	if(!mMesh->has_edge(i)) continue;
	edge_ptr ei = edges[i];
	bool pinned =
	  ei->v1()->vertex()->pinned == true &&
	  ei->v2()->vertex()->pinned == true;
	if (pinned) 	  continue;
	T l = ei->length();
	if(l > 1.75*minLength){
	  edgesToSplit.push_back(ei);
	  continue;
	}
      }
#endif

      std::cout << "sorting "<< edgesToSplit.size() << " edges" << std::endl;
      std::sort (edgesToSplit.begin(), edgesToSplit.end(), mEdgeSorter);
      std::cout << "splitting "<< edgesToSplit.size() << " edges" << std::endl;
      for(int i = edgesToSplit.size(); i > 0; i--){
	this->split_edge(edgesToSplit[i-1]);
      }

      return topology_change;
    }

    bool collapse_edges(){
      TIMER functionTimer(__FUNCTION__);
      bool topology_change = false;
      edge_array  collectedEdges;
      mMesh->reset_flags();
      vector<edge_ptr>& edges = mMesh->get_edges();
      for (int i = 0; i < edges.size(); i++){
	
	if (mMesh->has_edge(i) > 0) {
	  edge_ptr e = edges[i];
	  if(e->v1()->vertex()->flag == 1 ||
	     e->v2()->vertex()->flag == 1 
	     // e->v1()->vertex()->size() < 4 ||
	     // e->v2()->vertex()->size() < 4
	     ){
	    continue;
	  }
	  T dist = e->dist();
	  T distNext = e->v1()->next()->edge()->dist();
	  vertex_ptr v = e->v1()->vertex();
	  if (dist < minCollapseLength){
	    e->v1()->vertex()->flag = 1;
	    e->v2()->vertex()->flag = 1;
	    collectedEdges.push_back(e);
	  }
	}
      }

      std::cout << "deleting: " << collectedEdges.size() 
		<< " Tiny edges" << std::endl;

      for(int i = 0; i < collectedEdges.size(); i++){
	if(mMesh->has_edge(i)){
	  construct<SPACE> cons;
	  edge_ptr e = collectedEdges[i];
	  coordinate_type dp = 0.5*(e->v1()->vertex()->data +
				    e->v2()->vertex()->data);
	  T dp2 = 0.5*(e->v1()->vertex()->data2 +
		       e->v2()->vertex()->data2);

	  vertex_ptr vi = cons.collapse_edge(mMesh,e);
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

    bool collapse_low_curve(){
      TIMER functionTimer(__FUNCTION__);
      vector<T> edgeWeights;
      vector<T> vertexWeights;
      m2::mesh_calculator<SPACE> calc;
      calc.calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);
      bool topology_change = false;
      edge_array  collectedEdges;
      mMesh->reset_flags();
      vector<edge_ptr>& edges = mMesh->get_edges();
      for (int i = 0; i < edges.size(); i++){
	
	if (mMesh->has_edge(i) > 0) {
	  edge_ptr e = edges[i];
	  if(e->v1()->vertex()->flag == 1 ||
	     e->v2()->vertex()->flag == 1){
	    continue;
	  }
	  T dist = e->dist();
	  T distNext = e->v1()->next()->edge()->dist();
	  vertex_ptr v = e->v1()->vertex();
	  int iv0 = e->v1()->vertex()->position_in_set();
	  int iv1 = e->v2()->vertex()->position_in_set();
	  T kAvg = 0.5*(vertexWeights[iv0] + vertexWeights[iv1]);
	  if (kAvg < 200.0){
	    e->v1()->vertex()->flag = 1;
	    e->v2()->vertex()->flag = 1;
	    collectedEdges.push_back(e);
	  }
	}
      }

      std::cout << "collapsing: " << collectedEdges.size() 
		<< " low curvature edges" << std::endl;

      for(int i = 0; i < collectedEdges.size(); i++){
	if(mMesh->has_edge(i)){
	  construct<SPACE> cons;
	  edge_ptr e = collectedEdges[i];
	  coordinate_type dp = 0.5*(e->v1()->vertex()->data +
				    e->v2()->vertex()->data);
	  T dp2 = 0.5*(e->v1()->vertex()->data2 +
		       e->v2()->vertex()->data2);

	  vertex_ptr vi = cons.collapse_edge(mMesh,e);
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

    bool prune_edges(){
      TIMER functionTimer(__FUNCTION__);
      bool topology_change = false;
      //vertex_array  collectedEdges;
      // vector<T> edgeWeights;
      // vector<T> vertexWeights;
      mesh_calculator<SPACE> calc;
      // calc.calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);
      mMesh->reset_flags();
      vector<edge_ptr>& edges = mMesh->get_edges();
      vector<vertex_ptr>& verts = mMesh->get_vertices();
      vector<edge_ptr> edgesToCollapse;
      vector<T> vertThinness; vertThinness.resize(verts.size(),0.0);
      vector<coordinate_type> newCoord; newCoord.resize(verts.size());
      for (int i = 0; i < verts.size(); i++){  
	if(!verts[i]) continue;
	vertex_ptr v = verts[i];	
	T vthin = v->thinness();
	vertThinness[i] = vthin;
      }

      calc.template calcDiffuseQuantity<T>(*mMesh, vertThinness, 0.2);
      for (int i = 0; i < verts.size(); i++){  
	if(!verts[i]) continue;
	//std::cout << vertThinness[i] << std::endl;;
      }
#if 1
      for (int i = 0; i < edges.size(); i++){	
	if (mMesh->has_edge(i) > 0) {
	  edge_ptr e = edges[i];
	  T dist = e->dist();
	  T distNext = e->v1()->next()->edge()->dist();
	  vertex_ptr v1 = e->v1()->vertex();
	  vertex_ptr v2 = e->v2()->vertex();
	  if(v1->pinned) continue;
	  if(v2->pinned) continue;
	  coordinate_type n1 = e->v1()->face()->normal();
	  coordinate_type n2 = e->v2()->face()->normal();
	  // T k1 = vertexWeights[e->v1()->vertex()->position_in_set()];
	  // T k2 = vertexWeights[e->v2()->vertex()->position_in_set()];
	  // T k1p = vertexWeights[e->v1()->prev()->vertex()->position_in_set()];
	  // T k2p = vertexWeights[e->v2()->prev()->vertex()->position_in_set()];
	  //T rat = (k1p+k2p)/(k1+k2);
	  //if (rat < 2.0e-1){
	  if(0.5*(e->v1()->face()->thinness() + e->v2()->face()->thinness()) < 0.98){
	    e->v1()->vertex()->flag = 1;
	    e->v2()->vertex()->flag = 1;
	  }
	  // if (dot(n1,n2) < -0.85){
	  //   // edgesToCollapse.push_back(e->v2()->prev()->edge());

	  // }
	}

      }
#endif      

      for(int i = 0; i < verts.size(); i++){
	if(!mMesh->has_vertex(i))continue;
	if(verts[i]->size() == 1)continue;
	if(verts[i]->flag == 1){
	  topology_change = true;
	  vertex_ptr v = verts[i];
	  m2::mesh_filter<SPACE> filter;
	  newCoord[i] = 0.5*(1.0-vertThinness[i])*filter.laplacianFilterVertex(verts[i]);
	}
      }
      for(int i = 0; i < verts.size(); i++){
	if(!mMesh->has_vertex(i))continue;
	if(verts[i]->size() == 1)continue;
	if(verts[i]->flag == 1){
	  topology_change = true;
	  vertex_ptr v = verts[i];
	  m2::mesh_filter<SPACE> filter;
	  verts[i]->coordinate() += newCoord[i];
	}
      }

      return topology_change;
    }

    bool delete_degenerate(){
      TIMER functionTimer(__FUNCTION__);
      bool topology_change = false;
      vector<vertex_ptr>& verts = mMesh->get_vertices();
      vector<edge_ptr>& edges = mMesh->get_edges();
      int zeroEdges = 0;
      int twoFaces = 0;
      int flatVol = 0;
      vector<edge_ptr> edgesToDelete;
      for (int i = 0; i < edges.size(); i++){
	if (mMesh->has_edge(i)) {
	  edge_ptr e = edges[i];
	  if(e->flag == 1) continue;
	  int s1 = e->v1()->face()->size();
	  int s2 = e->v2()->face()->size();
	  T a1 = e->v1()->face()->calc_area();
	  T a2 = e->v2()->face()->calc_area();

	  if(a1 < 1e-12 || a2 < 1e-12) { 
	    e->flag = 1;
	    edgesToDelete.push_back(e);
	    continue;
	  }

	  if(s1 < 3 || s2 < 3) { 
	    twoFaces++;
	    e->flag = 1;
	    edgesToDelete.push_back(e);
	    continue;
	  }

	  if(e->v1()->face() == e->v2()->face()) { 
	    edgesToDelete.push_back(e);
	    e->flag = 1;
	    continue;
	  }

	  if(e->v1()->vertex() == e->v2()->vertex()){
	    zeroEdges++;
	    e->flag = 1;
	    edgesToDelete.push_back(e);
	    continue;
	  };

	  if(e->v1()->prev()->vertex() == e->v2()->prev()->vertex()){
	    flatVol++;
	    edgesToDelete.push_back(e);
	    e->flag = 1;
	    continue;
	  };
	}
      }
      topology_change = edgesToDelete.size() > 0;
      construct<SPACE> cons;
      for(int i = 0; i < edgesToDelete.size(); i++){
	edge_ptr e = edgesToDelete[i];
	face_ptr nf = cons.delete_edge(mMesh,e);
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
      for (int i = 0; i < verts.size(); i++){
	if(!mMesh->has_vertex(i)) continue;
	if(verts[i]->pinned) continue;
	vertex_ptr v = verts[i];

	if(verts[i]->size() < 3){
	  construct<SPACE> cons;
	  face_ptr nf = cons.delete_vertex_primitive(mMesh,v);
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
	//rem.triangulate(mMesh);
	mMesh->pack();
	mMesh->update_all();
	mMesh->reset_flags();
      }      
      return topology_change;
    }

    bool checkCompatibility(edge_ptr ei, edge_ptr ej){
      int vertIds[4];
      vertIds[0] = ei->v1()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[2] = ei->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

      vertIds[1] = ei->v2()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[5] = ei->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

      vertIds[2] = ej->v1()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[8] = ej->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

      vertIds[3] = ej->v2()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[11] = ej->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

      for(int i = 0; i < 4; i++)
	for(int j = 0; j < 4; j++){
	  if(i != j)
	    if (vertIds[i] == vertIds[j]) return false;
	}

      return true;
    };

    //---------------------------------------------------------------------------
    //sorting operators
    //---------------------------------------------------------------------------    

    struct edge_sort {
      bool operator() (edge_ptr ei, edge_ptr ej) { return (ei->length() < ej->length());}
    } mEdgeSorter;


    vector<T> calcAvgSdf(vector<coordinate_type> & coordinates,
			 vector<coordinate_type> & norms){
      vector<T> sdf; 
      vector<T> sdfWeights; 

      sdf.resize(coordinates.size(), 0.0);
      sdfWeights.resize(coordinates.size(), 0.0);

      vector<face_ptr> & faces = mMesh->get_faces();
      vector<coordinate_type> chargeCenters;
      vector<T> charges;
      vector<T> chargeMags;
      vector<T> chargeK;
      vector<coordinate_type> chargeNorms;
      for(int i = 0; i < faces.size(); i++){
	T area = faces[i]->area();
	coordinate_type norm = faces[i]->normal();
	if(area < 1e-9) continue;
	T charge =  area;

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
      std::stack<int> stack; stack.push(0);

      while(stack.size() > 0){
	int pId = stack.top(); stack.pop();
	Node & pNode = octree.nodes[pId];
	coordinate_type chargeCenter(0,0,0);
	T chargeNet = 0;
	T netChargeMag = 0;
	coordinate_type normalNet(0,0,0);
	int N = pNode.size;
	T iN = 1.0/(T)N;
	int beg = pNode.begin;
	for(int i = beg; i < beg+N;i++){
	  int ii = octree.permutation[i];
	  T chargeMag = chargeMags[ii];
	  //T chargeMag = 1.0;
	  coordinate_type chargeLoc = chargeCenters[ii];
	  coordinate_type chargeNorm = chargeNorms[ii];
	  T charge = charges[ii];
	  chargeNet    += charges[ii];
	  netChargeMag += chargeMag;
	  chargeCenter += chargeMag*chargeLoc;
	}
	chargeCenter /= netChargeMag;
	T netDistMag = 0;
	for(int i = beg; i < beg+N;i++){
	  int ii = octree.permutation[i];
	  T chargeMag = chargeMags[ii];
	  coordinate_type chargeLoc = chargeCenters[ii];
	  coordinate_type chargeNorm = chargeNorms[ii];
	  coordinate_type dp = chargeCenters[ii] - chargeCenter;
	  T dist = norm(dp);
	  T regLength = 1.5*minLength;
	  T w = 1.0/powf(dist*dist + regLength*regLength,1.0);
	  T charge = charges[ii];
	  netDistMag   += w*chargeMag;
	  normalNet    += w*chargeMag*chargeNorm;
	}

	normalNet    /= netDistMag;
	T netWeight = 0;

	normalNet.normalize();
	nodeCharges[pId]       = chargeNet; 
	nodeChargeCenters[pId] = chargeCenter; 
	nodeChargeNorms[pId]   = normalNet; 

	for(int j = 0; j < 8; j++){
	  if(pNode.children[j] != -1) stack.push(pNode.children[j]);
	}

      }
      T thresh = 0.5;

      //#pragma omp parallel for
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i = 0; i < coordinates.size(); i++){
	int count = 0;
	//if(!verts[i]) continue;
	coordinate_type pi = coordinates[i];
	std::stack<int> stack1; stack1.push(0);	

	double sdfi = sdf[i];
	double sdfWeightsi = sdfWeights[i];
	while(stack1.size() > 0){

	  int pId = stack1.top(); stack1.pop();
	  Node & pNode = octree.nodes[pId];
	  coordinate_type pj = nodeChargeCenters[pId];
	  coordinate_type Nj = nodeChargeNorms[pId];
	  Nj.normalize();
	  Debugger& debug = Debugger::get_instance();
	  //debug.DebugLines.push_back(pj);
	  //debug.DebugLines.push_back(pj+0.025*Nj);

	  coordinate_type dp = pi-pj;
	  T dc = norm(dp);
	  T sc = norm(pNode.half);

	  if(sc/dc > thresh || pNode.level < 2){
	    for(int j = 0; j < 8; j++){
	      if(pNode.children[j] != -1) {
		stack1.push(pNode.children[j]);
	      }
	    }
	  }

	  else{

	    int ii = octree.permutation[pNode.begin];
	    T ci = nodeCharges[pId];

	    //coordinate_type pj = pNode.centerOfMass;
	    if(norm(dp) < 1e-12) continue;

	    //Debugger& debug = Debugger::get_instance();
	    //debug.DebugBoxes.push_back(pNode.center);
	    //debug.DebugBoxes.push_back(pNode.half);
	    //debug.DebugLines.push_back(pi);
	    //debug.DebugLines.push_back(pj);

	    T dist = norm(dp);
	    T Px = dot(dp,Nj);

	    T i4pi = 0.25/M_PI;
	    //T regLength = 2.0*minLength;
	    //T expon = -(dist*dist)/(regLength*regLength);
	    //T w = exp(expon);
	    T regLength = 0.5*minLength;
	    T w = 1.0/powf(dist*dist + regLength*regLength,1.0);
	    //T w = 1.0/powf(dist*dist + regLength*regLength,2.0);
	    /*
	      std::cout << dp << " " << Nj 
	      << " Px: " << Px << " w: " << w 
	      << " " << dist << std::endl;	    
	      std::cout << i << " " << sdfWeightsi << " " << sdfi << std::endl;
	    */
	    sdfWeightsi += w;
	    if(pNode.size < 512)
	      sdfi += w*Px;
	    else
	      sdfi += w*dist;
	  }
	}
	sdfWeights[i] = sdfWeightsi;
	sdf[i] = sdfi;
      }

      for(int i = 0; i < sdf.size(); i++){
	//std::cout << " sdf: " << sdf[i] << " " << sdfWeights[i];
	sdf[i] *= 1.0/sdfWeights[i];
      }
      return sdf; 
    }

    
    vector<coordinate_type> integrateChargesBarnesHut(vector<coordinate_type> & coordinates){
      vector<coordinate_type> u; 
      u.resize(coordinates.size(), coordinate_type(0,0,0));

      vector<face_ptr> & faces = mMesh->get_faces();
      vector<coordinate_type> chargeCenters;
      vector<T> charges;
      vector<T> chargeMags;
      for(int i = 0; i < faces.size(); i++){
	T area = faces[i]->area();
	coordinate_type norm = faces[i]->normal();
	if(area < 1e-9) continue;
	T charge =  area;
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
      std::stack<int> stack; stack.push(0);

      while(stack.size() > 0){
	int pId = stack.top(); stack.pop();
	Node & pNode = octree.nodes[pId];
	coordinate_type chargeCenter(0,0,0);
	T chargeNet = 0;
	T netChargeMag = 0;
	int N = pNode.size;
	T iN = 1.0/(T)N;
	int beg = pNode.begin;
	for(int i = beg; i < beg+N;i++){
	  int ii = octree.permutation[i];
	  T chargeMag = chargeMags[ii];
	  //T chargeMag = 1.0;
	  coordinate_type chargeLoc = chargeCenters[ii];
	  T charge = charges[ii];
	  chargeNet    += charges[ii];
	  netChargeMag += chargeMag;
	  chargeCenter += chargeMag*chargeLoc;
	}
	chargeCenter /= netChargeMag;
	T netWeight = 0;

	nodeCharges[pId]       = chargeNet; 
	nodeChargeCenters[pId] = chargeCenter; 
	for(int j = 0; j < 8; j++){
	  if(pNode.children[j] != -1) stack.push(pNode.children[j]);
	}
      }
      T thresh = 0.5;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i = 0; i < coordinates.size(); i++){
	int count = 0;
	//if(!verts[i]) continue;
	coordinate_type pi = coordinates[i];

	std::stack<int> stack1; stack1.push(0);	
	while(stack1.size() > 0){

	  int pId = stack1.top(); stack1.pop();
	  Node & pNode = octree.nodes[pId];
	  coordinate_type pj = nodeChargeCenters[pId];
	  coordinate_type dp = pi-pj;
	  T dc = norm(dp);
	  T sc = norm(pNode.half);

	  if(sc/dc > thresh){
	    for(int j = 0; j < 8; j++){
	      if(pNode.children[j] != -1) {
		stack1.push(pNode.children[j]);
	      }
	    }
	  }

	  else{

	    int ii = octree.permutation[pNode.begin];
	    T ci = nodeCharges[pId];

	    //coordinate_type pj = pNode.centerOfMass;
	    if(norm(dp) < 1e-12) continue;
	    T dist = norm(dp);
	    dp.normalize();
	    T i4pi = 0.25/M_PI;
	    T regLength = 4.0*minLength;
	    T w = 1.0/powf(dist*dist + regLength*regLength,1.0);
	    coordinate_type ui = i4pi*ci*w*dp;
	    //std::cout << ui << std::endl;
	    u[i] += ui;
	  }
	}

      }
      return u; 
    }


    vector<coordinate_type> integrateOrientedChargesBarnesHut(vector<coordinate_type> & coordinates){
      vector<coordinate_type> u = integrateChargesBarnesHut(coordinates);
      vector<coordinate_type> N = integrateNormalsBarnesHut(coordinates);
      for(int i = 0; i < u.size(); i++){
      	T dotU = dot(u[i],N[i]);
      	if(dotU < 0) u[i] *= -1.0;
      	//u[i] *= 1.0/uSum[i];
      	//u[i].normalize();
      }
      return u;
    }

    vector<coordinate_type> integrateNormalsBarnesHut(const vector<coordinate_type> & coordinates){
      mMesh->update_all();      
      vector<face_ptr> & faces = mMesh->get_faces();

      vector<coordinate_type> chargeCenters;
      vector<coordinate_type> charges;
      vector<T> chargeMags;
      for(int i = 0; i < faces.size(); i++){
	T area = faces[i]->area();

	coordinate_type norm = faces[i]->normal();
	norm.normalize();
	if(area < 1e-9) continue;
	coordinate_type charge =  area*norm;
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
      std::stack<int> stack; stack.push(0);
      T netWeight = 0;
      while(stack.size() > 0){
	int pId = stack.top(); stack.pop();
	Node & pNode = octree.nodes[pId];
	coordinate_type chargeCenter(0,0,0);
	coordinate_type chargeNet = 0;
	T netChargeMag = 0;
	int N = pNode.size;
	T iN = 1.0/(T)N;
	int beg = pNode.begin;
	for(int i = beg; i < beg+N;i++){
	  int ii = octree.permutation[i];
	  T chargeMag = chargeMags[ii];
	  //T chargeMag = 1.0;
	  coordinate_type chargeLoc = chargeCenters[ii];
	  coordinate_type charge = charges[ii];
	  netChargeMag += chargeMag;
	  chargeCenter += chargeMag*chargeLoc;
	  T w = 1.0;
	  chargeNet += w*charges[ii];

	  netWeight += w;
	}
	chargeCenter /= netChargeMag;


	nodeCharges[pId]       = chargeNet; 
	nodeChargeCenters[pId] = chargeCenter; 
	for(int j = 0; j < 8; j++){
	  if(pNode.children[j] != -1) stack.push(pNode.children[j]);
	}
      }
      T thresh = 0.5;
      vector<coordinate_type> u; 
      vector<T> uSum; 
      u.resize(coordinates.size(), coordinate_type(0,0,0));
      uSum.resize(coordinates.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i = 0; i < coordinates.size(); i++){
	int count = 0;
	//if(!verts[i]) continue;
	coordinate_type pi = coordinates[i];

	std::stack<int> stack1; stack1.push(0);	
	while(stack1.size() > 0){

	  int pId = stack1.top(); stack1.pop();
	  Node & pNode = octree.nodes[pId];
	  coordinate_type pj = nodeChargeCenters[pId];
	  coordinate_type dp = pi-pj;
	  T dc = norm(dp);
	  T sc = norm(pNode.half);

	  if(sc/dc > thresh){
	    for(int j = 0; j < 8; j++){
	      if(pNode.children[j] != -1) {
		stack1.push(pNode.children[j]);
	      }
	    }
	  }

	  else{

	    int ii = octree.permutation[pNode.begin];
	    coordinate_type ci = nodeCharges[pId];
	    if(norm(dp) < 1e-16) continue;

	    if(norm(ci) < 1e-16) continue;
	    ci.normalize();
	    //coordinate_type pj = pNode.centerOfMass;
	    //coordinate_type ndp = dp; ndp.normalize();
	    //if(dot(ci,ndp) < 0.0) continue;
	    //coordinate_type ciN = ci; ciN.normalize();
	    //coordinate_type dpp = dot(dp,ciN)/dot(dp,dp)*dp;
	    T dist = norm(dp);

	    T i4pi = 0.25/M_PI;
	    T regLength = 1.0*minLength;
	    T m = 2.0;
	    T w = 1.0/powf(dist*dist + regLength*regLength,m);
	    //T w = exp(-dist*dist/regLength*regLength);
	    u[i] += w*ci;
	    //u[i] += w*coordinate_type(0.01,0.01,0.01);
	    uSum[i] += w;
	  }
	}
      }

      for(int i = 0; i < u.size(); i++){
	if(uSum[i] < 1e-9) continue;
	u[i] *= 1.0/uSum[i];
	u[i].normalize();
      }
      return u; 
    }

    vector<coordinate_type> integrateNormalDerivativeBarnesHut(vector<coordinate_type> & coordinates){
      mMesh->update_all();      
      vector<face_ptr> & faces = mMesh->get_faces();

      vector<coordinate_type> chargeCenters;
      vector<coordinate_type> charges;
      vector<T> chargeMags;
      for(int i = 0; i < faces.size(); i++){
	T area = faces[i]->area();

	coordinate_type norm = faces[i]->normal();
	norm.normalize();
	//if(area < 1e-9) continue;
	coordinate_type charge =  area*norm;
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
      std::stack<int> stack; stack.push(0);

      while(stack.size() > 0){
	int pId = stack.top(); stack.pop();
	Node & pNode = octree.nodes[pId];
	coordinate_type chargeCenter(0,0,0);
	coordinate_type chargeNet = 0;
	T netChargeMag = 0;
	int N = pNode.size;
	T iN = 1.0/(T)N;
	int beg = pNode.begin;
	for(int i = beg; i < beg+N;i++){
	  int ii = octree.permutation[i];
	  T chargeMag = chargeMags[ii];
	  //T chargeMag = 1.0;
	  coordinate_type chargeLoc = chargeCenters[ii];
	  coordinate_type charge = charges[ii];
	  netChargeMag += chargeMag;
	  chargeCenter += chargeMag*chargeLoc;
	  chargeNet += charges[ii];
	}
	chargeCenter /= netChargeMag;
	nodeCharges[pId]       = chargeNet; 
	nodeChargeCenters[pId] = chargeCenter; 
	for(int j = 0; j < 8; j++){
	  if(pNode.children[j] != -1) stack.push(pNode.children[j]);
	}
      }
      T thresh = 0.5;
      vector<coordinate_type> u; 
      vector<T> uSum; 
      u.resize(coordinates.size(), coordinate_type(0,0,0));
      uSum.resize(coordinates.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for
#endif 
      for(int i = 0; i < coordinates.size(); i++){
	int count = 0;
	//if(!verts[i]) continue;
	coordinate_type pi = coordinates[i];

	std::stack<int> stack1; stack1.push(0);	
	while(stack1.size() > 0){

	  int pId = stack1.top(); stack1.pop();
	  Node & pNode = octree.nodes[pId];
	  coordinate_type pij = nodeChargeCenters[pId];
	  coordinate_type dpi = pi-pij;
	  T dc = norm(dpi);
	  T sc = norm(pNode.half);

	  if(sc/dc > 0.5){
	    for(int j = 0; j < 8; j++){
	      if(pNode.children[j] != -1) {
		stack1.push(pNode.children[j]);
	      }
	    }
	  }

	  else{

	    int ii = octree.permutation[pNode.begin];
	    coordinate_type ci = nodeCharges[pId];

	    //coordinate_type pj = pNode.centerOfMass;
	    if(norm(dpi) < 1e-16) continue;
	    //coordinate_type ciN = ci; ciN.normalize();
	    //coordinate_type dpp = dot(dp,ciN)/dot(dp,dp)*dp;
	    T dist = norm(dpi);

	    T i4pi = 0.25/M_PI;
	    T regLength = 1.0*minLength;
	    T m = 2.0;
	    T wi = 1.0/powf(dist*dist + regLength*regLength,m);

	    int count = 0;
	    //if(!verts[i]) continue;
	    coordinate_type pi = coordinates[i];

	    std::stack<int> stack2; stack2.push(0);	
	    while(stack2.size() > 0){

	      int pJd = stack2.top(); stack2.pop();
	      Node & pNodej = octree.nodes[pJd];
	      coordinate_type pjj = nodeChargeCenters[pJd];
	      coordinate_type dpj = pi-pjj;
	      T dcj = norm(dpj);
	      T scj = norm(pNodej.half);

	      if(scj/dcj > 0.5){
		for(int j = 0; j < 8; j++){
		  if(pNodej.children[j] != -1) {
		    stack2.push(pNodej.children[j]);
		  }
		}
	      }

	      else{
		T distj = norm(dpj);
		  
		T i4pi = 0.25/M_PI;
		T m = 2.0;
		T wj = 1.0/powf(distj*distj + regLength*regLength,m);
		u[i] += wi*wj*ci;

	      }
	    }
	    uSum[i] += wi;

	  }
	}
      }
      for(int i = 0; i < u.size(); i++){
	u[i] *= 1.0/uSum[i]/uSum[i];
      }
      return u; 
    }

    void integrateNormalsRK2(T dt){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      vector<coordinate_type> y0;
      mMesh->update_all(); mMesh->reset_flags(); mMesh->pack();
      y0.resize(verts.size());
      vector<T> edgeWeights;
      vector<T> vertexWeights;
      for(int i = 0; i < verts.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;
	y0[i] = verts[i]->coordinate();
      }

      vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
      T norm2 = 0;
      vector<coordinate_type> y1(y0);
      for(int i = 0; i < u0.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;

	y1[i] += 0.5*dt*u0[i];
	norm2 += dot(u0[i],u0[i]);
      }
      std::cout <<" first stage norm: "<< sqrt(norm2) << std::endl;
      vector<coordinate_type> u1 = integrateNormalsBarnesHut(y1);
      norm2 = 0;
      for(int i = 0; i < verts.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;	
	verts[i]->coordinate() += dt*u1[i];
	norm2 += dot(u1[i],u1[i]);
      }
      std::cout <<" second stage norm: "<< sqrt(norm2) << std::endl;
    }

    vector<coordinate_type>  integrateNormals(T dt, vector<coordinate_type> & coordinates){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      vector<coordinate_type> y0(coordinates);
      mMesh->update_all(); mMesh->reset_flags(); mMesh->pack();

      vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
      T norm2 = 0;
      for(int i = 0; i < u0.size(); i++){
	y0[i] += dt*u0[i];
	norm2 += dot(u0[i],u0[i]);
      }
      std::cout <<" first stage norm: "<< sqrt(norm2) << std::endl;
      return y0;
    }

    vector<coordinate_type>  integrateNormals(vector<T> dt, vector<coordinate_type> & coordinates){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      vector<coordinate_type> y0(coordinates);
      mMesh->update_all(); mMesh->reset_flags(); mMesh->pack();

      vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
      T norm2 = 0;
      for(int i = 0; i < u0.size(); i++){
	y0[i] += dt[i]*u0[i];
	norm2 += dot(u0[i],u0[i]);
      }
      std::cout <<" first stage norm: "<< sqrt(norm2) << std::endl;
      return y0;
    }

    vector<coordinate_type>  integrateCharges(T dt, vector<coordinate_type> & coordinates){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      vector<coordinate_type> y0(coordinates);
      mMesh->update_all(); mMesh->reset_flags(); mMesh->pack();

      vector<coordinate_type> u0 = integrateOrientedChargesBarnesHut(y0);
      T norm2 = 0;
      for(int i = 0; i < u0.size(); i++){
	y0[i] += dt*u0[i];
	norm2 += dot(u0[i],u0[i]);
      }
      std::cout <<" first stage norm: "<< sqrt(norm2) << std::endl;
      return y0;
    }

    vector<coordinate_type>  integrateMorphing(T dt,
					       vector<coordinate_type> & coordinates,
					       vector<coordinate_type> & normals){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      vector<coordinate_type>  n0 = mMesh->get_normals();

      vector<coordinate_type> y0(coordinates);
      mMesh->update_all(); mMesh->reset_flags(); mMesh->pack();

      vector<coordinate_type> u0 = integrateMorphingCharges(y0, normals);
      //vector<coordinate_type> u1 = integrateChargesBarnesHut(y0);
      T norm2 = 0;
      for(int i = 0; i < u0.size(); i++){
	y0[i] +=  dt*u0[i];
	norm2 += dot(u0[i],u0[i]);
      }
      std::cout <<" first stage norm: "<< sqrt(norm2) << std::endl;
      return y0;
    }

    vector<coordinate_type> integrateNormalsRK2(T dt, vector<coordinate_type> & coordinates){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      vector<coordinate_type> y0(coordinates);
      vector<coordinate_type> y1(y0);
      mMesh->update_all(); mMesh->reset_flags(); mMesh->pack();

      vector<coordinate_type> u0 = integrateNormalsBarnesHut(y0);
      T norm2 = 0;
      for(int i = 0; i < u0.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;

	y0[i] += 0.5*dt*u0[i];
	norm2 += dot(u0[i],u0[i]);
      }
      std::cout <<" first stage norm: "<< sqrt(norm2) << std::endl;
      vector<coordinate_type> u1 = integrateNormalsBarnesHut(y0);
      norm2 = 0;
      for(int i = 0; i < coordinates.size(); i++){
	y1[i] += dt*u1[i];
	norm2 += dot(u1[i],u1[i]);
      }
      std::cout <<" second stage norm: "<< sqrt(norm2) << std::endl;
      return y1;
    }


    void hashFaces(){
      face_bin<SPACE> newHash(*mMesh, 2.0*minLength);
      std::cout << "hashed faces: " << newHash.faces().size() << std::endl;
      mFaceHash = newHash;
    }

    // void staticAddHandles(){
    //   vector<vertex_ptr>& verts = mMesh->get_vertices();
    //   for (int i = 0; i < verts.size(); i++){
    //   }
    // }

    void randomize(){
      vertex_array& vl = mMesh->get_vertices();
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
    //aabb_tree<SPACE,2> & aabb(){return mAaBb;}
    face_bin<SPACE> & faceHash(){return mFaceHash;}
    face_bin<SPACE> mFaceHash;
    //aabb_tree<SPACE,2> mAaBb;
    vector<coordinate_type> mDirectionWeights;
    vector<coordinate_type> mDebug;
    vector<coordinate_type*> mDirectionField;
    T maxCurvature, minCurvature, minLength, minCollapseLength, maxLength;
    control_ptr		 mMesh;

  };

  template <typename SPACE>
  class vortex_sheet{
    M2_TYPEDEFS;
  public:
    vortex_sheet(control_ptr obj_in){

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
      coordinate_type lengths = max-min;
      int minRes = 64;
      T minLength = lengths[0];
      minLength = minLength < lengths[1] ? minLength : lengths[1];
      minLength = minLength < lengths[2] ? minLength : lengths[2];

      T dx = minLength/(T)(minRes-1);
      // while(relaxing && k < 4){
      // 	set_operations<SPACE> setOps;
      // 	relaxing = setOps.flip_edges(mMesh);
      // 	k++;
      // }
      for(int i = 0; i < 1; i++){
	this->remesh();
       	set_operations<SPACE> setOps;
	bool relaxing = true;
	while(relaxing)
	  relaxing = setOps.flip_edges(mMesh);
      }
      //pin_bottom();
      //pin_half();

      this->hashFaces();
      //rem.triangulate(mMesh);
      //construct<SPACE> cons;
      //face_ptr nf = cons.delete_vertex(mMesh,mMesh->get_vertices()[0]);
      //vector<T> edgeWeights;
      //vector<T> vertexWeights;
      //calcCurveFlowNormal(*mMesh, vertexWeights, edgeWeights);

    }

    coordinate_type getButterflyVertexWeight(face_vertex_ptr v){
      // \/\./
      //vn = next 
      //vnn = vertext_next next
      //pvpp = prev vertex_prev next
      //coordinate_type c = v->coordinate();
      // coordinate_type cn = v->next()->coordinate();
      // coordinate_type cvnn = v->vnext()->next()->coordinate();
      // coordinate_type cpvpp = v->prev()->vprev()->prev()->coordinate();
      // return (8.0*c + 2.0*cn - 1.0*(cvnn + cpvpp));
      face_vertex_ptr itb = v;
      face_vertex_ptr ite = v->vprev();
      bool iterating = true;
      int K = 0, j=0;
      while(iterating){
	iterating = itb != ite;
	itb = itb->vnext();
	K++;
      }
      itb = v;
      iterating = true;
      T s = 0;
      coordinate_type cs(0,0,0);
      while(iterating){
	iterating = itb != ite;
	itb = itb->vnext();
	T frac = (T)j/(T)K;
	T sj = (0.25 + cos(2.0*M_PI*frac) + 0.5*cos(4.0*M_PI*frac))/(T)K;
	coordinate_type ci = itb->next()->coordinate();
	cs += sj*ci;
	s += sj;
	j++;
      }      
      coordinate_type c = v->coordinate();
      return 0.75*c + cs;
    } 

    coordinate_type getButterflyWeight(edge_ptr e){
      face_vertex_ptr v1 = e->v1();
      face_vertex_ptr v2 = e->v2();
      coordinate_type c1 = getButterflyVertexWeight(v1);
      coordinate_type c2 = getButterflyVertexWeight(v2);
      return 0.5*(c1 + c2);
    }


    void flip_edges(){
      int k = 0;
      bool relaxing = true;
      while(relaxing && k < 1){
       	set_operations<SPACE> setOps;
       	relaxing = setOps.flip_edges(mMesh);
       	k++;
      }
    }

    void remesh(){      
      bool relaxing = true;
      int k = 0;
      TIMER functionTimer(__FUNCTION__);
      std::cout << " mesh size: " << mMesh->get_vertices().size() << std::endl;
      // relaxing = this->delete_vertices_curve();

      mMesh->colorVertices();
      std::cout << " max graph color: " << mMesh->maxGraphColor << std::endl;

      relaxing = true;
      while(relaxing){	
       	relaxing = delete_degenerate();
      }

      relaxing = this->collapse_edges();

      relaxing = true;
      while(relaxing){
	relaxing = delete_degenerate();
      }

      relaxing = true;
      while(relaxing){
      	relaxing = this->insert_edges_curve();
      }

      relaxing = true;
      while(relaxing){
	relaxing = delete_degenerate();
      }
      /*
	mesh_filter<SPACE> filter;
	for(int i = 0; i < 1; i++){
      	mMesh->update_all();
      	//filter.taubinFilter(*mMesh,0.1,0.11);
      	filter.nullLaplacianFilter(*mMesh,0.25);
	}
      */
      
      //this->hashFaces();
      mesh_calculator<SPACE> calc;
      //calc.shadeVerticesWinding(*mMesh);
      calc.shadeVertices(*mMesh);
    }

    //face_bin<SPACE> & getHash(){return mFaceHash;}

    pole_tree<SPACE> & octree(){return mOctree;}
    //aabb_tree<SPACE,3> & aabb(){return mAaBb;}
    vertex_bin<SPACE> & vertexHash(){return mVertexHash;}
    edge_bin<SPACE> & edgeHash(){return mEdgeHash;}
    face_bin<SPACE> & faceHash(){return mFaceHash;}


    void hashVertices(){
      vertex_bin<SPACE> newHash(*mMesh, 1.0*regLength);
      mVertexHash = newHash;
    }

    void hashEdges(){
      edge_bin<SPACE> newHash(*mMesh, 1.0*regLength);
      mEdgeHash = newHash;
    }

    void hashFaces(){
      face_bin<SPACE> newHash(*mMesh, 2.0*minLength);
      mFaceHash = newHash;
    }

    void pin_half(){
      mMesh->pack();
      vertex_array& vl = mMesh->get_vertices();
      coordinate_type cen = mMesh->calc_center();
      coordinate_type min = mMesh->calc_min();
      coordinate_type max = mMesh->calc_max();
      for (int j = 0; j < vl.size(); j++) {
	vl[j]->pinned = false;
      }

      std::cout << " cen: " << cen << std::endl;
      for (int j = 0; j < vl.size(); j++) {
	if(vl[j]->coordinate()[1] < min[1] + 0.25*(max[1] - min[1])){
	  vl[j]->pinned = true;
	}
      }
    }

    void pin_bottom(){
      mMesh->pack();
      vertex_array& vl = mMesh->get_vertices();
      coordinate_type cen = mMesh->calc_center();
      coordinate_type min = mMesh->calc_min();
      coordinate_type max = mMesh->calc_max();
      for (int j = 0; j < vl.size(); j++) {
	vl[j]->pinned = false;
      }

      std::cout << " cen: " << cen << std::endl;
      for (int j = 0; j < vl.size(); j++) {
	if(vl[j]->coordinate()[1] < min[1] + 0.05*(max[1] - min[1])){
	  vl[j]->pinned = true;
	}
      }
    }

    void pin_bottom_facing(){
      mMesh->pack();
      vertex_array& vl = mMesh->get_vertices();
      coordinate_type cen = mMesh->calc_center();
      coordinate_type min = mMesh->calc_min();
      coordinate_type max = mMesh->calc_max();
      for (int j = 0; j < vl.size(); j++) {
	vl[j]->pinned = false;
      }

      std::cout << " cen: " << cen << std::endl;
      for (int j = 0; j < vl.size(); j++) {
	coordinate_type N = vl[j]->normal();
	T downAngle = dot(N,coordinate_type(0,-1.0,0));
	T bottomThresh = min[1] + 0.2*(max[1] - min[1]);
	if(vl[j]->coordinate()[1] < bottomThresh && downAngle > -0.5){
	  vl[j]->pinned = true;
	}
      }
    }

    void addCurveVorticity(T c, int dim){
      mMesh->update_all();
      mMesh->pack();
      vector<vertex_ptr>& verts = mMesh->get_vertices();
      vector<face_ptr>&   faces = mMesh->get_faces();
      vector<mat3> tensorField;
      vector<coordinate_type> weights;
      tensorField.resize(verts.size());
      weights.resize(verts.size());
      mesh_calculator<SPACE> calc;
      calc.calculateBiDirectionField(*mMesh, tensorField, weights, 2.0*minLength);
      for(int i = 0; i < faces.size(); i++){
	int i1 = faces[i]->fbegin()->vertex()->position_in_set();
	int i2 = faces[i]->fbegin()->next()->vertex()->position_in_set();
	int i3 = faces[i]->fbegin()->next()->next()->vertex()->position_in_set();
	
	coordinate_type v1(tensorField[i1](dim,0),tensorField[i1](dim,1),tensorField[i1](dim,2)); 
	coordinate_type v2(tensorField[i2](dim,0),tensorField[i2](dim,1),tensorField[i2](dim,2)); 
	coordinate_type v3(tensorField[i3](dim,0),tensorField[i3](dim,1),tensorField[i3](dim,2)); 

	T wa = 1.0/3.0*(weights[i1][dim] + weights[i2][dim] + weights[i3][dim]);
	coordinate_type va = 1.0/3.0*(v1 + v2 + v3); 
	coordinate_type g = coordinate_type(0,1,0); //make this globally consistent
	T s = dot(g,va)/dot(va,va);
	coordinate_type gonw = s*wa; gonw.normalize();
	faces[i]->data += c*gonw*wa;
      }
      
    }
   
    void updateVorticity(){
      T normd = 0;
      T maxMag = 0;
      vector<face_ptr> & faces = mMesh->get_faces();
      
      for(int i = 0; i < faces.size(); i++){
	if(!faces[i]) continue;
	circulationToVorticity(faces[i]);
	normd += dot(faces[i]->data,faces[i]->data);
	T mag = norm(faces[i]->data);
	maxMag = maxMag > mag ? maxMag:mag;
      }
      std::cout <<" vorticity update norm: "
		<< sqrt(normd) << " max magnitude: " << maxMag << std::endl;
    }

    void updateCirculation(){
      TIMER functionTimer(__FUNCTION__);
      T norm = 0;
      vector<face_ptr> & faces = mMesh->get_faces();
#ifdef _OPENMP
#pragma omp parallel for
#endif 
      for(int i = 0; i < faces.size(); i++){
	if(!faces[i]) continue;
	vorticityToCirculation(faces[i]);
	norm += dot(faces[i]->data,faces[i]->data);
      }
      std::cout <<" circulation update norm: "<< sqrt(norm) << std::endl;
    }    

    void integrateBaroclinity(vector<coordinate_type> g,
			      vector<T> s, T dt ){
      TIMER functionTimer(__FUNCTION__);
      mMesh->update_all();
      vector<face_ptr> & faces = mMesh->get_faces();
      T norm = 0;
      for(int i = 0; i < faces.size(); i++){
	if(!faces[i]) continue;
	if(faces[i]->fbegin()->vertex()->pinned) continue;
	if(faces[i]->fbegin()->next()->vertex()->pinned) continue;
	if(faces[i]->fbegin()->next()->next()->vertex()->pinned) continue;
	T a = faces[i]->calc_area();
	if(a > 1e-16){
	  //T Ba  = 2.0*faces[i]->data3*0.0005;
	  T Ba  = 2.0*10.0;
	  faces[i]->update_normal();
	  coordinate_type normal = faces[i]->normal();
	  coordinate_type dB = -dt*Ba*s[i]*cross(g[i],normal);
	  faces[i]->data += dB;
	  faces[i]->data2 = dB;
	  norm += dot(faces[i]->data,faces[i]->data);
	}
      }
      std::cout <<" vorticity norm: "<< sqrt(norm) << std::endl;
    }

    void integrateBaroclinity(T dt){
      TIMER functionTimer(__FUNCTION__);
      mMesh->update_all();
      vector<face_ptr> & faces = mMesh->get_faces();
      T norm = 0;
      for(int i = 0; i < faces.size(); i++){
	if(!faces[i]) continue;
	if(faces[i]->fbegin()->vertex()->pinned) continue;
	if(faces[i]->fbegin()->next()->vertex()->pinned) continue;
	if(faces[i]->fbegin()->next()->next()->vertex()->pinned) continue;
	T a = faces[i]->calc_area();
	if(a > 1e-16){
	  //T Ba  = 2.0*faces[i]->data3*0.0005;
	  T Ba  = 2.0*0.005;
	  faces[i]->update_normal();
	  coordinate_type normal = faces[i]->normal();
	  coordinate_type dB = -dt*Ba*cross(coordinate_type(0,-9.8,0),normal);
	  faces[i]->data += dB;
	  faces[i]->data2 = dB;
	  norm += dot(faces[i]->data,faces[i]->data);
	}
      }
      std::cout <<" vorticity norm: "<< sqrt(norm) << std::endl;
    }


    void integrateBaroclinity(T dt,T C){
      TIMER functionTimer(__FUNCTION__);
      mMesh->update_all();
      vector<face_ptr> & faces = mMesh->get_faces();
      T norm = 0;
      for(int i = 0; i < faces.size(); i++){
	if(!faces[i]) continue;
	if(faces[i]->fbegin()->vertex()->pinned) continue;
	if(faces[i]->fbegin()->next()->vertex()->pinned) continue;
	if(faces[i]->fbegin()->next()->next()->vertex()->pinned) continue;
	T a = faces[i]->calc_area();
	if(a > 1e-16){
	  //T Ba  = 2.0*faces[i]->data3*0.0005;
	  T Ba  = 2.0*C;
	  faces[i]->update_normal();
	  coordinate_type normal = faces[i]->normal();
	  coordinate_type dB = -dt*Ba*cross(coordinate_type(0,-9.8,0),normal);
	  faces[i]->data += dB;
	  faces[i]->data2 = dB;
	  norm += dot(faces[i]->data,faces[i]->data);
	}
      }
      std::cout <<" vorticity norm: "<< sqrt(norm) << std::endl;
    }

    void integrateVelocityEuler(T dt){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      //vector<coordinate_type> u0 = integrateVelocity(0.5*dt);
      vector<coordinate_type> u0 = integrateVelocityBarnesHut();
      //vector<coordinate_type> u0 = integrateVelocityBruteForce(0.5*dt);
      for(int i = 0; i < verts.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;
	verts[i]->coordinate() += dt*u0[i];
	verts[i]->data = u0[i];
      }
    }

    void integrateVelocityRK2(T dt){
      TIMER functionTimer(__FUNCTION__);
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      vector<coordinate_type> y0;
      y0.resize(verts.size());
      mMesh->update_all(); mMesh->reset_flags(); mMesh->pack();

      for(int i = 0; i < verts.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;
	y0[i] = verts[i]->coordinate();
      }
      //vector<coordinate_type> uheat = integrateHeat(0.5*dt);

      //vector<coordinate_type> u0 = integrateVelocity(0.5*dt);
      vector<coordinate_type> u0 = integrateVelocityBarnesHut();
      //vector<coordinate_type> u0 = integrateVelocityBruteForce(0.5*dt);
      T norm2 = 0;
      for(int i = 0; i < verts.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;

	verts[i]->coordinate() += 0.5*dt*(u0[i]);
	//verts[i]->data = u0[i] + uheat[i];
	norm2 += dot(u0[i],u0[i]);
      }
      std::cout <<" first stage norm: "<< sqrt(norm2) << std::endl;
      norm2 = 0;
      //vector<coordinate_type> u1 = integrateVelocity(0.5*dt);
      vector<coordinate_type> u1 = integrateVelocityBarnesHut();
      //vector<coordinate_type> u1 = integrateVelocityBruteForce(0.5*dt);
      T maxMag = 0;
      for(int i = 0; i < verts.size(); i++){
	if(!verts[i]) continue;
	if(verts[i]->pinned) continue;	

	coordinate_type du = dt*(u1[i]);
	verts[i]->coordinate() = y0[i] + du;
	//verts[i]->data = u1[i] + uheat[i];
	norm2 += dot(u1[i],u1[i]);
	T mag = norm(u1[i]);
	maxMag = maxMag > mag ? maxMag : mag;
      }
      std::cout << " second stage velocity norm: " << sqrt(norm2)
		<< " max velocity mag: " << maxMag <<  std::endl;
    }

    vector<coordinate_type> integrateHeat(T dt){
      vector<coordinate_type> u; 
      u.resize(mMesh->get_vertices().size(), coordinate_type(0,0,0));
      vector<vertex_ptr>& tverts = mMesh->get_vertices();
      T heatC = 0.05;
      for (long i = 0; i < tverts.size(); i++) {

	if(!mMesh->has_vertex(i)) continue;
	vertex_ptr v = tverts[i];
	if(v->size() ==  0) continue;

	face_vertex_ptr itb = v->fbegin();
	face_vertex_ptr ite = v->fend();
	bool at_head = false;
	T cumulativeDens = 0;
	int k = 0;
	while (!at_head) {
	  at_head = itb==ite;
	  mesh_calculator<SPACE> calc;
	  cumulativeDens += itb->face()->data3;
	  itb = itb->vnext();
	  k++;
	}

	T avgDens = cumulativeDens/(T)k;
	if(avgDens > 1e-10)
	  u[i] = coordinate_type(0.0,dt*heatC*avgDens,0.0);
	else
	  u[i] = coordinate_type(0,0,0);
      }
      return u;
    }

    vector<coordinate_type> integrateVelocityBruteForce(T dt){
      vector<coordinate_type> u; 
      u.resize(mMesh->get_vertices().size(), coordinate_type(0,0,0));

      vector<face_ptr> & faces = mMesh->get_faces();
      vector<coordinate_type> chargeCenters;
      vector<coordinate_type> charges;
      vector<T> chargeMags;

      for(int i = 0; i < faces.size(); i++){
	T area = faces[i]->area();
	coordinate_type vort = faces[i]->data;
	coordinate_type charge =  area*vort;
	chargeCenters.push_back(faces[i]->center());
	charges.push_back(charge);
	chargeMags.push_back(norm(charge));
      }

      vector<vertex_ptr> & verts = mMesh->get_vertices();
      for(int i = 0; i < verts.size(); i++){
	for(int j = 0; j < charges.size(); j++){
	  if(!verts[i]) continue;
	  coordinate_type pi = verts[i]->coordinate();
	  coordinate_type pj = chargeCenters[j];

	  coordinate_type ci = charges[j];
	  //coordinate_type pi = chargeCenters[ii];
	  coordinate_type dp = pi-pj;
	  if(norm(dp) < 1e-12) continue;
	  T dist = norm(dp);
	  T i4pi = 0.25/M_PI;

	  // T dist3 = dist*dist*dist;
	  // T l3   = regLength*regLength*regLength;	     
	  // T kappa = (1.0-exp(-dist3/l3))/dist3;
	  
	  T l2 = regLength*regLength;
	  T denom = powf((dist*dist + l2),1.5);
	  coordinate_type fR = dp/denom;
	  //coordinate_type fR = dp*kappa;
	  coordinate_type ui = i4pi*cross(ci,fR);
	  //std::cout << ui << std::endl;
	  u[i] += ui;
	}
      }
      return u;
    }

    vector<coordinate_type> integrateVelocityBarnesHut(){
      vector<coordinate_type> u; 
      u.resize(mMesh->get_vertices().size(), coordinate_type(0,0,0));

      vector<face_ptr> & faces = mMesh->get_faces();
      vector<coordinate_type> chargeCenters;
      vector<coordinate_type> charges;
      vector<T> chargeMags;
      TIMER buildTimer("building tree hierarchy");
      for(int i = 0; i < faces.size(); i++){
	T area = faces[i]->area();
	coordinate_type vort = faces[i]->data;
	coordinate_type charge =  area*vort;
	chargeCenters.push_back(faces[i]->center());
	charges.push_back(charge);
	chargeMags.push_back(norm(charge));
      }

      typedef pole_tree<SPACE> Tree;
      typedef pole_node<SPACE> Node;
      Tree octree(chargeCenters);

      vector<coordinate_type> nodeCharges;
      vector<coordinate_type> nodeChargeCenters;
      nodeCharges.resize(octree.nodes.size());
      nodeChargeCenters.resize(octree.nodes.size());
      std::stack<int> stack; stack.push(0);
      T netWeight = 0;
      while(stack.size() > 0){
	int pId = stack.top(); stack.pop();
	Node & pNode = octree.nodes[pId];
	coordinate_type chargeCenter(0,0,0);
	coordinate_type chargeNet(0,0,0);
	T netChargeMag = 0;
	int N = pNode.size;
	T iN = 1.0/(T)N;
	int beg = pNode.begin;
	for(int i = beg; i < beg+N;i++){
	  int ii = octree.permutation[i];
	  T chargeMag = chargeMags[ii] + .001;
	  //T chargeMag = 1.0;
	  coordinate_type chargeLoc = chargeCenters[ii];
	  coordinate_type charge = charges[ii];
	  netChargeMag += chargeMag;
	  chargeCenter += chargeMag*chargeLoc;

	  chargeNet += charges[ii];
	}
	chargeCenter /= netChargeMag;

	nodeCharges[pId]       = chargeNet; 
	nodeChargeCenters[pId] = chargeCenter; 
	for(int j = 0; j < 8; j++){
	  if(pNode.children[j] != -1) stack.push(pNode.children[j]);
	}
      }
      T thresh = 0.5;
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      TIMER traversalTimer("traversing tree");
#ifdef _OPENMP
#pragma omp parallel for
#endif 
      for(int i = 0; i < verts.size(); i++){
	int count = 0;
	if(!verts[i]) continue;
	coordinate_type pi = verts[i]->coordinate();

	std::stack<int> stack1; stack1.push(0);	
	while(stack1.size() > 0){

	  int pId = stack1.top(); stack1.pop();
	  Node & pNode = octree.nodes[pId];
	  coordinate_type pj = nodeChargeCenters[pId];
	  coordinate_type ci = nodeCharges[pId];
	  coordinate_type dp = pi-pj;
	  T dc = norm(dp);
	  T sc = norm(pNode.half);


	  //int ii = octree.permutation[pNode.begin];


	  //coordinate_type pj = pNode.centerOfMass;
	  T dist = dc;

	  // T l2 = regLength*regLength;
	  // T denom = powf((dist*dist + l2),1.5);
	  //coordinate_type fR = dp/denom;
	  T dist3 = dist*dist*dist;
	  T l3   = regLength*regLength*regLength;	     
	  T kappa = (1.0-exp(-dist3/l3))/dist3;
	  T i4pi = 0.25/M_PI;
	  if(sc/dc > thresh && kappa*dc*norm(ci) > 1e-16){
	    //if(sc/dc > thresh){
	    for(int j = 0; j < 8; j++){
	      if(pNode.children[j] != -1) {
		stack1.push(pNode.children[j]);
	      }
	    }
	  }

	  else{
	    //std::cout << ui << std::endl;

	    coordinate_type fR = dp*kappa;
	    coordinate_type ui = i4pi*cross(ci,fR);
	    u[i] += ui;
	  }
	}

      }
      return u; 
    }

    vector<coordinate_type> integrateVelocity(T dt){
      TIMER functionTimer(__FUNCTION__);
      //integrateBaroclinity(dt);
      //second order runge kutta
      this->hashVertices();
      this->hashFaces();
      T dx = mVertexHash.dx();
      int xRes = mFaceHash.xRes();
      int yRes = mFaceHash.yRes();
      int zRes = mFaceHash.zRes();
      T i4pi = 1.0/(4.0*M_PI);
      std::cout << dx << " " << regLength << std::endl;
      int off = 3;
      int ipN = 2*off + 1;
      int sz = 0;
      T cutOff = 3.5;
      for(int k = -off; k < off + 1; k++)
	for(int j = -off; j < off + 1; j++)
	  for(int i = -off; i < off + 1; i++){
	    T rad = norm(coordinate_type(i*dx,j*dx,k*dx));
	    if(rad <= cutOff*regLength)
	      sz++;
	  }
      std::cout << "query sphere size: "<< sz << std::endl;
      int * ip = new int[3*sz];
      int ind = 0;
      for(int k = -off; k < off + 1; k++)
	for(int j = -off; j < off + 1; j++)
	  for(int i = -off; i < off + 1; i++){
	    T rad = norm(coordinate_type(i*dx,j*dx,k*dx));
	    if(rad <= cutOff*regLength){
	      int index = (i + 1) + (j+1)*ipN + (k+1)*ipN*ipN;
	      ip[3*ind + 0] = i;
	      ip[3*ind + 1] = j;
	      ip[3*ind + 2] = k;
	      ind++;
	    }
	  }
      
      T imax = dx/xRes;

      std::cout << " integrating: " << mMesh->get_vertices().size() << std::endl;
      vector<coordinate_type> u; 
      u.resize(mMesh->get_vertices().size(), coordinate_type(0,0,0));
      int blockAvg = 0;
      int blockCount = 0;
      vector<int> &  binStart    = mFaceHash.binStart();
      vector<int> & binnedFaces = mFaceHash.binnedFaces();
      vector<face_ptr>   & faces = mMesh->get_faces();
      vector<vertex_ptr> & verts = mMesh->get_vertices();
      int i;
#ifdef _OPENMP
#pragma omp parallel for 
#endif 
      for(i = 0; i < verts.size(); i++){
	if(!verts[i]) continue;
	int b[3];
	vertex_ptr v = verts[i];
	coordinate_type p0 = v->coordinate();

	mFaceHash.nearestBin(p0, b);
	// int index = b[0] + b[1]*xRes + b[2]*xRes*yRes;	
	// int bStart = binStart[index];
	// int bEnd   = binStart[index + 1];
	// for(int i = bStart; i < bEnd; i++){
	for(int k = 0; k < sz; k++){
	  int xo = ip[3*k+0]; int yo = ip[3*k+1]; int zo = ip[3*k+2];
	  if(b[0] + xo < xRes && b[0] + xo > -1 &&
	     b[1] + yo < yRes && b[1] + yo > -1 &&
	     b[2] + zo < zRes && b[2] + zo > -1 ){
	    int offsetIndex = b[0]+xo + (b[1]+yo)*xRes + (b[2]+zo)*xRes*yRes;
	    int bStart0 = binStart[offsetIndex];
	    int bEnd0   = binStart[offsetIndex + 1];
	    blockCount += bStart0 - bEnd0;
	    for(int j = bStart0; j < bEnd0; j++){
	      int j0 = binnedFaces[j];
	      if(!mMesh->has_face(j0)) continue;
	      face_ptr f = faces[j0];
	      T        area = f->area();
	      int setPos = v->position_in_set();
	      coordinate_type p1   = f->center();
	      //std::cout << area << std::endl;
	      coordinate_type dp = (p0-p1);
	      T dist = norm(dp);
	      if(dist > 0.0){
		T dist3 = dist*dist*dist;

		T l3   = regLength*regLength*regLength;	     
		T kappa = (1.0-exp(-dist3/l3))/dist3;
		T l2 = regLength*regLength;
		T denom = powf((dist*dist + l2),1.5);
		//dp.normalize();
		coordinate_type fR = dp/denom;
		//coordinate_type fR = dp*kappa;
		coordinate_type ui = i4pi*area*cross(f->data,fR);
		u[setPos] += ui;
	      }
	    }
	  }
	}
	//	} //omp
      }
      return u;
    }
    void circulationToVorticity(face_ptr f){
      TIMER functionTimer(__FUNCTION__);
      T a = f->calc_area();
      if(a > 1e-16){
	face_vertex_ptr fv1 = f->fbegin();
	face_vertex_ptr fv2 = fv1->next();
	face_vertex_ptr fv3 = fv2->next();
	// coordinate_type e13 = fv1->vertex()->coordinate() - fv3->vertex()->coordinate();
	// coordinate_type e21 = fv2->vertex()->coordinate() - fv1->vertex()->coordinate();
	// coordinate_type e32 = fv3->vertex()->coordinate() - fv2->vertex()->coordinate();
	coordinate_type e13 = 
	  fv1->edge()->v1()->vertex()->coordinate() -
	  fv1->edge()->v2()->vertex()->coordinate();
	coordinate_type e21 = 
	  fv2->edge()->v1()->vertex()->coordinate() -
	  fv2->edge()->v2()->vertex()->coordinate();
	coordinate_type e32 = 
	  fv3->edge()->v1()->vertex()->coordinate() -
	  fv3->edge()->v2()->vertex()->coordinate();

	T e1 = fv1->data;
	T e2 = fv2->data;
	T e3 = fv3->data;
	f->data = 1.0/a*(e1*e13 + e2*e21 + e3*e32); 
      }
    }

    void vorticityToCirculation(face_ptr f){
      TIMER functionTimer(__FUNCTION__);
      T a = f->calc_area();
      if(a > 1e-16){
	face_vertex_ptr fv1 = f->fbegin();
	face_vertex_ptr fv2 = fv1->next();
	face_vertex_ptr fv3 = fv2->next();
	// coordinate_type e13 = fv1->vertex()->coordinate() - fv3->vertex()->coordinate();
	// coordinate_type e21 = fv2->vertex()->coordinate() - fv1->vertex()->coordinate();
	// coordinate_type e32 = fv3->vertex()->coordinate() - fv2->vertex()->coordinate();
	coordinate_type e13 = 
	  fv1->edge()->v1()->vertex()->coordinate() -
	  fv1->edge()->v2()->vertex()->coordinate();
	coordinate_type e21 = 
	  fv2->edge()->v1()->vertex()->coordinate() -
	  fv2->edge()->v2()->vertex()->coordinate();
	coordinate_type e32 = 
	  fv3->edge()->v1()->vertex()->coordinate() -
	  fv3->edge()->v2()->vertex()->coordinate();

	coordinate_type gamma = a*f->data2;
	if(norm(e13) < 1e-12 || norm(e21) < 1e-12 ||  norm(e32) < 1e-12){
	  fv1->data = 0.0;
	  fv2->data = 0.0;
	  fv3->data = 0.0;
	  return;
	}

	Eigen::MatrixXd E(4,3);  E = Eigen::MatrixXd::Zero(4,3);
	for(int i = 0; i < 3; i++){
	  E(i,0) = e13[i]; E(i,1) = e21[i]; E(i,2) = e32[i];
	}
	E(3,0) = 1;  E(3,1) = 1; E(3,2) = 1;
	Eigen::MatrixXd EtE(3,3);
	EtE = E.transpose()*E;
	Eigen::VectorXd b(4);
	b(0) = gamma[0];
	b(1) = gamma[1];
	b(2) = gamma[2];
	b(3) = 0;
	b = E.transpose()*b;
	Eigen::VectorXd c = EtE.ldlt().solve(b);
	//Eigen::VectorXd c = E.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
	m2::construct<SPACE> cons;
	if(c.squaredNorm() > 100.00){
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

    edge_ptr subdivideFace(face_vertex_ptr fv){
      T data3 = fv->face()->data3;
      TIMER functionTimer(__FUNCTION__);
      //assuming that the face vertex is the newly inserted one.
      face_vertex_ptr fv1 = fv; //this is the 
      face_vertex_ptr fv2 = fv->next();
      face_vertex_ptr fv3 = fv2->next();
      face_vertex_ptr fv4 = fv3->next();
#if 0
      coordinate_type e21 = fv2->vertex()->coordinate() - fv1->vertex()->coordinate();
      coordinate_type e32 = fv3->vertex()->coordinate() - fv2->vertex()->coordinate();
      coordinate_type e43 = fv4->vertex()->coordinate() - fv3->vertex()->coordinate();
      coordinate_type e14 = fv1->vertex()->coordinate() - fv4->vertex()->coordinate();
      coordinate_type e13 = fv1->vertex()->coordinate() - fv3->vertex()->coordinate();

      T e1 = fv1->data;      
      T e2 = fv2->data;
      T e3 = fv3->data;
      T e4 = fv4->data;
      coordinate_type gamma0 = e2*e21 + e3*e32 + e4*e43 + e1*e14; 
      Eigen::MatrixXd E(5,6);  E = Eigen::MatrixXd::Zero(5,6);
      for(int i = 0; i < 3; i++){
	E(i,0) = e21[i]; E(i,1) = e13[i]; E(i,2) = e32[i];
	E(i,3) = e14[i]; E(i,4) = e43[i]; E(i,5) = -e13[i];
      }
      E(3,0) = 1;  E(3,1) = 1; E(3,2) = 1;
      E(4,3) = 1;  E(4,4) = 1; E(4,5) = 1;
      Eigen::VectorXd b(5);
      b(0) = gamma0[0];
      b(1) = gamma0[1];
      b(2) = gamma0[2];
      b(3) = 0;
      b(4) = 0;
      b = E.transpose()*b;
      Eigen::MatrixXd EtE(3,3);
      EtE = E.transpose()*E;
      Eigen::VectorXd g = EtE.ldlt().solve(b);

      m2::construct<SPACE> cons;
      edge_ptr enew = cons.inser_edge(mMesh, fv1, fv3);
      enew->v1()->data = g(1);
      enew->v2()->data = g(5);
      fv3->data = g(0);
      fv2->data = g(2);
      fv1->data = g(3);
      fv4->data = g(4);
#endif
      m2::construct<SPACE> cons;
      edge_ptr enew = cons.insert_edge(mMesh, fv1, fv3);

      enew->v1()->face()->data3 = data3;
      enew->v2()->face()->data3 = data3;
      enew->v1()->data = 0.0;
      enew->v2()->data = 0.0;
      return enew;
    }

    void split_edge(edge_ptr e){
      TIMER functionTimer(__FUNCTION__);
      face_vertex_ptr fv1 = e->v1();
      face_vertex_ptr fv2 = e->v2();
      T circ1 = fv1->data;
      T circ2 = fv2->data;
      subdivide<SPACE> subd;
      //coordinate_type c = getButterflyWeight(e);
      vertex_ptr nv = subd.subdivide_edge(mMesh,e);
      //nv->coordinate() = c;
      fv1->data = circ2;
      fv2->data = circ1;
      fv1->next()->data = circ2;
      fv2->next()->data = circ1;

      nv->winding = 0.5*(e->v1()->vertex()->winding + 
			 e->v2()->vertex()->winding);
      subdivideFace(fv1->next());
      subdivideFace(fv2->next());
    }

    struct edge_sort {
      bool operator() (edge_ptr ei, edge_ptr ej) { return (ei->length() < ej->length());}
    } mEdgeSorter;

    bool insert_edges_curve(){
      TIMER functionTimer(__FUNCTION__);
      edge_array & edges = mMesh->get_edges();
      vector<vertex_ptr>& verts = mMesh->get_vertices();

      face_array& faces = mMesh->get_faces();
      int N = faces.size();

      bool topology_change = false;
      vector<T> edgeWeights;
      vector<T> vertexWeights;
      mMesh->update_all();
      mMesh->reset_flags();
      mMesh->pack();
      mesh_calculator<SPACE> calc;

      vector<edge_ptr> edgesToSplit;
#if 1
      for (int i = 0; i < edges.size(); i++){
      	if(!mMesh->has_edge(i)) continue;
	edge_ptr ei = edges[i];
	bool pinned =
	  ei->v1()->vertex()->pinned == true &&
	  ei->v2()->vertex()->pinned == true;
	if (pinned) 	  continue;
	T l = ei->length();
	if(l > 1.75*minLength){
	  edgesToSplit.push_back(ei);
	  continue;
	}
	// T l1n = ei->v1()->next()->edge()->length();
	// T l2n = ei->v2()->next()->edge()->length();
	// if(l/l1n > 3.0 || l/l2n > 3.0){
	//   if(l > 0.5*minLength)
	//     edgesToSplit.push_back(edges[i]);
	//   continue;
	// }

	// T l1p = ei->v1()->prev()->edge()->length();
	// T l2p = ei->v2()->prev()->edge()->length();
	// if(l/(l1n+l1p) > 0.85 || l/(l2n+l2p) > 0.85){
	//   if(l > 0.5*minLength)
	//     edgesToSplit.push_back(edges[i]);
	//   continue;
	// }
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

      std::cout << " - sorting "<< edgesToSplit.size() << " edges, ";
      std::sort (edgesToSplit.begin(), edgesToSplit.end(), mEdgeSorter);
      std::cout << "splitting "<< edgesToSplit.size() << " edges" << std::endl;
      for(int i = edgesToSplit.size(); i > 0; i--){
	this->split_edge(edgesToSplit[i-1]);
      }

      return topology_change;
    }

    bool collapse_edges(){
      TIMER functionTimer(__FUNCTION__);
      bool topology_change = false;
      edge_array  collectedEdges;
      vector<edge_ptr>& edges = mMesh->get_edges();
      mMesh->reset_flags();
      for (int i = 0; i < edges.size(); i++){

	if(!mMesh->has_edge(i)) continue;
	edge_ptr e = edges[i];
	if(e->flag == 1) continue;
	if(e->v1()->vertex()->flag == 1) continue;
	if(e->v2()->vertex()->flag == 1) continue;
	if(e->v2()->vertex()->pinned) continue;
	T dist = e->dist();
	if (dist < minCollapseLength){
	  e->v1()->vertex()->flag = 1;
	  e->v2()->vertex()->flag = 1;
	  collectedEdges.push_back(e);
	}
      }

      std::cout << " - deleting: " << collectedEdges.size() 
		<< " Tiny edges" << std::endl;

      for(int i = 0; i < collectedEdges.size(); i++){
	construct<SPACE> cons;
	if(!collectedEdges[i]) continue;
	edge_ptr e = collectedEdges[i];
	coordinate_type avg = 0.5*(e->v1()->vertex()->coordinate() +
				   e->v2()->vertex()->coordinate());
	e->v1()->vertex()->coordinate() = avg;
	//face_ptr nf = cons.delete_vertex_primitive(mMesh,e->v2()->vertex());
	//cons.collapse_edge_primitive(mMesh,e);	
	cons.collapse_edge(mMesh,e);
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

    /*
      bool delete_degenerate(){
      TIMER functionTimer(__FUNCTION__);
      bool topology_change = false;
      mMesh->reset_flags();
      vector<vertex_ptr>& verts = mMesh->get_vertices();
      vector<edge_ptr>& edges = mMesh->get_edges();
      int numEdges = 0;
      for (int i = 0; i < edges.size(); i++){
      if (!mMesh->has_edge(i)) continue;
      edge_ptr e = edges[i];
      int s1 = e->v1()->face()->size();
      int s2 = e->v2()->face()->size();

      if(s1 < 3    || s2 < 3     || 
      e->v1()->vertex() == e->v2()->vertex()){
      construct<SPACE> cons;
      edge_ptr e = edges[i];
      face_ptr nf = cons.delete_edge(mMesh,e);

      numEdges++;
      topology_change = true;
      };
      }

      //if(verts.size() > 0) topology_change = true;
      int numVerts = 0;
      for (int i = 0; i < verts.size(); i++){
      if(!mMesh->has_vertex(i)) continue;
      if(verts[i]->pinned) continue;
      vertex_ptr v = verts[i];
      vertex_ptr v1 = v->fbegin()->next()->vertex();
      vertex_ptr v2 = v->fbegin()->prev()->vertex();

      if(verts[i]->size() == 2){
      //std::cout << "degenerate vert: " << verts[i]->size() << std::endl;
      int sz = verts[i]->size();
      construct<SPACE> cons;

      edge_ptr e = v->fbegin()->next()->edge();
      cons.delete_edge(mMesh,e);	  

      face_ptr nf = cons.delete_vertex_primitive(mMesh,v);
      numVerts++;
      topology_change = true;
      continue;
      } 

      if(verts[i]->size() == 1){
      //std::cout << "degenerate vert: " << verts[i]->size() << std::endl;
      int sz = verts[i]->size();
      construct<SPACE> cons;
      face_ptr nf = cons.delete_vertex_primitive(mMesh,v);
      // v1->print();
      // v2->print();
      numVerts++;
      topology_change = true;
      continue;
      } 
      }
      std::cout << "deleted " << numEdges << " degenerate edges and " << numVerts << " 2-verts"<< std::endl;
      if (topology_change) {
      mMesh->pack();
      mMesh->update_all();
      mMesh->reset_flags();
      }      
      return topology_change;
      }
    */

    bool delete_degenerate(){
      TIMER functionTimer(__FUNCTION__);
      bool topology_change = false;
      vector<vertex_ptr>& verts = mMesh->get_vertices();
      vector<edge_ptr>& edges = mMesh->get_edges();
      int zeroEdges = 0;
      int twoFaces = 0;
      int flatVol = 0;
      vector<edge_ptr> edgesToDelete;
      for (int i = 0; i < edges.size(); i++){
	if (mMesh->has_edge(i)) {
	  edge_ptr e = edges[i];
	  if(e->flag == 1) continue;
	  int s1 = e->v1()->face()->size();
	  int s2 = e->v2()->face()->size();
	  T a1 = e->v1()->face()->calc_area();
	  T a2 = e->v2()->face()->calc_area();

	  if(a1 < 1e-12 || a2 < 1e-12) { 
	    e->flag = 1;
	    edgesToDelete.push_back(e);
	    continue;
	  }

	  if(s1 < 3 || s2 < 3) { 
	    twoFaces++;
	    e->flag = 1;
	    edgesToDelete.push_back(e);
	    continue;
	  }

	  if(e->v1()->face() == e->v2()->face()) { 
	    edgesToDelete.push_back(e);
	    e->flag = 1;
	    continue;
	  }

	  if(e->v1()->vertex() == e->v2()->vertex()){
	    zeroEdges++;
	    e->flag = 1;
	    edgesToDelete.push_back(e);
	    continue;
	  };

	  if(e->v1()->prev()->vertex() == e->v2()->prev()->vertex()){
	    flatVol++;
	    edgesToDelete.push_back(e);
	    e->flag = 1;
	    continue;
	  };
	}
      }

      topology_change = edgesToDelete.size() > 0;
      construct<SPACE> cons;
      for(int i = 0; i < edgesToDelete.size(); i++){
	edge_ptr e = edgesToDelete[i];
	face_ptr nf = cons.delete_edge(mMesh,e);
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
      for (int i = 0; i < verts.size(); i++){
	if(!mMesh->has_vertex(i)) continue;
	if(verts[i]->pinned) continue;
	vertex_ptr v = verts[i];

	if(verts[i]->size() < 3){
	  construct<SPACE> cons;
	  face_ptr nf = cons.delete_vertex_primitive(mMesh,v);
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
#endif

      if (topology_change) {
	m2::remesh<space3> rem;
	//rem.triangulate(mMesh);
	mMesh->pack();
	mMesh->update_all();
	mMesh->reset_flags();
      }      
      return topology_change;
    }

    bool checkCompatibility(edge_ptr ei, edge_ptr ej){
      int vertIds[4];
      vertIds[0] = ei->v1()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[2] = ei->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

      vertIds[1] = ei->v2()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[5] = ei->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

      vertIds[2] = ej->v1()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[8] = ej->v1()->face()->fbegin()->prev()->vertex()->position_in_set();

      vertIds[3] = ej->v2()->face()->fbegin()->vertex()->position_in_set();
      //vertIds[11] = ej->v2()->face()->fbegin()->prev()->vertex()->position_in_set();

      for(int i = 0; i < 4; i++)
	for(int j = 0; j < 4; j++){
	  if(i != j)
	    if (vertIds[i] == vertIds[j]) return false;
	}

      return true;
    };
    
    vector<edge_ptr> getPairList(){
      TIMER functionTimer(__FUNCTION__);

      mMesh->reset_flags();
      this->hashEdges();
      T dx = mEdgeHash.dx();
      int xRes = mEdgeHash.xRes();
      int yRes = mEdgeHash.yRes();
      int zRes = mEdgeHash.zRes();
      int sz = 8;
      std::cout << " edge pair query sphere size: "<< sz << std::endl;
      int * ip = new int[3*sz];
      int ind = 0;
      for(int k = 0; k < 2; k++)
	for(int j = 0; j < 2; j++)
	  for(int i = 0; i < 2; i++){
	    int index = i + j*2 + k*4;
	    ip[3*index + 0] = i;
	    ip[3*index + 1] = j;
	    ip[3*index + 2] = k;
	    ind++;
	  }

      T imax = dx/xRes;
      vector<int> & binStart    = mEdgeHash.binStart();
      vector<int> & binnedVerts = mEdgeHash.binnedEdges();
      vector<edge_ptr>   & edges = mMesh->get_edges();
      vector<edge_ptr> edgePairs;

      std::cout << " getting edge pairs: " << edges.size() << std::endl;

      for(int z = 0; z < zRes; z++)
	for(int y = 0; y < yRes; y++)
	  for(int x = 0; x < xRes; x++){
	    int index  = x + y*xRes + z*xRes*yRes;
	    int bStart0 = binStart[index];
	    int bEnd0   = binStart[index + 1];

	    for(int i = bStart0; i < bEnd0; i++){
	      int i0 = binnedVerts[i];
	      int jm = -1;
	      
	      if (!mMesh->has_edge(i0))   continue;
	      if (edges[i0]-> flag == 1) continue;
	      bool breakLoop = false;
	      edge_ptr ei = edges[i0];
	      //if((ei->v1()->vertex()->winding + ei->v1()->vertex()->winding) < 2.0*M_PI) continue;
	      for(int k = 0; k < sz; k++){
		//for(int k = 0; k < 0; k++){
		int xo = ip[3*k+0]; int yo = ip[3*k+1]; int zo = ip[3*k+2];

		if(x + xo < xRes && x + xo > -1 &&
		   y + yo < yRes && y + yo > -1 &&
		   z + zo < zRes && z + zo > -1 ){
		  int offsetIndex = x+xo + (y+yo)*xRes + (z+zo)*xRes*yRes;
		  int bStart1 = binStart[offsetIndex];
		  int bEnd1   = binStart[offsetIndex + 1];
		  for(int j = bStart1; j < bEnd1; j++){
		    int j0 = binnedVerts[j];
		    if(j0 == i0)             continue;
		    if(!mMesh->has_edge(j0)) continue;
		    edge_ptr ej = edges[j0];

		    int Numi1 = ei->v1()->face()->size();
		    int Numi2 = ei->v2()->face()->size();
		    int Numj1 = ei->v1()->face()->size();
		    int Numj2 = ei->v2()->face()->size();

		    if(Numi1 != 3) continue;
		    if(Numi2 != 3) continue;
		    if(Numj1 != 3) continue;
		    if(Numj2 != 3) continue;

		    int setPos = ej->position_in_set();
		    //battery of degeneracy tests:
		    
		    if(ej->flag == 1 || ei->flag == 1) continue;
		    bool compatible = checkCompatibility(ei,ej);
		    if(!compatible) continue;
		    if(ei->v1()->face()->calc_area() < 0.1*minLength*minLength) continue;
		    if(ei->v2()->face()->calc_area() < 0.1*minLength*minLength) continue;
		    if(ej->v1()->face()->calc_area() < 0.1*minLength*minLength) continue;
		    if(ej->v2()->face()->calc_area() < 0.1*minLength*minLength) continue;

		    distance_calculator<SPACE> dcalc;
		    T d = dcalc.calcEdgeEdgeDistance(ei,ej);
		    coordinate_type Ni = 
		      0.5*ei->v1()->face()->normal() +
		      0.5*ei->v2()->face()->normal();
		    coordinate_type Nj = 
		      0.5*ej->v1()->face()->normal() + 
		      0.5*ej->v2()->face()->normal();
		    Ni.normalize();
		    Nj.normalize();
		    T dNij = dot(Ni,Nj);
		    coordinate_type Vi = 0.5*ei->v1()->face()->data + 
		      0.5*ei->v2()->face()->data;
		    coordinate_type Vj = 0.5*ej->v1()->face()->data + 
		      0.5*ej->v2()->face()->data;
		    T dVij = norm(Vi+Vj)/norm(Vi);
		    T avgWinding = 0.25*(ei->v1()->vertex()->winding +
					 ei->v1()->vertex()->winding + 
					 ej->v1()->vertex()->winding +
					 ej->v1()->vertex()->winding); 
		    //   avgWinding > 2.0*M_PI
		    if(d < edgeJoinThresh && dNij < -0.98){
		      jm = j0;
		    }
		  }
		}
	      }
	      if(jm > 0){
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
    
    void pipeFace(face_ptr f0, 
		  face_ptr f1,
		  T d00, T d01,
		  T d10, T d11){
      
      face_vertex_ptr v0b = f0->fbegin();
      face_vertex_ptr v0e = f0->fend();
      face_vertex_ptr v1b = f1->fbegin();
      face_vertex_ptr v1e = f1->fend();
      T mine = 9999;
      face_vertex_ptr v0s = v0b;
      face_vertex_ptr v1s = v1b;
      bool iterating0 = true;
      while(iterating0){
	iterating0 = v0b != v0e;
	bool iterating1 = true;
	while(iterating1){
	  iterating1 = v1b != v1e;
	  T d0 = norm(v0b->coordinate() - v1b->coordinate());
	  T d1 = norm(v0b->next()->next()->coordinate() - 
			 v1b->next()->next()->coordinate());
	  T  e = d0*d0 + d1*d1;
	  if(e < mine){
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
      vector<vertex_ptr>      verts;
      while(iterating0){
	iterating0 = v0b != v0e;
	if(v0s->vertex() != v1s->vertex()){
	  pairs.push_back(v0s);
	  pairs.push_back(v1s);
	  verts.push_back(v0s->vertex());
	  verts.push_back(v1s->vertex());
	}
	else verts.push_back(v0s->vertex());
	v0s = v0s->next();
	v1s = v1s->prev();
	v0b = v0b->next();
      }

      for(int i = 0; i < pairs.size(); i+=2){
	construct<SPACE> cons;
	edge_ptr ei = cons.insert_edge(mMesh, pairs[i], pairs[i+1]);
	ei->v1()->face()->color.g = 1.0;
	ei->v2()->face()->color.g = 1.0;
	if(i==0) ei->v1()->data = d00; ei->v2()->data = d01;
	if(i==2) ei->v1()->data = d10; ei->v2()->data = d11;
      }

      // for(int k = 0; k < 10; k++){
      // 	mesh_filter<SPACE> filter;
      // 	for(int i = 0; i < verts.size(); i++){
      // 	  verts[i]->coordinate() +=
      // 	    0.1*filter.laplacianFilterVertex(verts[i]);
      // 	}
      // }
    }

    bool joinEdges(){
      bool topology_change = false;
      vector<edge_ptr> edgePairs = this->getPairList();
      if(edgePairs.size() > 0) topology_change = true;
      
      std::cout << " - joining: " << edgePairs.size() << " pairs" << std::endl;
      for(int i = 0; i < edgePairs.size(); i+=2){
	edge_ptr e0 = edgePairs[i+0];
	edge_ptr e1 = edgePairs[i+1];

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
	if(v00->vertex()->pinned) continue;
	if(v01->vertex()->pinned) continue;
	if(v10->vertex()->pinned) continue;
	if(v11->vertex()->pinned) continue;
	face_ptr f0 = cons.delete_edge(mMesh,e0);
	face_ptr f1 = cons.delete_edge(mMesh,e1);
	//std::cout << f0->size() << " " << f1->size() << std::endl;
	pipeFace(f0,f1,dt00,dt01,dt10,dt11);

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
    
    

    void drawVorticity(){
      face_array faces = mMesh->get_faces(); 
      for(int i = 0; i < faces.size(); i++){
	face_ptr f = faces[i];
	if(f){
	  coordinate_type w = f->data;
	  coordinate_type cen = f->center();
	  glBegin(GL_LINES);
	  T c = 0.005;
	  glColor3f(0.0,0.0,1.0);
	  glVertex3d(cen[0],cen[1],cen[2]);
	  glVertex3d(cen[0]+c*w[0],cen[1]+c*w[1],cen[2]+c*w[2]);
	  glEnd();
	}
      }
    }

    void drawVelocity(){
      vertex_array verts = mMesh->get_vertices(); 
      for(int i = 0; i < verts.size(); i++){
	vertex_ptr v = verts[i];
	if(v){
	  coordinate_type vel = v->data;
	  coordinate_type cen = v->coordinate();
	  glBegin(GL_LINES);
	  T c = 0.01;
	  glColor3f(0.0,1.0,0.0);
	  glVertex3d(cen[0],cen[1],cen[2]);
	  glVertex3d(cen[0]+c*vel[0],cen[1]+c*vel[1],cen[2]+c*vel[2]);
	  glEnd();
	}
      }
    }

    edge_bin<SPACE>   mEdgeHash;
    pole_tree<SPACE>  mOctree;
    //aabb_tree<SPACE,3>  mAaBb;
    face_bin<SPACE>   mFaceHash;
    vertex_bin<SPACE> mVertexHash;
    T maxCurvature, minCurvature, minLength, minCollapseLength, maxLength, regLength, edgeJoinThresh;
    vector<coordinate_type> mEdgeCirculation;
    vector<coordinate_type> mFaceVorticity;
    control_ptr mMesh;
  };
}//end namespace
#endif

