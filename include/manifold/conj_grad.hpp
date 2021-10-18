/*
 *  conj_grad.hpp
 *  Manifold
 *
 *  Created by John Delaney on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __CONJ_GRAD__
#define __CONJ_GRAD__

#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>

template <typename MAT, typename T>
void	cgsolve(MAT& A, vector<T>& x, vector<T>& b){
	//conjugate gradient solver... not preconditioned
	//depending on the data type, reqires copy into vector x
	
	int nthreads, tid, i, chunk;
	chunk = 10;
	T tol = 1; T alpha,beta,r0,r1,eps;
	eps = 0.000001;
	long N = x.size();
	vector<T> R;  R.resize(N);
	vector<T> p;  p.resize(N);
	vector<T> Ap; Ap.resize(N);
	
	T RR = 0;
	
	//we'll have to let the user parallelize the multiplication... though could do it via a kernel, for now this works :(
	A.mult(x,Ap);
//#pragma omp parallel private(i,tid)	
	//#pragma omp parallel shared(b,R,xt,x,p,nthreads) private(i,tid)
	{
	  //	tid = omp_get_thread_num();
	  //	if (tid == 0)
	  //	{
	  //		nthreads = omp_get_num_threads();
	  //	}
//#pragma omp for schedule(dynamic,chunk)
		for (i = 0; i < N; i++){
			R[i] = b[i] - Ap[i];
			p[i] = R[i];
			RR  += R[i]*R[i];
		}		
	}
	
	r0 = RR;
	r1 = r0;
	/*
	 ----------------------------------------------------------------------------------------
	 begin iterating	
	 ----------------------------------------------------------------------------------------*/	
	int counter = 0;

	while (counter < N){		
		T pAp = 0;			
		RR    = 0;
		alpha = 0;
		beta  = 0;
		
		A.mult(p,Ap);
#ifdef _OPENMP		
//#pragma omp parallel private(i,tid)			
#endif
		{	
#ifdef _OPENMP
			tid = omp_get_thread_num();
			if (tid == 0) nthreads = omp_get_num_threads();	

//#pragma omp for schedule(dynamic,chunk)
#endif		
			for (i = 0; i < N; i++) {			
				pAp += p[i]*Ap[i];
			}
		};
		
		alpha = r1/(pAp);

//#pragma omp parallel private(i,tid)			

		{
#ifdef _OPENMP
			tid = omp_get_thread_num();
			if (tid == 0) nthreads = omp_get_num_threads();	
#endif				
//#pragma omp for schedule(dynamic,chunk)
			for (i = 0; i < N; i++) {				
				x[i] +=  p[i]*alpha;			
				R[i] -=  Ap[i]*alpha;
				RR   +=  R[i]*R[i];				
			}
		}
		r1 = RR;
		beta = r1/r0;
		double rnorm = r1.norm2();
		//if (r1 <= r0*0.1) counter = N;
		if (rnorm == 0) counter = N;		
		if (rnorm < 0.00001) counter = N;				
		//if (r1 > r0) counter = N;
		
//		cout << "r:       " << counter << ": " << rnorm << endl;
		//cout << "beta:    " << counter << ": " << beta[0] << endl;		

//		#pragma omp parallel private(i,tid)	
		////#pragma omp parallel shared(b,R,xt,x,p,nthreads) private(i,tid)
		{
#ifdef _OPENMP
			tid = omp_get_thread_num();
			if (tid == 0) nthreads = omp_get_num_threads();			
//			#pragma omp for schedule(dynamic,chunk)
#endif
			for (int i = 0; i < N; i++){
				p[i] = R[i] + p[i]*beta;
				//cout << "p:    " << counter << ": " << p[i][0] << " " << p[i][1] << " " << p[i][2] << endl;
			}
		}
		counter++;
		r0 = r1;
	}	
}

#endif

