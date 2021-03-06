/*
 Solving the Nearest Point-on-Curve Problem 
 and
 A Bezier Curve-Based Root-Finder
 by Philip J. Schneider
 from "Graphics Gems", Academic Press, 1990
 */

/*	point_on_curve.c	*/		

#ifndef __BEZIER_NEAREST_POINT__
#define __BEZIER_NEAREST_POINT__

#include <iostream>


/* find minimum of a and b */
#define MIN(a,b)	(((a)<(b))?(a):(b))	

/* find maximum of a and b */
#define MAX(a,b)	(((a)>(b))?(a):(b))	

/* swap a and b (see Gem by Wyvill) */
#define SWAP(a,b)	{ a^=b; b^=a; a^=b; }

/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) */
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

/* clamp the input to the specified range */
#define CLAMP(v,l,h)	((v)<(l) ? (l) : (v) > (h) ? (h) : v)

#define	EPSILON	(ldexp(1.0,-MAXDEPTH-1)) /*Flatness control value */
#define	DEGREE	3			/*  Cubic Bezier curve		*/
#define	W_DEGREE 5			/*  Degree of eqn to find roots of */

template <typename T>
inline T distSqr(Eigen::Matrix<T,2,1> p0,
		 Eigen::Matrix<T,2,1> p1){
	return (p1 -p0).magSqr();
};

template <typename T>
class bezier_nearest_point{
	
	typedef Eigen::Matrix<T,2,1> ctype;
public:
	bezier_nearest_point(){
		MAXDEPTH = 64;
	}
	/*
	 *  NearestPointOnCurve :
	 *  	Compute the parameter value of the point on a Bezier
	 *		curve segment closest to some arbtitrary, user-input point.
	 *		Return the point on the curve at that parameter value.
	 *
	 */
	
	ctype NearestPointOnCurve(ctype P, ctype * V, T& dist, T& t)
	// 	P;			/* The user-supplied point	  */
	//  V;			/* Control points of cubic Bezier */
	{
		ctype	*w;			/* Ctl pts for 5th-degree eqn	*/
		T		t_candidate[W_DEGREE];	/* Possible roots		*/     
		int 	n_solutions;		/* Number of roots found	*/	
		
		/*  Convert problem to 5th-degree Bezier form	*/
		w = ConvertToBezierForm(P, V);
		
		/* Find all possible roots of 5th-degree equation */
		n_solutions = FindRoots(w, W_DEGREE, t_candidate, 0);
		free((char *)w);
		
		/* Compare distances of P to all candidates, and to t=0, and t=1 */
		{
			T new_dist;
			ctype 	p;
			ctype  v;
			int		i;
			
			
			/* Check distance to beginning of curve, where t = 0	*/
			dist = distSqr(P,V[0]);
			t = 0.0;
			
			/* Find distances for candidate points	*/
			for (i = 0; i < n_solutions; i++) {
				p = Bezier(V, DEGREE, t_candidate[i],
						   (ctype *)NULL, (ctype *)NULL);
				new_dist = distSqr(P,p);
				if (new_dist < dist) {
					dist = new_dist;
					t = t_candidate[i];
				}
			}
			
			/* Finally, look at distance to end point, where t = 1.0 */
			new_dist = distSqr(P, V[DEGREE]);
			if (new_dist < dist) {
				dist = new_dist;
				t = 1.0;
			}
		}
		
		/*  Return the point on the curve at parameter value t */
		//printf("t : %4.12f\n", t);
		return (Bezier(V, DEGREE, t, (ctype *)NULL, (ctype *)NULL));
	}
	
	
	/*
	 *  ConvertToBezierForm :
	 *		Given a point and a Bezier curve, generate a 5th-degree
	 *		Bezier-format equation whose solution finds the point on the
	 *      curve nearest the user-defined point.
	 */

	ctype *ConvertToBezierForm(ctype P, ctype  *V)
	// 	P;			/* The point to find t for	*/
	//	V;			/* The control points		*/
	{
		int 	i, j, k, m, n, ub, lb;	
		int 	row, column;		/* Table indices		*/
		ctype 	c[DEGREE+1];		/* V(i)'s - P			*/
		ctype 	d[DEGREE];		/* V(i+1) - V(i)		*/
		ctype 	*w;			/* Ctl pts of 5th-degree curve  */
		T 	cdTable[3][4];		/* Dot product of c, d		*/
		//TODO: make this a member
		T z[3][4] = {	/* Precomputed "z" for cubics	*/
			{1.0, 0.6, 0.3, 0.1},
			{0.4, 0.6, 0.6, 0.4},
			{0.1, 0.3, 0.6, 1.0},
		};
		
		
		/*Determine the c's -- these are vectors created by subtracting*/
		/* point P from each of the control points				*/
		for (i = 0; i <= DEGREE; i++) {
			c[i] = V[i] - P;

		}
		/* Determine the d's -- these are vectors created by subtracting*/
		/* each control point from the next					*/
		for (i = 0; i <= DEGREE - 1; i++) { 
			d[i] = V[i+1] - V[i];
			d[i] = V2ScaleII(d[i], 3.0);
		}
		
		/* Create the c,d table -- this is a table of dot products of the */
		/* c's and d's							*/
		for (row = 0; row <= DEGREE - 1; row++) {
			for (column = 0; column <= DEGREE; column++) {
				cdTable[row][column] = d[row][0]*c[column][0] + d[row][1]*c[column][1];
			}
		}
		
		/* Now, apply the z's to the dot products, on the skew diagonal*/
		/* Also, set up the x-values, making these "points"		*/
		w = (ctype *)malloc((unsigned)(W_DEGREE+1) * sizeof(ctype));
		for (i = 0; i <= W_DEGREE; i++) {
			w[i][1] = 0.0;
			w[i][0] = (T)(i) / W_DEGREE;
		}
		
		n = DEGREE;
		m = DEGREE-1;
		for (k = 0; k <= n + m; k++) {
			lb = MAX(0, k - m);
			ub = MIN(k, n);
			for (i = lb; i <= ub; i++) {
				j = k - i;
				w[i+j][1] += cdTable[j][i] * z[j][i];
			}
		}
		
		return (w);
	}
	
	
	/*
	 *  FindRoots :
	 *	Given a 5th-degree equation in Bernstein-Bezier form, find
	 *	all of the roots in the interval [0, 1].  Return the number
	 *	of roots found.
	 */

	int FindRoots(ctype  *w,
						 int			degree, 
						 T				*t,
						 int			depth)
//	w;			/* The control points		*/
//	degree;		/* The degree of the polynomial	*/
//	t;			/* RETURN candidate t-values	*/
//	depth;		/* The depth of the recursion	*/
	{  
		int 	i;
		ctype 	Left[W_DEGREE+1],	/* New left and right 		*/
				Right[W_DEGREE+1];	/* control polygons		*/
		int 	left_count,		/* Solution count from		*/
				right_count;		/* children			*/
		T 	left_t[W_DEGREE+1],	/* Solutions from kids		*/
			right_t[W_DEGREE+1];
		
		switch (CrossingCount(w, degree)) {
			case 0 : {	/* No solutions here	*/
				return 0;	
			}
			case 1 : {	/* Unique solution	*/
				/* Stop recursion when the tree is deep enough	*/
				/* if deep enough, return 1 solution at midpoint 	*/
				if (depth >= MAXDEPTH) {
					t[0] = (w[0][0] + w[W_DEGREE][0]) / 2.0;
					return 1;
				}
				if (ControlPolygonFlatEnough(w, degree)) {
					t[0] = ComputeXIntercept(w, degree);
					return 1;
				}
				break;
			}
		}
		
		/* Otherwise, solve recursively after	*/
		/* subdividing control polygon		*/
		Bezier(w, degree, 0.5, Left, Right);
		left_count  = FindRoots(Left,  degree, left_t, depth+1);
		right_count = FindRoots(Right, degree, right_t, depth+1);
		
		
		/* Gather solutions together	*/
		for (i = 0; i < left_count; i++) {
			t[i] = left_t[i];
		}
		for (i = 0; i < right_count; i++) {
			t[i+left_count] = right_t[i];
		}
		
		/* Send back total number of solutions	*/
		return (left_count+right_count);
	}
	
	
	/*
	 * CrossingCount :
	 *	Count the number of times a Bezier control polygon 
	 *	crosses the 0-axis. This number is >= the number of roots.
	 *
	 */

	int CrossingCount(ctype	*V,int degree)
	//V;			/*  Control pts of Bezier curve	*/
	//degree;		/*  Degreee of Bezier curve 	*/
	{
		int 	i;	
		int 	n_crossings = 0;	/*  Number of zero-crossings	*/
		int		sign, old_sign;		/*  Sign of coefficients	*/
		old_sign = V[0][1] == 0 ? 0:  V[0][1] > 0 ? 1:-1;
		for (i = 1; i <= degree; i++) {
			sign = V[i][1] == 0 ? 0:  V[i][1] > 0 ? 1:-1;
			if (sign != old_sign) n_crossings++;
			old_sign = sign;
		}
		return n_crossings;
	}
	
	
	
	/*
	 *  ControlPolygonFlatEnough :
	 *	Check if the control polygon of a Bezier curve is flat enough
	 *	for recursive subdivision to bottom out.
	 *
	 *  Corrections by James Walker, jw@jwwalker.com, as follows:
	 
	 There seem to be errors in the ControlPolygonFlatEnough function in the
	 Graphics Gems book and the repository (NearestPoint.c). This function
	 is briefly described on p. 413 of the text, and appears on pages 793-794.
	 I see two main problems with it.
	 
	 The idea is to find an upper bound for the error of approximating the x
	 intercept of the Bezier curve by the x intercept of the line through the
	 first and last control points. It is claimed on p. 413 that this error is
	 bounded by half of the difference between the intercepts of the bounding
	 box. I don't see why that should be true. The line joining the first and
	 last control points can be on one side of the bounding box, and the actual
	 curve can be near the opposite side, so the bound should be the difference
	 of the bounding box intercepts, not half of it.
	 
	 Second, we come to the implementation. The values distance[i] computed in
	 the first loop are not actual distances, but squares of distances. I
	 realize that minimizing or maximizing the squares is equivalent to
	 minimizing or maximizing the distances.  But when the code claims that
	 one of the sides of the bounding box has equation
	 a * x + b * y + c + max_distance_above, where max_distance_above is one of
	 those squared distances, that makes no sense to me.
	 
	 I have appended my version of the function. If you apply my code to the
	 cubic Bezier curve used to test NearestPoint.c,
	 
	 static ctype bezCurve[4] = {    /  A cubic Bezier curve    /
	 { 0.0, 0.0 },
	 { 1.0, 2.0 },
	 { 3.0, 3.0 },
	 { 4.0, 2.0 },
	 };
	 
	 my code computes left_intercept = -3.0 and right_intercept = 0.0, which you
	 can verify by sketching a graph. The original code computes
	 left_intercept = 0.0 and right_intercept = 0.9.
	 
	 */
	
	/* static int ControlPolygonFlatEnough( const ctype* V, int degree ) */

	int ControlPolygonFlatEnough(ctype	*V, int degree)
//	V;		/* Control points	*/
//	degree;		/* Degree of polynomial	*/
	{
		int     i;        /* Index variable        */
		T  value;
		T  max_distance_above;
		T  max_distance_below;
		T  error;        /* Precision of root        */
		T  intercept_1,
		intercept_2,
		left_intercept,
		right_intercept;
		T  a, b, c;    /* Coefficients of implicit    */
		/* eqn for line from V[0]-V[deg]*/
		T  det, dInv;
		T  a1, b1, c1, a2, b2, c2;
		
		/* Derive the implicit equation for line connecting first */
		 /*  and last control points */
		a = V[0][1] - V[degree][1];
		b = V[degree][0] - V[0][0];
		c = V[0][0] * V[degree][1] - V[degree][0] * V[0][1];
		
		max_distance_above = max_distance_below = 0.0;
		
		for (i = 1; i < degree; i++)
		{
			value = a * V[i][0] + b * V[i][1] + c;
			
			if (value > max_distance_above)
			{
				max_distance_above = value;
			}
			else if (value < max_distance_below)
			{
				max_distance_below = value;
			}
		}
		
		/*  Implicit equation for zero line */
		a1 = 0.0;
		b1 = 1.0;
		c1 = 0.0;
		
		/*  Implicit equation for "above" line */
		a2 = a;
		b2 = b;
		c2 = c - max_distance_above;
		
		det = a1 * b2 - a2 * b1;
		dInv = 1.0/det;
		
		intercept_1 = (b1 * c2 - b2 * c1) * dInv;
		
		/*  Implicit equation for "below" line */
		a2 = a;
		b2 = b;
		c2 = c - max_distance_below;
		
		det = a1 * b2 - a2 * b1;
		dInv = 1.0/det;
		
		intercept_2 = (b1 * c2 - b2 * c1) * dInv;
		
		/* Compute intercepts of bounding box    */
		left_intercept = MIN(intercept_1, intercept_2);
		right_intercept = MAX(intercept_1, intercept_2);
		
		error = right_intercept - left_intercept;
		
		return (error < EPSILON)? 1 : 0;
	}
	
	
	/*
	 *  ComputeXIntercept :
	 *	Compute intersection of chord from first control point to last
	 *  	with 0-axis.
	 * 
	 */
	/* NOTE: "T" and "Y" do not have to be computed, and there are many useless
	 * operations in the following (e.g. "0.0 - 0.0").
	 */
	T ComputeXIntercept(ctype 	*V, int		degree)
	//V;			/*  Control points	*/
	//degree; 		/*  Degree of curve	*/
	{
		T	XLK, YLK, XNM, YNM, XMK, YMK;
		T	det, detInv;
		T	S;
		T	_X, _Y;
		
		XLK = 1.0 - 0.0;
		YLK = 0.0 - 0.0;
		XNM = V[degree][0] - V[0][0];
		YNM = V[degree][1] - V[0][1];
		XMK = V[0][0] - 0.0;
		YMK = V[0][1] - 0.0;
		
		det = XNM*YLK - YNM*XLK;
		detInv = 1.0/det;
		
		S = (XNM*YMK - YNM*XMK) * detInv;
		/*  T = (XLK*YMK - YLK*XMK) * detInv; */
		
		_X = 0.0 + XLK * S;
		/*  Y = 0.0 + YLK * S; */
		
		return _X;
	}
	
	
	/*
	 *  Bezier : 
	 *	Evaluate a Bezier curve at a particular parameter value
	 *      Fill in control points for resulting sub-curves if "Left" and
	 *	"Right" are non-null.
	 * 
	 */
	static ctype Bezier(ctype 	*V, 
							   int				degree,
							   T				t, 
							   ctype 	*Left, 
							   ctype 	*Right)
	//degree;		/* Degree of bezier curve	*/
	//V;			/* Control pts			*/
	//t;			/* Parameter value		*/
	//Left;			/* RETURN left half ctl pts	*/
	//Right;		/* RETURN right half ctl pts	*/
	{
		int 	i, j;		/* Index variables	*/
		ctype 	Vtemp[W_DEGREE+1][W_DEGREE+1];
		
		
		/* Copy control points	*/
		for (j =0; j <= degree; j++) {
			Vtemp[0][j] = V[j];
		}
		
		/* Triangle computation	*/
		for (i = 1; i <= degree; i++) {	
			for (j =0 ; j <= degree - i; j++) {
				Vtemp[i][j][0] =
				(1.0 - t) * Vtemp[i-1][j][0] + t * Vtemp[i-1][j+1][0];
				Vtemp[i][j][1] =
				(1.0 - t) * Vtemp[i-1][j][1] + t * Vtemp[i-1][j+1][1];
			}
		}
		
		if (Left != NULL) {
			for (j = 0; j <= degree; j++) {
				Left[j]  = Vtemp[j][0];
			}
		}
		if (Right != NULL) {
			for (j = 0; j <= degree; j++) {
				Right[j] = Vtemp[degree-j][j];
			}
		}
		
		return (Vtemp[degree][0]);
	}
	
	ctype V2ScaleII(ctype	v, T	s)
	//v;
	//s;
	{
		ctype result;
		result[0] = v[0] * s;
		result[1] = v[1] * s;
		return (result);
	}
	
	int		MAXDEPTH;	/*  Maximum depth for recursion */
};
#endif
