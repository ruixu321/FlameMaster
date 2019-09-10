#ifndef __Spline__
#define __Spline__

#if defined (applec) || defined (powerc)
#pragma once
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif


/* typedefs */
typedef struct Spline {
	int		n;			/* number of points									*/
	int		storage;	/* signals external storage							*/
	int		xStorage;	/* signals external storage for x only				*/
	Double	*a;			/* coefficient a of ( s = a + b x + c x^2 + d x^3 )	*/
	Double	*b;			/* coefficient b of ( s = a + b x + c x^2 + d x^3 )	*/
	Double	*c;			/* coefficient c of ( s = a + b x + c x^2 + d x^3 )	*/
	Double	*d;			/* coefficient d of ( s = a + b x + c x^2 + d x^3 )	*/
	Double	*x;			/* points											*/
} Spline, *SplinePtr;


typedef struct PSpline {
	int		n;			/* number of points									*/
	int		storage;	/* signals external storage							*/
	Spline	*xs;		/* x-direction spline								*/
	Spline	*ys;		/* y-direction spline								*/
	Double	*u;			/* parameter of both splines						*/
	Double	arcLength;	/* arc length of the parametric spline				*/
} PSpline, *PSplinePtr;


/* prototypes */
void FreeSpline( SplinePtr sp );
SplinePtr ComputeSimpleSpline( Double *x, Double *f, int nPoints, 
	int freeLeft, Double sL, int freeRight, Double sR, void *storage, int bSpline );
void SplineInterpolate( SplinePtr s, Double *x, Double *f, int points );
int Locate( const Double xx[], int n, Double x );
int SplineRandomAccess( SplinePtr s, Double x, Double *value );
int SplineSequentialAccess( SplinePtr s, Double x, Double *value );
int SplineDerivative( SplinePtr s, Double x, Double *deriv );
void PrintSpline( SplinePtr s, char *format, FILE *fp );
void WriteSimpleSpline( SplinePtr sp, FILE *fp );
SplinePtr ReadSimpleSpline( FILE *fp, void *storage, int size );
void FreePSpline( PSplinePtr ps );
PSplinePtr ComputePSpline( Double *x, Double *y, int nPoints, void *storage, int bSpline );
void PSplineInterpolate( PSplinePtr ps, Double *x, Double *y, int points, void *storage );
Double PSplineArcLength( PSplinePtr ps );
void PrintPSpline( PSplinePtr ps, char *format, FILE *fp );

#ifdef __cplusplus
}
#endif

#endif /* __Spline__ */
