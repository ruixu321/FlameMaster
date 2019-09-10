/*
	Newton.h
	20-01-92	( pt )
*/

#ifndef __Newton__
#define __Newton__


#ifdef applec
#pragma once
#endif

struct NewtonInfo;

typedef int (*NewtonFunc)( const VectorPtr vec, VectorPtr f, void *object );
typedef int (*NewtonDeriv)( const VectorPtr vec, MatrixPtr m, struct NewtonInfo *ni, void *object );
typedef Double (*NewtonNorm)( VectorPtr x );


typedef struct NewtonInfo {
	int				n;				/* number of equations					*/
	int				step;			/* current step							*/
	int				maxSteps;		/* maximum number of steps				*/
	int				modified;		/* signals modified newton's method		*/
	int				converged;		/* signals convergence					*/
	Double			lastNorm;		/* last norm evaluated					*/
	Double			tol;			/* required tolerance					*/
	Double			minScale;		/* minimum value for scaling			*/
	IntVectorPtr	index;			/* integer workspace					*/
	VectorPtr		x;				/* solution vector						*/
	int				freeX;
	VectorPtr		inc;			/* solution increment vector			*/
	VectorPtr		f;				/* rhs vector							*/
	VectorPtr		ws;				/* workspace (dimensioned n)			*/
	VectorPtr		scale;			/* scaling vector						*/
	VectorPtr		relErr;			/* relative error vector				*/
	MatrixPtr		jac;			/* jacobian								*/
	NewtonFunc		func;			/* function to evaluate the rhs			*/
	NewtonDeriv		deriv;			/* function to evaluate the jacobian	*/
	NewtonNorm		norm;			/* function to evaluate the norm		*/
} NewtonInfo, *NewtonInfoPtr;


#ifdef __cplusplus
extern "C" {
#endif
NewtonInfoPtr NewNewtonInfo( int n, VectorPtr xIn );
void FreeNewtonInfo( NewtonInfoPtr ni );
void SetNewtonFuncs( NewtonInfoPtr ni, int (*func)( const VectorPtr vec, VectorPtr f, void *object ), NewtonDeriv deriv, 
	NewtonNorm norm );
void ResetNewtonScale( NewtonInfoPtr ni, Double val );
void SetNewtonScale( NewtonInfoPtr ni, VectorPtr scale );
int NewtonSolve( NewtonInfoPtr ni, int resetScale, int newJacobian, void *object );
void PrintNewtonInfo( NewtonInfoPtr ni, FILE *fp, const char *label );

void warning_fl( const char *errorString, const char *file, int line );
#define WARNING( s ) warning_fl( (s), __FILE__, __LINE__ )
#ifdef __cplusplus
}
#endif


#endif	/* __Newton__ */
