/*
	ArrayManager.h: header file for ArrayManager routines
	
	© Copyright 1989-1993 by Josef Goettgens. All rights reserved.
*/

#ifndef __ArrayManager__
#define __ArrayManager__

#ifndef getc
#include <stdio.h>
#endif

/* we always use the custom data type Double */
typedef double Double;
#ifndef __alligator__
# include "alligator.h"
#endif

/*	Default definition of FPUType	*/
#define FPUType	Double

#ifdef THINK_C
# include <Memory.h>
# undef		FPUType
# define	FPUType	long double
# define	COPYSIGN
#endif

#ifdef applec
# ifndef __MEMORY__
#  include	<Memory.h>
# endif
# include	<SANE.h>	/* copysign() is declared here */
# undef		FPUType
# define	FPUType	long double
#endif

/*	define a special symbol for Symantec's MPW compiler
#if defined(__SC__) && defined(macintosh)
#define MPWSC 1
#endif
#ifdef MPWSC
# undef		FPUType
# define	FPUType	long double
#endif
*/

#ifdef HP
# define COPYSIGN
#endif


#ifndef FALSE
#define FALSE	0
#endif
#ifndef TRUE
#define TRUE	1
#endif


typedef struct Vector {
	int		len;			/* logical length of the vector		*/
	int		phys_len;		/* physical length of the vector	*/
	Double	*vec;			/* pointer to the elements			*/
} Vector;
typedef Vector *VectorPtr;

typedef struct IntVector {
	int		len;			/* logical length of the vector		*/
	int		phys_len;		/* physical length of the vector	*/
	int		*vec;			/* pointer to the elements			*/
} IntVector;
typedef IntVector *IntVectorPtr;

typedef struct Matrix {
	int		rows, cols;				/* number of rows and columns	*/
	int		phys_rows, phys_cols;	/* physical dimensions			*/
	int		partition;				/* partition type				*/
	Double	**mat;					/* pointer to matrix-elements	*/
} Matrix;
typedef Matrix *MatrixPtr;

typedef struct IntMatrix {
	int		rows, cols;				/* number of rows and columns	*/
	int		phys_rows, phys_cols;	/* physical dimensions			*/
	int		partition;				/* partition type				*/
	int		**mat;					/* pointer to matrix-elements	*/
} IntMatrix;
typedef IntMatrix *IntMatrixPtr;

enum { kRowPointers, kColumnPointers };


typedef struct Tensor {		/* That's our term for a 3D array.	*/
	int		rows, cols;				/* logical dimensions			*/
	int		phys_rows, phys_cols;	/* physical dimensions			*/
	int		partition;				/* partition type				*/
	Double	**mat;					/* for handling sub-matrices	*/
	int		planes, phys_planes;	/* number of planes				*/
	Double	***tensor;				/* pointer to tensor elements	*/
} Tensor;
typedef Tensor *TensorPtr;

enum {
	kLUDecomposition = 1,	/* performs lu decomposition destructively	*/
	kBackSubstitution,		/* performs back substitution				*/
	kSolve					/* performs options 1 & 2					*/
};


typedef struct AMPrintOptions {
	int		lineWidth;		/* max. # of chars per line; default is 80	*/
	const char	*format;		/* format string for numbers; default is %g	*/
	const char	*intFormat;		/* format string for ints; default is %d	*/
	const char	*sep;			/* separator string; default is " " (space)	*/
	const char	*title;			/* title string; default is NULL			*/
	const char	*colLabel;		/* column label; default is NULL			*/
	const char	*rowLabel;		/* row label; default is "Row #%d"			*/
	const char	*planeLabel;	/* plane label; default is "Plane #%d"		*/
} AMPrintOptions;
typedef AMPrintOptions *AMPrintOptionsPtr;


/*
	Prototypes
*/
#ifdef __cplusplus
extern "C" {
#endif

/*	Allocation/Deallocation functions
*/

extern Double *_New1DArray( int );
extern Double **_New2DArray( int, int );
extern Double **_Make2DArray( Double *, int, int );
extern Double ***_New3DArray( int, int, int );
extern Double ***_Make3DArray( Double *, int, int, int );
extern void _Free1DArray( Double * );
extern void _Free2DArray( Double ** );
extern void _Free3DArray( Double *** );
extern void _Release2DArray( Double ** );
extern void _Release3DArray( Double *** );
extern VectorPtr _NewVector( int );
extern void _DisposeVector( VectorPtr );
extern MatrixPtr _NewMatrix( int, int, int );
extern void _DisposeMatrix( MatrixPtr );
extern TensorPtr _NewTensor( int, int, int, int );
extern void _DisposeTensor( TensorPtr );
extern int *_New1DIntArray( int );
extern int **_New2DIntArray( int, int );
extern void _Free1DIntArray( int * );
extern void _Free2DIntArray( int ** );
extern IntVectorPtr _NewIntVector( int );
extern void _DisposeIntVector( IntVectorPtr );
extern IntMatrixPtr _NewIntMatrix( int, int, int );
extern void _DisposeIntMatrix( IntMatrixPtr );

extern Double *New1DArray_dbg( int, const char *, int );
extern Double **New2DArray_dbg( int, int, const char *, int );
extern Double **Make2DArray_dbg( Double *, int, int, const char *, int );
extern Double ***New3DArray_dbg( int, int, int, const char *, int );
extern Double ***Make3DArray_dbg( Double *, int, int, int, const char *, int );
extern void Free1DArray_dbg( Double *, const char *, int );
extern void Free2DArray_dbg( Double **, const char *, int );
extern void Free3DArray_dbg( Double ***, const char *, int );
extern void Release2DArray_dbg( Double **, const char *, int );
extern void Release3DArray_dbg( Double ***, const char *, int );
extern VectorPtr NewVector_dbg( int, const char *, int );
extern void DisposeVector_dbg( VectorPtr, const char *, int );
extern MatrixPtr NewMatrix_dbg( int, int, int, const char *, int );
extern void DisposeMatrix_dbg( MatrixPtr, const char *, int );
extern TensorPtr NewTensor_dbg( int, int, int, int, const char *, int );
extern void DisposeTensor_dbg( TensorPtr, const char *, int );
extern int *New1DIntArray_dbg( int, const char *, int );
extern int **New2DIntArray_dbg( int, int, const char *, int );
extern void Free1DIntArray_dbg( int *, const char *, int );
extern void Free2DIntArray_dbg( int **, const char *, int );
extern IntVectorPtr NewIntVector_dbg( int, const char *, int );
extern void DisposeIntVector_dbg( IntVectorPtr, const char *, int );
extern IntMatrixPtr NewIntMatrix_dbg( int, int, int, const char *, int );
extern void DisposeIntMatrix_dbg( IntMatrixPtr, const char *, int );


#ifdef qMemDebug

#define New1DArray(n)			New1DArray_dbg((n),__FILE__,__LINE__)
#define New2DArray(n,m)			New2DArray_dbg((n),(m),__FILE__,__LINE__)
#define Make2DArray(a,n,m)		Make2DArray_dbg((a),(n),(m),__FILE__,__LINE__)
#define New3DArray(n,m,l)		New3DArray_dbg((n),(m),(l),__FILE__,__LINE__)
#define Make3DArray(a,n,m,l)	Make3DArray_dbg((a),(n),(m),(l),__FILE__,__LINE__)
#define Free1DArray(a)			Free1DArray_dbg((a),__FILE__,__LINE__)
#define Free2DArray(a)			Free2DArray_dbg((a),__FILE__,__LINE__)
#define Free3DArray(a)			Free3DArray_dbg((a),__FILE__,__LINE__)
#define Release2DArray(a)		Release2DArray_dbg((a),__FILE__,__LINE__)
#define Release3DArray(a)		Release3DArray_dbg((a),__FILE__,__LINE__)
#define NewVector(n)			NewVector_dbg((n),__FILE__,__LINE__)
#define DisposeVector(a)		DisposeVector_dbg((a),__FILE__,__LINE__)
#define NewMatrix(n,m,l)		NewMatrix_dbg((n),(m),(l),__FILE__,__LINE__)
#define DisposeMatrix(a)		DisposeMatrix_dbg((a),__FILE__,__LINE__)
#define NewTensor(n,m,l,o)		NewTensor_dbg((n),(m),(l),(o),__FILE__,__LINE__)
#define DisposeTensor(a)		DisposeTensor_dbg((a),__FILE__,__LINE__)
#define New1DIntArray(n)		New1DIntArray_dbg((n),__FILE__,__LINE__)
#define New2DIntArray(n,m)		New2DIntArray_dbg((n),(m),__FILE__,__LINE__)
#define Free1DIntArray(a)		Free1DIntArray_dbg((a),__FILE__,__LINE__)
#define Free2DIntArray(a)		Free2DIntArray_dbg((a),__FILE__,__LINE__)
#define NewIntVector(n)			NewIntVector_dbg((n),__FILE__,__LINE__)
#define DisposeIntVector(a)		DisposeIntVector_dbg(a,__FILE__,__LINE__)

#else

#define New1DArray(n)			_New1DArray((n))
#define New2DArray(n,m)			_New2DArray((n),(m))
#define Make2DArray(a,n,m)		_Make2DArray((a),(n),(m))
#define New3DArray(n,m,l)		_New3DArray((n),(m),(l))
#define Make3DArray(a,n,m,l)	_Make3DArray((a),(n),(m),(l))
#define Free1DArray(a)			_Free1DArray((a))
#define Free2DArray(a)			_Free2DArray((a))
#define Free3DArray(a)			_Free3DArray((a))
#define Release2DArray(a)		_Release2DArray((a))
#define Release3DArray(a)		_Release3DArray((a))
#define NewVector(n)			_NewVector((n))
#define DisposeVector(a)		_DisposeVector((a))
#define NewMatrix(n,m,l)		_NewMatrix((n),(m),(l))
#define DisposeMatrix(a)		_DisposeMatrix((a))
#define NewTensor(n,m,l,o)		_NewTensor((n),(m),(l),(o))
#define DisposeTensor(a)		_DisposeTensor((a))
#define New1DIntArray(n)		_New1DIntArray((n))
#define New2DIntArray(n,m)		_New2DIntArray((n),(m))
#define Free1DIntArray(a)		_Free1DIntArray((a))
#define Free2DIntArray(a)		_Free2DIntArray((a))
#define NewIntVector(n)			_NewIntVector((n))
#define DisposeIntVector(a)		_DisposeIntVector(a)

#endif	/* qMemDebug */


void MinMax( const Double *x, int n, Double *maxi, Double *mini );
void ClearArray( Double *ptr, int len );
void ClearIntArray( int *ptr, int len );

void Clear1DArray( Double *v, int len );
void Clear2DArray( Double **m, int rows, int cols );
void Clear3DArray( Double ***t, int planes, int rows, int cols );
void ClearVector( VectorPtr v );
void ClearMatrix( MatrixPtr m );
void ClearTensor( TensorPtr t );

/*	Access functions
*/
MatrixPtr SubMatrix( int plane, TensorPtr t );
int GetMElement( Double *value, int row, int column, MatrixPtr m );
int PutMElement( Double value, int row, int column, MatrixPtr m );
int GetTElement( Double *value, int plane, int row, int column, TensorPtr t );
int PutTElement( Double value, int plane, int row, int column, TensorPtr t );

enum {
	kIndexOK,
	kInvalidPartition,
	kNoData,
	kInvalidLogicalRowIndex,
	kInvalidLogicalColIndex,
	kInvalidLogicalPlaneIndex,
	kInvalidPhysicalRowIndex,
	kInvalidPhysicalColIndex,
	kInvalidPhysicalPlaneIndex
};


/*	Printing functions
*/
void AccessError( int err, int plane, int row, int col, char *label );
char *AMInfo( char * );

AMPrintOptionsPtr NewAMPrintOptions( void );
void FreeAMPrintOptions( AMPrintOptionsPtr prnt );
void DefaultAMPOpts( AMPrintOptionsPtr prnt );

void Print1DArray( Double *v, int len, AMPrintOptionsPtr prnt, FILE *fp );
void Print2DArray( Double **m, int rows, int cols, int partition, 
				   AMPrintOptionsPtr prnt, FILE *fp );
void Print3DArray( Double ***t, int planes, int rows, int cols,
				   int partition, AMPrintOptionsPtr prnt, FILE *fp );
void PrintVector( VectorPtr v, AMPrintOptionsPtr prnt, FILE *fp );
void PrintMatrix( MatrixPtr m, AMPrintOptionsPtr prnt, FILE *fp );
void PrintTensor( TensorPtr t, AMPrintOptionsPtr prnt, FILE *fp );

void Print1DIntArray( int *v, int len, AMPrintOptionsPtr prnt, FILE *fp );
void Print2DIntArray( int **m, int rows, int cols, int partition, 
				   AMPrintOptionsPtr prnt, FILE *fp );
void PrintIntVector( IntVectorPtr v, AMPrintOptionsPtr prnt, FILE *fp );
void PrintIntMatrix( IntMatrixPtr m, AMPrintOptionsPtr prnt, FILE *fp );

/*
	Basic Linear Algebra stuff
*/
void ludcmp( Double **a, int n, int *indx, Double *d, Double *vv );
void lubksb( Double **a, int n, const int *indx, Double *b );
void thomas( Double *a, Double *b, Double *c, Double *r, int n );

void DoolittleC( Double *r, Double **b, int nrows, int ncols, int ihalfb, int kkk );
void DoolittleR( Double *r, Double **b, int nrows, int ncols, int ihalfb, int kkk );

/* AM_blas.c */
int iamax(int n, const Double *sx, int incx);
Double asum(int n, const Double *sx, int incx);
void saxpy(int n, Double sa, Double *sx, int incx, Double *sy, int incy);
void saxpyx(int n, Double sa, Double *sx, int incx, Double *sy, int incy);
void copy(int n, Double *sx, int incx, Double *sy, int incy);
Double dot(int n, const Double *sx, int incx, const Double *sy, int incy);
Double nrm2(int n, const Double *sx, int incx);
Double norm( Double *x, int n );
Double r1mach(void);
void scal(int n, Double sa, Double *sx, int incx);
void vexopy(int n, Double *v, Double *x, Double *y, int itype);
void vfill(int n, Double *v, Double val);
#ifdef COPYSIGN
Double copysign(Double x, Double y);
#endif

/* AM_Gauss.c */
int gefa( MatrixPtr a, int *ipvt );
int gesl( MatrixPtr a, int *ipvt, Double b[], int job );
int geco( MatrixPtr a, int ipvt[], Double *rcond, Double z[] );

/* AM_BlockTriDiagSolver.c */
int decbt( TensorPtr a, TensorPtr b, TensorPtr c, int *ip );
int solbt( TensorPtr a, TensorPtr b, TensorPtr c, MatrixPtr y, int *ip );

#ifdef __cplusplus
}
#endif

#endif	/* __ArrayManager__ */
