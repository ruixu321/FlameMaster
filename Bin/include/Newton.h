#ifndef NEWTON_H__
#define NEWTON_H__

#include "ArrayManager.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath> 
#include <cstring>

#include "Interrupt.h"
#if defined (applec) || defined (powerc)
#include <CursorCtl.h>
#endif

//		macros and types
#ifndef __MYIOCONST__
#define __MYIOCONST__
#define TAB "\t"
#define NEWL "\n"
#endif // __MYIOCONST__

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

// GB : New C++ Standard
using namespace std;
// GB : New C++ Standard

enum				{ kPath, kFileName };
enum RHSMode		{ kUpdate, kDoNothing, kSave, kRestore, kTest };	// mbo

typedef char Flag;

inline
int SIGN( int a )
{
	return ( a > 0 ) ? 1 : -1;
}

inline
Double SIGN( Double a )
{
	return ( a > 0 ) ? 1.0 : -1.0;
}

inline
Double MAX( Double a, Double b )
{
	return ( a > b ) ? a : b;
}

inline
Double MIN( Double a, Double b )
{
	return ( a < b ) ? a : b;
}

inline
int maxint( int a, int b )
{
	return ( a > b ) ? a : b;
}

inline
int minint( int a, int b )
{
	return ( a < b ) ? a : b;
}

//		class declaration
class NodeInfo;
class TBVPInput;
class TContinuation;
class TTime;
class TDamp;
class TGrid;
class TAdaptiveGrid;
class TNewton;
class TBVPSolver;

typedef NodeInfo *NodeInfoPtr;
typedef TBVPInput *TBVPInputPtr;
typedef TContinuation *TContinuationPtr;
typedef TTime *TTimePtr;
typedef TDamp *TDampPtr;
typedef TGrid *TGridPtr;
typedef TAdaptiveGrid *TAdaptiveGridPtr;
typedef TNewton *TNewtonPtr;
typedef TBVPSolver *TBVPSolverPtr;

typedef Double (*CalcFuncComponentPtr)( int equation, NodeInfoPtr nodeInfo, void *object, Flag theFlag );

// prototypes
void	Test_convergence( TBVPSolverPtr solver, void *object );
void	Monitor( TNewtonPtr bt, TContinuationPtr contin, TDampPtr damp, TTimePtr time );
Double	dfdyCentral( int i, int j, NodeInfoPtr nodeInfo );
Double	dfdyUpwind( int i, int j, CalcFuncComponentPtr func, NodeInfoPtr nodeInfo, void *object );
Double	dfdyUpwind( int i, int j, NodeInfoPtr nodeInfo );
void	DoExit( TNewtonPtr bt, void *object );
char 	*GetFullPath( const char *relPath, Flag type );
Double	Integrate( Double *x, Double *y, Double x1, Double x2, int nPoints, int &i1, int &i2 );
Double	Integrate( Double *x, Double *y, int i1, int i2 );
void	NewtonFatalError(const char *str);
void	CopyMatrix( MatrixPtr source, MatrixPtr dest );


int 		copy_mat( MatrixPtr ptr1, MatrixPtr ptr2);
int 		copy_vec( VectorPtr ptr1, VectorPtr ptr2);
int 		add_mat( MatrixPtr ptr1, MatrixPtr ptr2);
int 		sub_mat( int flag , MatrixPtr ptr1, MatrixPtr ptr2 );
int			mult_mat( MatrixPtr ptr, Double fact );
int			div_mat( int flag ,MatrixPtr ptr1, MatrixPtr ptr2 );
int			MaxOfFabsVec( int len, Double *vec, int offset );
Double 		Norm_mat( MatrixPtr ptr );
Double		MaxNorm_mat( MatrixPtr ptr, int *co, int *ro );
Double		MaxChange( MatrixPtr oldPtr, MatrixPtr newPtr, int *co, int *ro );
Double		ScaledNorm( MatrixPtr ptr, VectorPtr vecPtr );
Double		ScaledMaxNorm( MatrixPtr ptr, VectorPtr vecPtr );

// globals

extern AMPrintOptionsPtr gPrnt;	// initialized in TBVPSolver::Solve()
	
#ifdef powerc
typedef char **ConstStringArray;
#else
typedef const char *const *ConstStringArray;
#endif

typedef void ( *JacFuncPtr )( void *, NodeInfoPtr ); 
typedef void ( *RHSFuncPtr )( void *, NodeInfoPtr, RHSMode ); 
typedef void ( *OutputFuncPtr )( void *, FILE *, const char * ); 
typedef int  ( *PostIterFuncPtr )( void * ); 
typedef void ( *SetObjNodeInfoPtr )( int k, void *object ); 
typedef void ( *PostConvFuncPtr )( void * ); 
typedef void ( *UpdateBoundaryPtr )( void * ); 
typedef ConstStringArray ( *GetVariableNamesPtr )( void * ); 

typedef struct UTFuncs {
	UTFuncs( void );

	JacFuncPtr			JacobiFirstPoint;
	JacFuncPtr			JacobiRestPoint;
	JacFuncPtr			JacobiLastPoint;
	RHSFuncPtr			RHSFirstPoint;
	RHSFuncPtr			RHSRestPoint;
	RHSFuncPtr			RHSLastPoint;
	OutputFuncPtr		WriteOutput;
	PostIterFuncPtr		PostIter;
	SetObjNodeInfoPtr	SetObjNodeInfo;
	PostConvFuncPtr		PostConvergence;
	GetVariableNamesPtr	GetVariableNames;
	UpdateBoundaryPtr	UpdateLeftBoundary;
	UpdateBoundaryPtr	UpdateRightBoundary;
} UTFuncs, *UTFuncsPtr;

// classes

//*************************************************************************************************************
class NodeInfo {
public:
	NodeInfo( void ) { ; };

	int		nOfEquations;
	Double	parameter;
	Double	**a;
	Double	**b;
	Double	**c;
	Double	**bPrev;
	Double	**cNext;
	Double	h;
	Double	hm;
	Double	hNext;
	Double	hmPrev;
	Double	hnenn;
	Double	*yPrevPrev;
	Double	*yPrev;
	Double	*y;
	Double	*yNext;
	Double	*yNextNext;
	Double	*yMax;
	Double	*rhs;
	Double	*rhsSaved;
	int		*bcFlagLeft;
	int		*bcFlagRight;
	Double	*bcLeft;
	Double	*bcRight;
	Flag	firstPoint;
	Flag	lastPoint;
	Double	*x;
	int		gridPoint;
};

class TContinuation {
public:
	TContinuation( TNewtonPtr bt, TBVPInputPtr input ) { InitContin( bt, input ); };
	~TContinuation( void );

	int 	Contin_step( TBVPSolverPtr solver, TDampPtr damp );
	Double	GetParameter( void ) { return parameter; };
	void	AdjustContinDimension( int nOfPoints );
	
protected:
	Double		parameter;
	int			n_cont;
	int			cont_steps;
	Double		inc;
	VectorPtr	x_last_cont;
	MatrixPtr	y_last_cont;
	
private:
	void	InitContin( TNewtonPtr bt, TBVPInputPtr input );
};

class TBVPSolver {
public:
	TBVPSolver( TBVPInputPtr input );
	~TBVPSolver( void );

	TContinuationPtr	contin;		
	TTimePtr			time;		
	TDampPtr			damp;
	TNewtonPtr			bt;
	int		Solve( void *object );
	void	UpdateAllDimensions( int nOfPoints );
	void	ReInit( void );
};

class TGrid {
public:
	TGrid( TNewtonPtr bt, TAdaptiveGridPtr adapGrid, TBVPInputPtr input ) { InitGrid( bt, adapGrid, input ); };	
    ~TGrid( void );

	void		UpdateOneGridDimension( int number );	
	int			*GetBcFlagLeft( void ) { return bcFlagLeft; };
	int			*GetBcFlagRight( void ) { return bcFlagRight; };
	VectorPtr	GetBcLeft( void ) { return fBcLeft; };
	VectorPtr	GetBcRight( void ) { return fBcRight; };
	int			GetNGridPoints( void ) { return n_gridpoints; };
	Double		GetLeft( void ) { return left; };
	Double		GetRight( void ) { return right; };
	void 		SetLeft( Double leftB ) { left = leftB; };
	void 		SetRight( Double rightB ) { right = rightB; };
	MatrixPtr	GetY( void ) { return y; };
	VectorPtr	GetX( void ) { return x; };
	VectorPtr	GetYLeft( void ) { return yleft; };
	VectorPtr	GetYRight( void ) { return yright; };
	void 		AdjustNGridPoints( int number ) { n_gridpoints = number; };
 	void 		Make_equi_Grid( void );
	void		ComputeFDWeights( void );
   
private:
	void		InitGrid( TNewtonPtr bt, TAdaptiveGridPtr adapGrid, TBVPInputPtr input );
	//	void		SetInitialBC( TNewtonPtr bt, TInputDataPtr inp );

	int         n_gridpoints;	//  current number of gridpoints
    VectorPtr   x;				//	pointer to the locations of the gridpoints
    MatrixPtr   y;				//	independent variables ; offset for 1d-flame is: 	
								//				0: U
								//				1: V
								//				2: Z
								//	3 - nEquation: Y_{k}
								//			 last: T
																	
    Double      left, right;		// boundary points
	int			*bcFlagLeft;	//  uses enum BoundaryCondition
	int			*bcFlagRight;	//  uses enum BoundaryCondition
	VectorPtr	fBcLeft;		//  value of boundary condition 
	VectorPtr	fBcRight;		//  value of boundary condition 
    VectorPtr   yleft;          //  value at left boundary
    VectorPtr   yright;         //	value at right boundary
	VectorPtr	fFDWeightsConvM;	//  finite-difference weights for convection minus
	VectorPtr	fFDWeightsConv;		//  finite-difference weights for convection
	VectorPtr	fFDWeightsDiffM;	//  finite-difference weights for diffusion minus
	VectorPtr	fFDWeightsDiff;		//  finite-difference weights for diffusion
};

class TAdaptiveGrid {
public:
   	TAdaptiveGrid( TNewtonPtr bt, TBVPInputPtr input ) { InitGridding( bt, input ); };
	~TAdaptiveGrid( void );

	void		new_grid( TBVPSolverPtr solver, void *object );
	Double		GetLeft( void ) { return grid->GetLeft(); };
	Double		GetRight( void ) { return grid->GetRight(); };
	void 		SetLeft( Double leftB ) { grid->SetLeft( leftB ); };
	void 		SetRight( Double rightB ) { grid->SetRight( rightB ); };
	int			IsFine( void ) { return ( grid == fine ); };
	TGridPtr	GetFine( void ) { return fine; };	
	TGridPtr	GetCoarse( void ) { return coarse; };
	TGridPtr	GetCurrentGrid( void ) { return grid; };
	Flag		GetOneSolOneGrid( void ) { return fOneSolOneGrid; };
	Flag		GetSolHasNewGrid( void ) { return fSolHasNewGrid; };
	int			GetIterOfGeneration( void ) { return fIterOfGeneration; };
	void		InitIterOfGeneration( void ) { fIterOfGeneration = 0; };
	void		UnSetSolHasNewGrid( void ) { fSolHasNewGrid = FALSE; };
	void 		Make_coarse_grid( TNewtonPtr bt );
	VectorPtr	GetWorkVector( void ) { return aux; };
	MatrixPtr	GetWorkMatrix( void ) { return auxmat; };
	void 		UpdateGridDimension( void );
	void		SetSolutionScaler( void );
	VectorPtr	GetSolutionScaler( void ) { return fSolutionScaler; };
	void		ScaleSolution( TGridPtr grid, Double **y );
	void		DeScaleSolution( TGridPtr grid, Double **y );
	void		ToOldSolution( TBVPSolverPtr solver, VectorPtr xOld, MatrixPtr yOld );
	void		LinearInterpolateMat( TNewtonPtr bt, Double *x_old, Double **y_old, int n_old,
						Double *x, Double **y );

private:
	void 		InitGridding( TNewtonPtr bt, TBVPInputPtr input );	
	void		switch_to_coarse( TBVPSolverPtr solver );
	void		switch_to_fine( TBVPSolverPtr solver );
	void 		compute_M( TNewtonPtr bt );
	void 		compute_F( int n, Double *x, Double *M, TNewtonPtr bt );
	void		max3Filter( Double *x, Double *aux, int n );
	void 		max5Filter( Double *x, Double *aux, int n );
	void 		smoothFilter3( Double *x, Double *aux, int n );
	void 		filter( int which );
	void 		smooth_x( Double *x, Double *aux, int n, TNewtonPtr bt );
	void 		SplineInterpolateY( TNewtonPtr bt, Double *x_old, Double **y_old, int n_old,
						Double *x, Double **y );
	void 		LinearInterpolateY( TNewtonPtr bt, Double *x_old, Double **y_old, int n_old,
						Double *x, Double **y );

	VectorPtr	fSolutionScaler;	// build from the maximum values of the initial guess
	TGridPtr	grid;				// pointer to the current grid class, points to 
									// 'TGridPtr fine' or 'TGridPtr coarse'
	TGridPtr	fine;				// pointer to the fine grid class
	TGridPtr	coarse;				// pointer to the coarse grid class
    VectorPtr   M;					// pointer to the vector of error measures
    Double      R;					// R = (Æs)max/(Æs)min, a global constraint on the mesh sizes
    Double      q;					// measure of allowed error
    MatrixPtr   auxmat;				//
    VectorPtr   aux;				//
    VectorPtr   aux2;				// used by spline package
	Flag		init;				//
	Flag		fOneSolOneGrid;		//	if TRUE, a new grid is generated for every solution
	Flag		fSolHasNewGrid;		//	if fOneSolOneGrid is TRUE, fSolHasNewGrid is FALSE before
									//	and TRUE after the first grid generation
	int			fIterOfGeneration;	//	contains the number of the iteration of last grid generation
	Double		*integral;			//	internal workspace of length nOfEquations,
									//	used by function compute_M
	Double		fStart;				//	left boundary for grid correction
	Double		fEnd;				//	right boundary for grid correction
	Double		fAlpha;				//	fractional part of gridpoints in correction area
};

class TDamp {
public:
   	TDamp( TNewtonPtr bt, TBVPInputPtr input ) { InitDamp( bt, input ); };
	~TDamp( void );
	Double		GetLambda( void ) { return Lambda; };
	MatrixPtr	GetDyOld( void ) { return dyOld; };
	int 		Update_Lambda( TNewtonPtr bt, TTimePtr time, void *object );
	void		SetConvergence( void ) { converg_Lambda = TRUE; };
	void		UnSetConvergence( void ) { converg_Lambda = FALSE; };
	Flag		IsLambdaConvergent( void ) { return converg_Lambda; };
	void		AdjustDampDimension( int nOfPoints );
	
private:
	void 		InitDamp( TNewtonPtr bt, TBVPInputPtr input );
	void 		Predict_Lambda( TNewtonPtr bt, MatrixPtr dy, Double norm );
    Double      Lambda;
    Double      LambdaMin;
    Double      Sigma;
    Double      Tau;
    Flag        converg_Lambda;
/*  the value of yold is just in use during the function Update_Lambda
	and not from one call of the function to the next  */  
	MatrixPtr   yOld; 
    MatrixPtr   dyOld;  
};

class TNewton {
friend int TDamp::Update_Lambda( TNewtonPtr bt, TTimePtr time, void *object );
friend void TAdaptiveGrid::new_grid( TBVPSolverPtr solver, void *object );

public:
    TNewton( TBVPInputPtr input ) { InitNewton( input ); };
	~TNewton( void );

	int					Newton_step( TDampPtr damp, TTimePtr time, void *object );
	void 				BackSolve( MatrixPtr locDy, void *object, Flag updateRHS = TRUE );
	void 				AdjustNewtonDimensions( int nOfPoints );
	int 				UpdateRHS( MatrixPtr locDyMat, void *object );
	void				ScaleRHS( MatrixPtr dys );
	void				PrintTheVector( VectorPtr vecPtr, char *header );
	void				PrintVectorOfMatrix( MatrixPtr matPtr, int number, char *header );
	void				WriteOutput( void *object, FILE *fp, const char* tail );
	void				PrintSolution( Double *x, Double **y, ConstStringArray names, FILE *fp );
	ConstStringArray	GetVariableNames( void *object );

	UTFuncsPtr			GetUTFuncs( void ) { return fUtFuncs; };
	void				SetNOfEquations( int number );
	int					GetNOfEquations( void ) { return n_equations; };
//	void 				SetResNorm( MatrixPtr ptr );
//	void 				SetDyNorm( MatrixPtr ptr );
	void 				SetLeaveContin( void ) { leave_contin = TRUE; };
	void 				UnSetLeaveContin( void ) { leave_contin = FALSE; };
	void 				SetLeaveNewton( void ) { leave_newton = TRUE; };
	void 				UnSetLeaveNewton( void ) { leave_newton = FALSE; };
	void 				SetConvergeNewton( void ) { converg_newton = TRUE; };
	void 				UnSetConvergeNewton( void ) { converg_newton = FALSE; };
	void				SetNodeInfo( void *object, int k );
	void				SetParameterPtr( Double *parameterPtr ) { parameter = parameterPtr; };
	void 				SetLeft( Double leftB ) { grid->SetLeft( leftB ); };
	void 				SetRight( Double rightB ) { grid->SetRight( rightB ); };
	NodeInfoPtr			GetNodeInfo( void ) { return fNodeInfo; };
	MatrixPtr			GetDy( void ) { return dy; };
	MatrixPtr			GetDySaved( void ) { return dySaved; };
    TAdaptiveGridPtr	GetGrid( void ) { return grid; };
	void				ScaleMatrix( MatrixPtr matPtr, VectorPtr vecPtr );
	void				DeScaleMatrix( MatrixPtr matPtr, VectorPtr vecPtr );
	void				IncNIter( void ) { ++n_iter; };
	void				InitNIter( void ) { n_iter = 0; grid->InitIterOfGeneration(); };
	Double				GetLeft( void ) { return grid->GetLeft(); };
	Double				GetRight( void ) { return grid->GetRight(); };
	Double				GetResTol( void ) { return TolRes; };
	Double				GetDyTol( void ) { return TolDy; };
	Double				GetResNorm( void ) { return normRes; };
	Double				GetDyNorm( void ) { return normDy; };
	int 				GetConvergeNewton( void ) { return converg_newton; };
	int					GetCurrentGridPoints( void ) { return grid->GetCurrentGrid()->GetNGridPoints(); };
	Double				GetDx( void ) { return ( GetRight() - GetLeft() ) / ( GetCurrentGridPoints() + 1 ); };
	int					GetMaxGridPoints( void ) { return max_gridpoints; };
	int					GetInitialGridpoints( void ) { return initialGridPoints; };
	int					GetNEquations( void ) { return fA->rows; };
	int					GetNVariables( void ) { return fA->phys_rows; };
	int					GetMaxIter( void ) { return max_iter; };
	Flag				GetDampFlag( void ) { return damp_flag; };
	Flag				GetTimedepFlag( void ) { return timedep_flag; };
	Flag				GetContinuationFlag( void ) { return continuation_flag; };
	int					GetDeltaNewGrid( void ) { return delta_new_grid; };
	int					GetNIter( void ) { return n_iter; };
	int					GetConvergNewton( void ) { return converg_newton; };
	int					GetLeaveNewton( void ) { return leave_newton; };
	int					GetLeaveContin( void ) { return leave_contin; };
	Flag				IsModified( void ) { return fModified; };
	Flag				PrintFullRes( void ) { return fWriteFullRes; };
	Flag				IsGriddingStep( void ) { return fIsGriddingStep; };
	Flag				IsWriteBT( void ) { return fWriteBT; };
	Flag				IsWriteResiduum( void ) { return fWriteResiduum; };
	Flag				IsWatchGridding( void ) { return fWatchGridding; };
	Flag				IsWriteEverySolution( void ) { return fWriteEverySolution; };
	void				SetDampFlag( Flag what ) { damp_flag = what; };
	enum				FileType	{ kNone, kData, kText };
	char				*GetOutputPath( void ) { return fOutputPath; };
	char				*GetOutFileBuff( void ) { return fOutFile; };
	char 				*GetOutfileName( const char *name, FileType type, Flag iter = FALSE );
	FILE 				*GetOutfile( const char *name, FileType type, Flag iter = FALSE );
	void				PrintResiduals( MatrixPtr rhs, ConstStringArray names );
	Double				GetParameter( void ) { 
	   						return ( continuation_flag )? *parameter : 1.0; 
						};
	void				CheckGridding( TTimePtr tim );
	int					isNewGrid( void ) { return !( ( GetNIter() - 1 ) % GetDeltaNewGrid() || GetNIter() == 2 ); };
	void				PostConvergence( void *object ) { fUtFuncs->PostConvergence( object ); };
	int					CheckModified( void );
	void				SetUtFuncs( JacFuncPtr jFirst, JacFuncPtr jRest, JacFuncPtr jLast
									, RHSFuncPtr rFirst, RHSFuncPtr rRest, RHSFuncPtr rLast
									, OutputFuncPtr writeOutput, PostIterFuncPtr postIter
									, SetObjNodeInfoPtr setObjNodeInfo
									, PostConvFuncPtr postConvergence
									, GetVariableNamesPtr GetVariableNames
									, UpdateBoundaryPtr = NULL
									, UpdateBoundaryPtr = NULL );

	Flag                            IsFullStep ( void ) { return fullstep; }; //GB
	void                            ForceFullStep ( void ) { forcefullstep = true; }; //GB

private:
	void				InitNewton( TBVPInputPtr input );
	int 				UpdateJacRHS( void *object );
	int				UpdateJacRHSNum( void *object );
//	void 				Update_timedepJac( TTimePtr time );
//	void				Update_timedepRHS( TTimePtr time );
	void				ScaleSystem( void );
	void				ScaledY( Double **dY, Double **scaler );
	void				NewDeScaleSolution( Double **dy );
	int					Update_Y( TDampPtr damp, TTimePtr time );
	void				TestSolution( void );
	void 				SaveResNorm( void ) { normOldRes = normRes; };
	void 				SaveDyNorm( void ) { if ( n_iter ) normOldDy = normDy; };
	int					UpdateSolution( void *object, TDampPtr damp, TTimePtr time ); 
	void				InitModifiedConstants( void );
	Flag				AcceptModified( void );
	void				SetMaximaOfSolution( void );


    Flag         		timedep_flag;		// flag signalling time dependence
    Flag         		continuation_flag;	// if TRUE, we solve a continuation problem
    Flag         		damp_flag;			// if TRUE, use damped newton
    Flag         		fullstep;			// if TRUE, Newton step is modified
    bool                forcefullstep;		// if TRUE, a fullstep Newton is done (GB)
    int         		max_gridpoints;		// allowed number of gridpoints
    int         		n_equations;		// number of equations ( or number of independent
											// variables or number of blocks of the jacobian )
	int					initialGridPoints;	// initial number of gridpoints
    int         		max_iter;			// allowed number of iterations
	int					delta_new_grid;		// generate new grid all delta_new_grid iterations
    int         		converg_newton;		// 
    int         		leave_newton;		// 
    int         		n_iter;				// current iteration number
    int         		leave_contin;		//
	Double				*parameter;			// points to the member parameter of class TContinuation
    Double      		normOldRes;			//
    Double      		normOldDy;			//
	Flag				fWithModified;		// defined by user
    Flag      			fModified;			// evaluated in CheckModified
    int      			fModCounter;		// counts modified steps for the calculation of delta_t
    Flag      			fUseNumericalJac;	// numerical evaluation of the Jacobian if true
    Flag      			fUseSecOrdJac;	// numerical evaluation of the Jacobian if true
	Double				fAMod;
	Double				fBMod;
	Double				fCMod;
	Flag				fIsGriddingStep;	//
	Flag				fWriteFullRes;		//
	Flag				fWriteBT;
	Flag				fWriteResiduum;
	Flag				fWatchGridding;
	Flag				fWriteEverySolution;
	char				*fOutputPath;
	char				*fOutFile;
    Double      		normRes;			//
    Double      		normDy;				//
    Double      		TolRes;				//
    Double      		TolDy;				//
    Double      		dx;					//
    TAdaptiveGridPtr   	grid;				// 
    MatrixPtr   		dy;					//
    MatrixPtr   		dySaved;			//  used for numerical Jacobian and modified Newton
    TensorPtr   		fA, fB, fC;			//
    int         		*ip;				//
	MatrixPtr			fJacobiScaler;
	VectorPtr			fResiduals;			//	contains residuum for each variable
	NodeInfoPtr			fNodeInfo;
	UTFuncsPtr			fUtFuncs;

public:
void SetResNorm( MatrixPtr ptr ) { SaveResNorm(); normRes = Norm_mat( ptr ); };
void SetDyNorm( MatrixPtr ptr ) { SaveDyNorm(); normDy = Norm_mat( ptr ); };

	// for interrupt handling TNewton could use the global variable gExit, which is set in InterruptHandler
	// to initialize the InterruptHandler, function InstallInterruptHP is called in function InitNewton

	VectorPtr			fYMax;				//	contains maximum for each variable

};

class TTime {
public:
	TTime( TNewtonPtr bt, TBVPInputPtr input ) { InitTime( bt, input ); };
	~TTime( void );
	void		Calc_Delta_t( TBVPSolverPtr solver );
//	void		InitDeltaTMax( TNewtonPtr bt ) { Delta_tMax = 100.0 * bt->GetDx(); };
	void		InitDeltaTMax( TNewtonPtr /*bt*/ ) { Delta_tMax = 10.0; };
    Double		GetDeltaT( void ) { return Delta_t; };
    MatrixPtr	GetYOld( void ) { return yOld; };
    VectorPtr	GetXOld( void ) { return x_old; };
   	Double		GetDeltaTMax( void ) { return Delta_tMax; };
	void		AdjustTimeDimension( int nOfPoints );
	void		SetTimeConverged( void ) { fTimeConverged = TRUE; };
	void		UnSetTimeConverged( void ) { fTimeConverged = FALSE; };
	void		InitDeltaT( void ) { Delta_t = fDeltaTStart; };
	void		UnsetAdjustTime( void ) { adjustTime = FALSE; };
	void		SetAdjustTime( void ) { adjustTime = TRUE; };
	Flag		GetAdjustTime( void ) { return adjustTime; };
	Flag		GetTimeConverged( void ) { return fTimeConverged; };
	MatrixPtr	GetYPrime( void ) { return yPrime; };
	void		AdjustTimeStep( TBVPSolverPtr solver, Flag converged, void *object );
	Double		UpdateYPrime( TNewtonPtr bt );
	void		ReInitTimeMats( TNewtonPtr bt, void *object, Double dt, Double dtOld );
	void	 	InitTimeMats( TNewtonPtr bt, void *object );
	void 		InitYPrime( TNewtonPtr bt, void *object );
	int			GetTimeStep( void ) { return fNTimeStep; };
	void			InitNTimeStep( void ) { fNTimeStep = 1; };
protected:
    int         stepsGood;
    int         stepsBad;
	Flag		adjustTime;
	Flag		fIsDtGuess;
	Double		fAllowedFact;
    Double      Delta_t;
    Double      Delta_tOld;
    Double      Delta_tOldOld;
    Double      fDeltaTStart;
    Double      normGoUp;
    Double      normOld;
    Double      Delta_tMin;
    Double      Delta_tMax;
    MatrixPtr   yPrime;   
    MatrixPtr   yOld;   
	VectorPtr	x_old;
	Flag		fTimeConverged;
    int         fNTimeStep;
private:
	void		InitTime( TNewtonPtr bt, TBVPInputPtr input );
};

class TBVPInput {
public:
	TBVPInput( int nVariables ) { InitTBVPInput( nVariables ); };
	~TBVPInput( void );

	Flag		fWriteBT;
	Flag		fWriteResiduum;
	Flag		fWatchGridding;
	Flag		fWriteEverySolution;
	char		*fOutputPath;

	Flag		fWriteFullRes;
	Flag		fUseModifiedNewton;
	Flag		fUseNumericalJac;
	Flag      	fUseSecOrdJac;

	int			*bcFlagLeft;
	int			*bcFlagRight;
	VectorPtr	fBcLeft;
	VectorPtr	fBcRight;
    VectorPtr   yleft;
    VectorPtr   yright;

	int			fNVariables;
	int			fInitialEquations;
	int			fMaxGridPoints;
	int			fInitialGridPoints;
	Flag		fDampFlag;
	Flag		fTimeDepFlag;
	Flag		fContinFlag;
	int			fDeltaNewGrid;
	Double		fTolRes;
    Double		fTolDy;
    int			fMaxIter;
    Double		fLeft;
    Double		fRight;
	Flag		fOneSolOneGrid;
    Double		fR;
    Double		fQ;
	Double		fStart;
	Double		fEnd;
	Double		fAlpha;
    Double		fLambdaMin;
    int			fContSteps;
    Double		fDeltaTStart;
	
private:
	void			InitTBVPInput( int nVariables );
};

#endif // NEWTON_H__
