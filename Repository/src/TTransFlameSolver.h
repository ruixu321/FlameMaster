#ifndef __TTRANS_FLAME_SOLVER_H__
#define __TTRANS_FLAME_SOLVER_H__

#include "ConfigFM.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>

#if SUNDIALS_VERSION_MAJOR > 2
#include <sunmatrix/sunmatrix_band.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_band.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>       /* access to CVDls interface            */
#else
#include <cvode/cvode_band.h>
#endif // SUNDIALS_VERSION_MAJOR > 2

#include <stddef.h> // for ptrdiff_t
#include <time.h>

#include "Constants.h"
#include "MapMan.h"
#include "Interrupt.h"
#include "BetaPDF.h"

#define CVODE
#define CBETAPDF

#ifndef CVODE
#include "dassl.h"
#endif

template<typename Flame>
void ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
template<typename Flame>
void ResTransFlameImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
				, int *iRes, Double *rPar, int *iPar );
template<typename Flame>
void JacTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *pd
			, Double *cj, Double *rPar, int *iPar );
template<typename Flame>
int ResTransFlameImpliSolverCV( realtype T, N_Vector u, N_Vector udot, void *data );

template<class Species>
class TTransFlameSolver : public T0DFlame<Species> {
template<typename Flame>
friend void ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
template<typename Flame>
friend void ResTransFlameImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
				, int *iRes, Double *rPar, int *iPar );
template<typename Flame>
friend void JacTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *pd
			, Double *cj, Double *rPar, int *iPar );
template<typename Flame>
friend int ResTransFlameImpliSolverCV( realtype T, N_Vector u, N_Vector udot, void *data );
public:
  using T0DFlame<Species>::fInputData;
  using T0DFlame<Species>::fSpecies;
  using T0DFlame<Species>::fReaction;
  using T0DFlame<Species>::fSoot;
  using T0DFlame<Species>::fProperties;
  using T0DFlame<Species>::GetProperties;
  using T0DFlame<Species>::GetOutfile;
  using T0DFlame<Species>::GetFuelIndexVec;
  using T0DFlame<Species>::GetSpecies;
  using T0DFlame<Species>::GetAuthor;
  using T0DFlame<Species>::GetElementMassFraction;
  using T0DFlame<Species>::GetNFuels;
  using T0DFlame<Species>::GetOutfileName;
  using T0DFlame<Species>::GetInputData;
  using T0DFlame<Species>::GetNu;
  using T0DFlame<Species>::GetSoot;
  using T0DFlame<Species>::ComputeZBilger;
  using T0DFlame<Species>::CompLewisNumbers;
  using T0DFlame<Species>::GetFuelIndex;
  using T0DFlame<Species>::PrintFlameletVector;

	TTransFlameSolver( const FirstInput& firstInp ) : 
            fGrid(0),
			fFirstSpecies(1), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
			fVariables(fInputData->GetCounter()->species+fInputData->fVariablesWithoutSpecies+1),
			fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies+1),
			fSootMoments(fTemperature+1),
			fNOfEquations(fInputData->GetCounter()->species+fVariablesWithoutSpecies - fInputData->GetCounter()->steadyStates),
			T0DFlame<Species>( firstInp ) { InitTTransFlameSolver();};	
	virtual ~TTransFlameSolver( void );

	int			GetOffsetFirstSpecies( void ) { return fFirstSpecies; };
	int			GetOffsetTemperature( void ) { return fTemperature; };
	int			GetVariablesWithoutSpecies( void ) { return fVariablesWithoutSpecies; };
	ConstStringArray	GetVariableNames( void ) { return fVariableNames; };

	void		Initialize( Double timeStart
					, ConstStringArray names
					, Double **startSolution, int vars
					, Double *grid, int gridPointsA
					, Double pressureStart, Double scalarDissRateStart
					, Double firstTimeStep
					, Double ZRef
					, Double ZRStart );
	void		GetSpeciesData( char **newNames, MatrixPtr highMat
					, MatrixPtr lowMat, Double *molarMass );
	void		GetSolution( ConstStringArray names, Double **outSol, Double *grid
					, int gridPointsA, int vars, Double *density = NULL );
	void		GetSolutionInfo( Double *timeStep );
	void		SetSootSources( int whichMoment, MatrixPtr sourcesMat
					, Double theTime, Double *temp, Double **Y
					, Double pressure, Double **moments );
	int			GetSootSources( int whichMoment, ConstStringArray names, Double **sources, Double *grid, int gridPointsA, int vars );
	void		MakeGrid( VectorPtr theGrid, Double left, Double right, Flag equidistant );
	void		MakeGrid( Double *grid, int gridPoints
							, Double left, Double right, Flag equidistant );
	void		PrintSolution( FILE *fp, int nGPoints, int nEq, Double *grid
									, Double **sol, ConstStringArray names );
	Flag		Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, int deltaStepsOut );
	Flag		Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd, int deltaStepsOut );
	void		SetEndValues( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd );
	void		PrintDissRate( Double t );
	void		WriteFlameletFile( FILE *fp, const char *head, const char *tail );
	void		WriteFlameletFileInitialize( FILE *fp, const char *head, const char *tail );
	void		WriteFlameletFile( Double time, Double **sol, Double *x, 
											FILE *fp, const char *head, const char *tail );
	FILE 		*GetOutputFile( Double time, const char *head, const char *tail, FileType type );
	FILE		*GetOutputFile( const char *head, const char *tail, const FileType type );
	Double		GetCurrentTime( void );
	Double 		GetDissRateReq( Double ZReq, Double dissRate, Double z );
	void 		PrintProdRateGlobalReac( Double time );
	Double		GetRandomNumber( void ) { return fRandomNumber; };
	Double		GetDissFunc( Double z );
	void		ReadStartProfiles( TInputDataPtr inp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure
				, Double *chi, Double *theTime, Double *currTimeStep
				, int tempOff, int speciesOff, int sootOff );
	Double		GetTotEnt( Double temp, Double *Y );

protected:
	Double	GetTempOfEnthalpy( Double ent, Double *Y, Double initialGuess );
	Double	Interpol( Double t, Double valOld, Double tOld, Double valNew, Double tNew );
	Double	ComputeMeanFromPDF( Double t, int speciesIndex, Double **PDF, Double *timePDF );
	Double	CheckPDF( Double t, Double **PDF, Double *timePDF );
	Double	ComputeEmissionIndex( Double time, int speciesIndex, Double **PDF
								, Double *timePDF );
	Double	ComputeEmissionIndexSoot( Double time, int which, Double **PDF
								, Double *timePDF );
	Double	TurbMeanX( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar = NULL );
	Double	TurbMeanY( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar = NULL );
	Double	TurbMeanTemp( Double t, Double ZMean, Double ZVar, Double *tempVar = NULL );
	Double	TurbMeanTotEnergy( Double t, Double ZMean, Double ZVar );
	Double	TurbMeanZBarlow( Double t, Double ZMean, Double ZVar );
	Double	TurbMeanSoot( Double t, int sootIndex, Double ZMean, Double ZVar, Double *sootVar );

	const int	fGrid;
	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fVariablesWithoutSpecies;
	const int	fNOfEquations;
	const int	fVariables;
	const int	fSootMoments;

private:
  	void	InitTTransFlameSolver( void );
	void	SetInitialConditions( Double *y, TInputDataPtr inp );
	void 	InitDassl( int nGridPoints );
	void 	InitDasslWorkSpace( int i );
	void 	SetInfo( int i );
	Flag	FirstStep( void );
	Flag	OneStep( int k );
	int		GetActualPoint( Double tEnd );
	void	ComputeGPDistributionWeights( Double *y );
	void	SetFDWeights( Double *y );
	void	SaveSolution( int k, Double time, Double *y );
	void	SetWeights( void );
	void	SetMaxVals( void );
	void	SetInitial( ConstStringArray names, Double **startSol, Double *grid
						, int gridPoints, int vars );
	Double	GetTotEnt( int k, Double *Z, Double t );
	Double	GetZStoi( void );
	Double	GetDissRate( Double time, Double z );
	Double	GetDissRate( Double t, Double z, Double rho );
	Double	ExactChi( Double Z );
	Double	DissRateFact( Double t, Double rho );
	Double	GetRefDissRate( Double time );
	Double	GetZR( void ) { return Interpol( GetCurrentTime(), fZRStart, fTStart, fZREnd, fTEnd ); };
	Double	GetDelQ( Double t, Double Z );
	void	SetMolarMassOverRInf( void );
	void 	ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
						, int *iRes, Double *rPar, int *iPar );
	void	ResTransFlameImpliSolver( Double * /*T*/, Double *y, Double *yPrime
				, Double *delta
				, int * /*iRes*/, Double * /*rPar*/, int * /*iPar*/ );
	void	JacTransFlameSolver( Double *T, Double *y, Double *yPrime
			, Double **pd, Double cj, Double *rPar, int *iPar );
	int		GetVariableIndex( const char *name );
	int		GetVariableIndex( const char *name, ConstStringArray array, int len );
	void	SetMaxTimeStep( int kAct, Double t );
	Double	GetTempSource( Double Z );
	Double	GetPressure( void ) 
			{ std::cerr << "#error: function GetPressure not allowed in class " << 
					"TTransFlameSolver" << NEWL; return -1.0; };
	Double	GetPressure( Double time );
	void	PrintProdRate( Double time, int speciesIndex, FILE *fp );
	void	WriteSootInfo( Double theTime, Double *temp, Double **Y
				, Double pressure, Double **moments, FILE *fp );

	void	DoExit( void );
	void	SetOutSolution( void );
	Flag	FirstImpliStep( void );
	Flag	OneImpliStep( void );
	void	SetRandomNumber( void );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure, Double *chi
				, Double *theTime, Double *currTimeStep
				, int tempOff, int speciesOff, int sootOff );
	void	SetInitialValues( int nGridPoints, int nVars, Double *x, Double **y );
	void	PostIter( Double t );
    Double  GetRosseRadiation( int k, Double *nTemp, Double **moments, Double rfPrev
                                            , Double rfCurr, Double rfNext );


	void	UpdateThermoProps( int k, Double *Y, Double temp, Double &pressure
								, Double &density, EqOfState what, Double *sootMoments );

	void	ComputeDeltaI( int k, Double *Y, Double temp );
	Double	CompOneDeltaI( int i, int k, Double *Y );
	void	CompDeltaIG( int k, Double temp );
	void	SetGrid( void );
	Double	GetU( Double t, Double Z );
	Flag	CheckOutput( Double *t );
	int		InitCVODE( void );
	void	SolutionToCVode( void );
	void	SolutionFromCVode( void );

	enum FortranInd	{ kF1, kF2, kF3, kF4, kF5, kF6, kF7, kF8, kF9, kF10, kF11
					, kF12, kF13, kF14, kF15, kF16, kF17, kF18, kF19, kF20 };

	// input

	string      fRadiationName;

	Double		fTStart;
	Double		fTEnd;
	Double		fZl;	// left boundary
	Double		fZr;	// right boundary
	Double		fFirstTimeStep;
	Double		*fRTol;
	Double		*fATol;
	int			fNOfSolver;
	int			fNGridPoints;
	Flag		fEquidistant;
	Flag		fPrintMolarFractions;
	Double		fDeltaTMax;
	VectorPtr	fZRin;
	VectorPtr	fTimeIn;
	VectorPtr	fZIn;
	MatrixPtr	fChiIn;
	MatrixPtr	fUstOverU;
	MatrixPtr	fDelQ;
	VectorPtr	fZCount;
	VectorPtr	fChiCount;
	

	Flag		fFirstCall;
	Double		fZRef;
	Double		fWOverRInf;

	FILE		*fOutFilePtr;
	FILE		*fCAinFile;

	Double		fPressStart;
	Double		fChiStart;
	Double		fTempOxStart;
	Double		fTempFuelStart;
	Double		fZRStart;
	Double		fPressEnd;
	Double		fChiEnd;
	Double		fTempOxEnd;
	Double		fTempFuelEnd;
	Double		fZREnd;
	Double		fRandomNumber;
	Double		fDTdtOx;
	Double		fDTdtFu;
	Double		fDPdt;
	Double		fdDeltaZdt;
	VectorPtr	fTotEntStart;
	VectorPtr	fTotEntEnd;

	MatrixPtr	fSpecHeatCp;
	MatrixPtr	fSpecEnthalpy;
	MatrixPtr	fProdRate;

	MatrixPtr	fSpecSource;   // Rui

	VectorPtr	fHeatCpMix;
	VectorPtr	fViscosity;
	VectorPtr	fMolarMassMix;
	VectorPtr	fDensity;
	VectorPtr	fLambdaOverCpMix;

	VectorPtr	fDiffTermY;
	VectorPtr	fDiffTermW;
	VectorPtr	fDiffCorrY;
	VectorPtr	fDiffCorrW;
	
	int			fActualPoint;
	int			fMaxOrd;
	char		**fVariableNames;
	int			**fNActualOrd;
	int			**fNActualStep;
	Double		*fActualTimeStepSize;
	Double		fTout;
	Double		*fh;
	Double		*fhm;
	Double		*fMonFct;
	Double      fKappa;
	Double		fTauGrid;
	Double      fRadLossPercent;   // Rui
	VectorPtr	fgpDens; 
	VectorPtr	fSolGrid;
	Double		*fFDWCurr;
	Double		*fFDWPlus;
	Double		*fFDWMinus;
	Double		*fWCurr;
	Double		*fWPlus;
	Double		*fWMinus;
	VectorPtr	fSolTime;
	MatrixPtr	fSolMassFracs;
	MatrixPtr	fSolSootMoments;
	VectorPtr	fSolTemp;
	VectorPtr	fSolOldTime;
	MatrixPtr	fSolOldMassFracs;
	MatrixPtr	fSolOldSootMoments;
	VectorPtr	fSolOldTemp;
	VectorPtr	fMaxVals;

	MatrixPtr	fMassFracsWork;
	MatrixPtr	fSootMomentsWork;
	VectorPtr	fTempWork;
	VectorPtr	fOutSolWork;
	Double		*fTempGSave;
	Double		**fYGSave;
	Double		**fDeltaI;
	Double		***fG_ij;

// cvode stuff
	int			fNActualOrdCV;
	int			fNActualStepCV;
	Double		fActualTimeStepSizeCV;
	void		*fMem;

#if SUNDIALS_VERSION_MAJOR > 2
    SUNMatrix A;
    SUNLinearSolver LS;
#endif // SUNDIALS_VERSION_MAJOR > 2

	N_Vector	fCVY;       // solution vector space
	realtype	*fCYdata;  // pointer to access data in N_Vector

//	variables used for dassl
	IntVectorPtr	*fInfo;		// length is fInfo[nGridPoints]->vec[16]
	VectorPtr		fTime;		// workspace for vector of solution; fTime->vec[fNOfEquations]
	MatrixPtr		fSolution;	// workspace for vector of solution
	MatrixPtr		fSolPrime;	// workspace for partial derivatives
								// lenght is fSolution->mat[nGridPoints][fNOfEquations]
	int				fIdid;
	int				fML;
	int				fLRW;
	VectorPtr		*fRWork;	// workspace for ddassl
								// lenght is fRWork[nGridPoints]->vec[fLRW]
	int				fLIW;
	IntVectorPtr	*fIWork;	// workspace for ddassl
								// lenght is fIWork[nGridPoints]->vec[fLIW]
	int				fDasslNEq;	// = fNOfEquations, but dassl needs non constant
	Double			**fDmDy;	// dMdY[i][j] = dM_j/dY_i; dMdY[nSpeciesIn+1][nSpeciesIn]
};

#include "TTransFlameSolver.hpp"

#endif // __TTRANS_FLAME_SOLVER_H__
