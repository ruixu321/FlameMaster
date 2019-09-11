#ifndef TTRANS_FLAMELET_H__
#define TTRANS_FLAMELET_H__

#include "Constants.h"
#include "TTransFlameSolver.h"

template<class Species>
class TTransFlamelet : public TTransFlameSolver<Species> {
public:
  using TTransFlameSolver<Species>::fInputData;
  using TTransFlameSolver<Species>::fVariables;
  using TTransFlameSolver<Species>::fSpecies;
  using TTransFlameSolver<Species>::fSoot;
  using TTransFlameSolver<Species>::GetFuelIndex;
  using TTransFlameSolver<Species>::GetSoot;
  using TTransFlameSolver<Species>::GetOutfile;
  using TTransFlameSolver<Species>::GetNFuels;
  using TTransFlameSolver<Species>::WriteFlameletFile;
  using TTransFlameSolver<Species>::WriteFlameletFileInitialize;
  using TTransFlameSolver<Species>::TurbMeanY;
  using TTransFlameSolver<Species>::TurbMeanX;
  using TTransFlameSolver<Species>::TurbMeanSoot;
  using TTransFlameSolver<Species>::PrintSolution;
  using TTransFlameSolver<Species>::MakeGrid;
  using TTransFlameSolver<Species>::CompLewisNumbers;
  using TTransFlameSolver<Species>::TurbMeanTotEnergy;
  using TTransFlameSolver<Species>::TurbMeanTemp;
  using TTransFlameSolver<Species>::ComputeEmissionIndex;
  using TTransFlameSolver<Species>::ComputeZBilger;
  using TTransFlameSolver<Species>::ComputeMeanFromPDF;
  using TTransFlameSolver<Species>::CheckPDF;
  using TTransFlameSolver<Species>::ComputeEmissionIndexSoot;
  
  using TTransFlameSolver<Species>::PrintFlameletVector;
  using TTransFlameSolver<Species>::GetOutputFile;
  using TTransFlameSolver<Species>::Interpol;
  using TTransFlameSolver<Species>::Initialize;
	TTransFlamelet( const FirstInput& firstInp ) : fTemperature(0), 
			fFirstSpecies(1),
			fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
			fNOfEquations(fInputData->GetCounter()->species+fVariablesWithoutSpecies - fInputData->GetCounter()->steadyStates),
			TTransFlameSolver<Species>( firstInp ) { InitTTransFlamelet();};	
	virtual ~TTransFlamelet( void );
	void				Solve( void );
	ConstStringArray	GetVariableNames( void ) { return fVariableNames; };

private:
  	void	InitTTransFlamelet( void );
	void	SetInitialBC( void );
	void	SetInitialValues( void );
	Double	GetTempEnd( Double time, Double timeStart, Double timeCut
									, Double tempStart, Double tempEnd );
	void	ReadStartProfiles( TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	Double	GetScalarDiss( Double theTime );
	void	ReadPDF( void );
	void	WriteScalars( Double time, int i );

	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fVariablesWithoutSpecies;
	const int	fNOfEquations;
	char		**fVariableNames;

	int fDeltaStepsOut;
	int			fPDFGridPoints;
	VectorPtr	fTimePDFIn;
	MatrixPtr	fPDF;
	FILE		*ffpEI;


	Double		fRPM;
	int			fVarsIn;
	VectorPtr	fTimeIn;
	VectorPtr	fTFuelIn;
	VectorPtr	fTOxIn;
	VectorPtr	fPressureIn;
	VectorPtr	fChiIn;
	VectorPtr	fZRIn;
	VectorPtr	fZMeanIn;
	VectorPtr	fZVarIn;
	VectorPtr	fXOverD;
	Flag		fUseInput;
	Double		fTStart;
	Double		fTEnd;
	Double		fTCurr;
	Double		fDeltaT;
	Double		fDeltaTStart;
	Double		fTimeCut;
	int			fNOutputs;
	Double		fScalarDissRate;
	
	Double		fPressStart;
	Double		fTempOxStart;
	Double		fTempFuelStart;
	Double		fPressEnd;
	Double		fChiEnd;
	Double		fTempOxEnd;
	Double		fTempFuelEnd;

	VectorPtr	fGridSol;
	MatrixPtr	fStartSol;
	Flag		fEquidistant;
	int			fNGridPoints;
	int			fMaxGridPoints;
	int			*fBCFlagLeft;
	int			*fBCFlagRight;
};

#include "TTransFlamelet.hpp"

#endif // TTRANS_FLAMELET_H__
