#ifndef TCOUNT_DIFF_FLAME_MIX_H__
#define TCOUNT_DIFF_FLAME_MIX_H__

#include <time.h>

#include "Constants.h"
#include "TFlame.h"
#include "RadiationModel.h"

template<class Flame>
void CountDiffMixJacRest( void *object, NodeInfoPtr nodeInfo );
template<class Flame>
void CountDiffMixRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
template<class Flame>
void CountDiffMixOutput( void *object, FILE *fp, const char* tail );
template<class Flame>
int CountDiffMixPostIter( void *object );
template<class Flame>
void SetCountDiffMixNodeInfo( int k, void *object );
template<class Flame>
void CountDiffMixPostConv( void *object );
template<class Flame>
ConstStringArray GetCountDiffMixVarNames( void *object );
template<class Flame>
void CountDiffMixUpdateBoundary( void  *object );

template<class Species>
class TCountDiffFlameMix : public T1DFlame<Species> {
    template<class Flame>
        friend void CountDiffMixJacRest( void *object, NodeInfoPtr nodeInfo );
    template<class Flame>
        friend void CountDiffMixRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
    template<class Flame>
        friend void CountDiffMixOutput( void *object, FILE *fp, const char* tail );
    template<class Flame>
        friend int CountDiffMixPostIter( void *object );
    template<class Flame>
        friend void CountDiffMixPostConv( void *object );

    public:
    using T1DFlame<Species>::fInputData;
    using T1DFlame<Species>::fSpecies;
    using T1DFlame<Species>::fReaction;
    using T1DFlame<Species>::GetSolver;
    using T1DFlame<Species>::fSoot;
    using T1DFlame<Species>::fUseNumericalJac;
    using T1DFlame<Species>::CheckBC;
    using T1DFlame<Species>::fSolver;
    using T1DFlame<Species>::CheckInitialGuess;
    using T1DFlame<Species>::fFlameNode;
    using T1DFlame<Species>::UseDiffCorr;
    using T1DFlame<Species>::GetMixtureSpecificationLeft;
    using T1DFlame<Species>::GetTemperature;
    using T1DFlame<Species>::GetMassFracs;
    using T1DFlame<Species>::GetSoot;
    using T1DFlame<Species>::GetPressure;
    using T1DFlame<Species>::GetInputData;
    using T1DFlame<Species>::GetSpecies;
    using T1DFlame<Species>::fStrainRate;
    using T1DFlame<Species>::fProperties;
    using T1DFlame<Species>::SetPressure;
    using T1DFlame<Species>::GetStrainRate;
    using T1DFlame<Species>::GetFuelIndex;
    using T1DFlame<Species>::SetMixtureSpecificationRight;
    using T1DFlame<Species>::SetMixtureSpecificationLeft;
    using T1DFlame<Species>::GetProperties;
    using T1DFlame<Species>::UpdateThermoProps;
    using T1DFlame<Species>::fSolTemp;
    using T1DFlame<Species>::GetOutfile;
    using T1DFlame<Species>::GetOutFileBuff;
    using T1DFlame<Species>::GetOutputPath;
    using T1DFlame<Species>::fThermoDiffusion;
    using T1DFlame<Species>::ThermoDiffusion;
    using T1DFlame<Species>::SetFlameNode;
    using T1DFlame<Species>::ReadStartProfiles;
    using T1DFlame<Species>::ComputeProperties;
    using T1DFlame<Species>::CompLewisNumbers;
    using T1DFlame<Species>::SetInitialBC;
    using T1DFlame<Species>::GetZStoich;
    using T1DFlame<Species>::fSavedTemp;
    using T1DFlame<Species>::UnSetDiffCorr;
    using T1DFlame<Species>::PrintFlameletVector;
    TCountDiffFlameMix( const FirstInput& firstInp ) : 
        fFirstSpecies(0), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
        fLnChi(fTemperature+1),
        fSootMoments(fLnChi+1),
        fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
        T1DFlame<Species>( firstInp ) { 

            InitTCountDiffFlameMix();};	
    ~TCountDiffFlameMix( void );

    int			GetOffsetVVelocity( void );
    int			GetOffsetUVelocity( void );
    int			GetOffsetMixFrac( void );
    int			GetOffsetFirstSpecies( void );
    int			GetOffsetTemperature( void );
    int			GetVariablesWithoutSpecies( void );
    int			GetOffsetLnChi( void );
    ConstStringArray	GetVariableNames( void );
    VectorPtr	GetDissRateVector( void ) { return fDissRate; };
    //	Double		GetDissRate( void ) { return fDissRate->vec[fDissRate->len]; };
    Double		GetDissRate( void ) { return (fArcLengthContin) ? exp( fSolLnChi->vec[1] ) : fDissRate->vec[fDissRate->len]; };
    Double		GetDissRate( Double z );
    Double		DissRateFact( Double rho );
    //	Double		GetDissRate( Double z, Double rho );
    Double		GetDissRate( Double z, Double rho, Double chiStoich );
    int			GetZRefLoc( VectorPtr xVec );
    Double		GetZRef( void );

    // continuation stuff
    //	void		IncTempContStart( Double inc ) { fTempContStart += inc; };
    void		SetTempContStart( Double val ) { fTempContStart = val; };
    void		SetChiContStart( Double val ) { fLnChiContStart = val; };
    void		SetDeltaArcLength( Double val ) { fDeltaArcLength = val; };
    void		SetTMaxLoc( int loc ) { fTmaxLoc = loc; };
    void		SetdTds( Double dTds ) { fdTds = dTds; };
    void		SetdlnChids( Double dlnChids ) { fdlnChids = dlnChids; };
    void		IncNFlameletCount( void ) { ++fNFlameletsCount; };
    void		ReInitArcLengthCont( void ) { fDeltaArcLength = 0.0; fdTds = 0.0; fdlnChids = 0.0; fNFlameletsCount = 0; };
    Flag		GetArcLengthCont( void ) { return fArcLengthContin; };
    Flag		GetArcUp( void ) { return fArcUp; };
    int			GetTMaxLoc( void ) { return fTmaxLoc; };
    int			GetMaxFlamelets( void ) { return fMaxFlamelets; };
    int			GetNFlameletsCount( void ) { return fNFlameletsCount; };
    Double		GetdTds( void ) { return fdTds; };
    Double		GetdlnChids( void ) { return fdlnChids; };
    Double		GetChiContStart( void ) { return fLnChiContStart; };
    Double		GetTempContStart( void ) { return fTempContStart; };
    Double		GetDeltaArcLength( void ) { return fDeltaArcLength; };
    Double		GetDeltaChiref( void ) { return fDeltaChiRef; };
    Double		GetDeltaTref( void ) { return fDeltaTRef; };
    int			GetLnChi( void ) { return fLnChi; };

    int   Niter = 0; // number that check if it is the first iter

    protected:
    void	PrintRHSTemp( TNewtonPtr /*bt*/ ) { cerr << "nothing happens" << NEWL; };

    private:


    const int	fFirstSpecies;
    const int	fTemperature;
    const int	fLnChi;
    const int	fSootMoments;
    const int	fVariablesWithoutSpecies;
    char		**fVariableNames;
    Flag		fPrintMolarFractions;
    VectorPtr	fSolLnChi;	// 
    VectorPtr	fSavedLnChi;		// 
    Double		fRhoInf;
    Double		fRhoRef;
    VectorPtr	fZCount;
    VectorPtr	fChiCount;
    VectorPtr	fDissRate;

    // continuation stuff
    int			fTmaxLoc;			//	location of maximum temperature
    int			fMaxFlamelets;
    int			fNFlameletsCount;
    Flag		fArcLengthContin;
    Flag		fArcUp;
    Double		fTempContStart;
    Double		fLnChiContStart;
    Double		fSavedTempContStart;
    Double		fSavedLnChiContStart;
    Double		fDeltaArcLength;
    Double		fdTds;				//	dT_max/ds, where s is the arclength
    Double		fDeltaChiRef;
    Double		fDeltaTRef;
    Double		fdlnChids;			//	da_inv/ds, where s is the arclength and a_inv the inverse of the strainrate

    // Liquid fuel stuff
    Double		fCl; // cp of liquid
    Double		fHv; // enthalpy of evaporation at 1atm and boiling conditions
    Double		fT_B;// boiling temperature 

    // radiation

    Double      fRadiativeFrac; 
    string      fRadiationName; 
    string      fZlZr; 

    void	InitTCountDiffFlameMix( void );
    void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
    FILE	*GetOutputFile( const char *head, const char *tail, const FileType type );
    void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
    void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign );
    void	UpdateDimensions( int len );
    void	SaveSolution( void );
    void	RestoreSolution( void );
    void	SolutionToSolver( void );
    void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
    void	UpdateSolutionOnePoint( Double *y, int gridPoint );
    Double	ExactChi( Double Z );
    Double	Interpol( Double x, Double xOld, Double yOld, Double xNew, Double yNew );
    void    PrintProdRate( int speciesIndex, FILE *fp );

};

#include "TCountDiffFlameMix.hpp"

#endif // TCOUNT_DIFF_FLAME_MIX_H__
