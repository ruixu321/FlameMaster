#ifndef TFLAME_H__
#define TFLAME_H__

#include <vector>

#include "Config.h"
#include "Constants.h"
#include "Input.h"
#include "TProperties.h"
#include "TReaction.h"
#include "TSoot.h"
#include "TRadiation.h"

#include "TSpecies.h"

template<class Species>
class TFlame {
public:
	TFlame( const FirstInput& firstInp ) { InitTFlame( &firstInp ); };
	virtual ~TFlame( void );

/*	virtual TSpecies<TransportModel>*		GetSpecies( void ) = 0;
	virtual TPropertiesPtr	GetProperties( void ) = 0;
	virtual TReactionPtr	GetReaction( void ) = 0;*/
//	TSootPtr		GetSoot( void ) { return fSoot; };
	TInputDataPtr	GetInputData( void ) { return fInputData; };
	
	Double			dmdYAnalyt( int mNumber, int yNumber, Double *Y, Double mixDensity
								, Double *reactionRate, Double *theReactionRate );
	void			dmdTAnalyt( int mNumber, Double temp, Double &dmdT, Double *reactionRate );
	virtual Double	GetPressure( void ) { return fPressure->vec[fPressure->len]; };
	void			SetPressure( Double press ) { fPressure->vec[fPressure->len] = press; };
	Flag			ClipNegativeConcs( void ) { return fClipNegativeConcs; };
	int				CheckSolution( Double &temp, Double *Y, int YLen );

	IntVectorPtr	GetFuelIndexVec( void ) { return fFuelIndex; };
	int				GetFuelIndex( int i ) { return fFuelIndex->vec[i]; };
	int				GetFuelIndex( void ) { return fFuelIndex->vec[0]; };
	int				GetNFuels( void ) { return fFuelIndex->len; };
	int				GetFromSpeciesIndex( void ) { return fFromSpeciesIndex; };
	int				GetToSpeciesIndex( void ) { return fToSpeciesIndex; };
	Double			GetNu( char *globalReaction, const char *name );
	Double			GetNuProduct( char *globalReaction, const char *name );
	char			*GetDateCreated( void ) { return fDateCreated; };
	char			*GetAuthor( void ) { return fAuthor; };
	char			*GetOutputPath( void ) { return fOutputPath; };
	char			*GetOutFileBuff( void ) { return fOutFile; };
	FILE 			*GetOutfile( const char *name, FileType type );
	char			*GetOutfileName( const char *name, FileType type );

    string          fRadiationName;
    Double          fRadiativeFrac;

protected:
	VectorPtr 		fPressure;
	TInputDataPtr	fInputData;
//	TReactionPtr	fReaction;
//	TSpecies<TransportModel>*		fSpecies;
//	TPropertiesPtr	fProperties;
//	TSootPtr		fSoot;
	ContinType		fContinType;
	Double			fContBound;
	Double			fContInc;
	Flag			fUseNumericalJac;
	Flag			fUseNumericalDM;

	// Sensitivity Analysis
	Flag			fSensAnal;
	Flag                    fSensObjAll;         // all species are SensObj
	Flag                    fSensAnalSpec;       // Sens Anal on Species
	Flag                    fSensAnalReac;       // Sens Anal on Reactions
	// cai: 08/24/2012
	Flag                    fSensAnalClas;       // Sens Anal on Classes
	Double                  fSensAnalFac;        // Factor by which A is multiplied
	Flag                    fSensMax;            // sensitivity on max and location of max
	Flag                    fSensFinal;          // sensitivity on final values
	char			**fSensObj;          // Name of species/variables used as target
	int			    fNSensObj;
	Flag			fReactionFluxes;

	Flag			fClipNegativeConcs;

	void			PrintFlameletVector( int len, Double *vec, const char *name, FILE *fp, int inc = 1 );
	void			NextPressure( void ) { fPressure->len++; };
	Flag			IsLastPressure( void ) 
					{ 
						return ( fPressure->len+1 == fPressure->phys_len ) 
								? TRUE : FALSE; 
					};

    private:
	void	InitTFlame( const FirstInput* firstInp ); 	// allocate and fill classes

	char			*fOutputPath;
	char			*fOutFile;
	IntVectorPtr	fFuelIndex;
	int				fFromSpeciesIndex;
	int				fToSpeciesIndex;
	char			*fDateCreated;
	char			*fAuthor;
};

template<class Species>
class T0DFlame : public TFlame<Species> {
public:
	using TFlame<Species>::fInputData;
	using TFlame<Species>::GetInputData;
	using TFlame<Species>::GetFuelIndex;
	T0DFlame( const FirstInput& firstInp ) 
		: TFlame<Species>( firstInp ),
		fReaction( fInputData ), 
		fSpecies( fInputData, &fReaction )
	{ 
		InitT0DFlame();
	};	
	virtual ~T0DFlame( void );

	// access functions
	Species*		    GetSpecies( void ) { return &fSpecies; };
	T0DPropertiesPtr	GetProperties( void ) { return fProperties; };
	T0DReactionPtr		GetReaction( void ) { return &fReaction; };
	T0DSootPtr			GetSoot( void ) { return fSoot; };
	string              fRadiationName;

	virtual int		    GetOffsetFirstSpecies( void ) = 0;
	virtual int		    GetOffsetTemperature( void ) = 0;
	virtual int		    GetVariablesWithoutSpecies( void ) = 0;

	void			    CompLewisNumbers( const char *lewisFile );

	void			    UpdateThermoProps( Double *Y, Double temp, Double &pressure, Double &density, EqOfState what, Double *sootMoments );

	void			    ComputeProperties( Double temp, Double *Y, Double &pressure, Double &density, EqOfState what, Double *sootMoments );
	Double 			    ComputeZBilger( Double *Y, Double *YFuelSide, Double *YOxSide );
	string              GetRadiationName(void){return fRadiationName;}  


protected:

	T0DReaction		    fReaction;
	Species	            fSpecies;
	T0DPropertiesPtr	fProperties;
	T0DSootPtr			fSoot;

	Double				GetElementMassFraction( Double *Y, const char * const atomName, Double atomMolarMass );

private:
	void 			InitT0DFlame( void );

	// solution
//	VectorPtr		fMassFractions;
//	Double			fTemp;
};

template<class Species>
class T1DFlame : public TFlame<Species> {
typedef void ( *JacFuncPtr )( TFlameNodePtr fFlameNode );
public:
	using TFlame<Species>::fInputData;
	using TFlame<Species>::fUseNumericalJac;
	using TFlame<Species>::fUseNumericalDM;
	using TFlame<Species>::CheckSolution;
	using TFlame<Species>::fContInc;
	using TFlame<Species>::fContBound;
	using TFlame<Species>::fContinType;
	using TFlame<Species>::fClipNegativeConcs;
	using TFlame<Species>::fPressure;
	using TFlame<Species>::fNSensObj;
	using TFlame<Species>::fSensObj;
	using TFlame<Species>::GetPressure;
	using TFlame<Species>::GetInputData;
	using TFlame<Species>::SetPressure;
	using TFlame<Species>::GetFuelIndex;
	using TFlame<Species>::GetToSpeciesIndex;
	using TFlame<Species>::GetFromSpeciesIndex;
	using TFlame<Species>::GetNFuels;
	using TFlame<Species>::GetFuelIndexVec;
	using TFlame<Species>::GetOutFileBuff;
	using TFlame<Species>::GetOutputPath;
	using TFlame<Species>::GetNu;
	T1DFlame( const FirstInput& firstInp ) 
		: TFlame<Species>( firstInp ), 
		fReaction( fInputData ), 
		fSpecies( fInputData, &fReaction )
	{ 
		InitT1DFlame();
	};	
	virtual ~T1DFlame( void );

	//void				UpdateThermoProps( vector<double> physicalGrid ); // nowhere use this
	void				UpdateThermoProps( vector<double> physicalGrid, double Delta, double ChiRef, vector<double> HRR, string RadiationModel);// By Alexis
	void				UpdateThermoProps( void );// Original
	void				UpdateThermoPropsRadiation( void );// By Alexis
	void				UpdateThermoProps( TFlameNodePtr flameNode, NodeInfoPtr nodeInfo ); // used in T1DFlame<Species>::RHSAction
	void				ComputeProperties( TFlameNodePtr flameNode, Double temp
									, Double *Y, Double pressure );
	void				ComputePropertiesMod( TFlameNodePtr flameNode, Double temp
									, Double *Y, Double pressure );
	TBVPSolverPtr		GetSolver( void ) { return fSolver; };//
	T1DReaction*			GetReaction( void ) { return &fReaction; };
	Species*				GetSpecies( void ) { return &fSpecies; };
	T1DPropertiesPtr	GetProperties( void ) { return fProperties; };
	T1DSootPtr			GetSoot( void ) { return fSoot; };
	VectorPtr			GetStrainRateVector( void ) { return fStrainRate; };
	VectorPtr			GetTemperature( void ) { return fSolTemp; };
	MatrixPtr			GetMassFracs( void ) { return fSolMassFracs; };
//	VectorPtr			GetSavedTemperature( void ) { return fSavedTemp; };
//	MatrixPtr			GetSavedMassFracs( void ) { return fSavedMassFracs; };
	Double				GetStrainRate( void ) { return fStrainRate->vec[fStrainRate->len]; };
	void				SetStrainRate( Double strainRate ) { fStrainRate->vec[fStrainRate->len] = strainRate; };
	Flag				UseDiffCorr( void ) { return !fNoDiffCorr; };//
	void				UnSetDiffCorr( void ) { fNoDiffCorr = TRUE; };//
	VectorPtr			GetDiffCorr( void ) { return fDiffusivityCorrection; };//
	Double				ThermoDiffusion( int speciesIndex, CoordType coordinate, NodeInfoPtr nodeInfo );
	void				FillJacThermoDiffusion( int nVariable, Double constCoeff, CoordType coordinate, NodeInfoPtr nodeInfo );
	void				ComputeDiffusivityCorrection( Double **Y, NodeInfoPtr nodeInfo );//
	virtual FILE		*GetOutputFile( const char *head, const char *tail, const FileType type ) = 0;
	void CompLewisNumbers( const char *lewisFile, const std::string localLewisFileName = "" );

    // radiation

    string              GetRadiationName(void){return fRadiationName;}
	
	/*11082014 Chemical equilibrium*/
	/*void			ChemEquilSolver(Double temp);*/
	
	void 				SetFlameNode( int k );
	void				PostConvergence( void *object );
	void				SetMixtureSpecificationLeft( int num ) { fMixtureSpecificationLeft = num; };
	void				SetMixtureSpecificationRight( int num ) { fMixtureSpecificationRight = num; };
	int					GetMixtureSpecificationLeft( void ) { return fMixtureSpecificationLeft; };
	int					GetMixtureSpecificationRight( void ) { return fMixtureSpecificationRight; };
	TFlameNodePtr		GetFlameNode( void ) { return fFlameNode; };
	TFlameNodePtr		GetFlameNodeSaved( void ) { return fFlameNodeSaved; };
	void				SetBVPInput( TBVPInputPtr input );
	void				FillJacDiffusion( int nVariable, int nEquation, Double constCoeff, Double *diffCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double				StandardDiffusion( int nVariable, Double *diffCoeff, NodeInfoPtr nodeInfo );
	void				FillJacMixFracDiffusion( int nVariable, int nEquation, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double				SecondDerivMixFracDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double				GetZStoich( void );
	Double				GetZStoich_mf( void );
	Double				GetGeometry( void ) { return fGeometry; };
	void				XToEta( TNewtonPtr bt, VectorPtr etaVec );
	virtual void		EtaToX( TNewtonPtr /*bt*/, VectorPtr /*xPhysVec*/ ) { cerr << "#error: no definition of virtual function EtaToX" << NEWL;};
	void				OriginToZstoich( VectorPtr xVec, VectorPtr mixtureFraction, Double zStoich );
	Double				ComputeEmissionIndex( int speciesIndex, Double *x );
	void				SaveSolution( void );
	void				RestoreSolution( void );
	Flag				AdjustGrid( PostIterFuncPtr PostIter );
	int					SensitivityAnalysis( Double coeffSpecies, Double coeffTemp, CoordType coordinate );
	void				ReactionFluxes( CoordType coordinate );
	Double				dmdYAnalyt( int mNumber, int yNumber, Double *Y, TFlameNodePtr flameNode );
	void				FilldMdYOnePointNum( TFlameNodePtr flameNode );
	void				FilldMdTOnePointNum( TFlameNodePtr flameNode );
	void				FilldMomdMomOnePointNum( TFlameNodePtr flameNode );
	void				FilldMdYOnePointAnal( TFlameNodePtr flameNode );
	void				FilldMdTOnePointAnal( TFlameNodePtr flameNode );
	void				FilldMdYOnePoint( TFlameNodePtr flameNode );
	void				FilldMdTOnePoint( TFlameNodePtr flameNode );
	void				FilldMomdMomOnePoint( TFlameNodePtr flameNode );
	void				(T1DFlame::*FilldMdYOnePointPtr)( TFlameNodePtr flameNode );
	void				(T1DFlame::*FilldMdTOnePointPtr)( TFlameNodePtr flameNode );
	void				(T1DFlame::*FilldMomdMomOnePointPtr)( TFlameNodePtr flameNode );
	void				WriteRoggFiles( TNewtonPtr bt );
	Double 				ComputeZBilger( Double *Y, Double *YFuelSide, Double *YOxSide );


	// radiation

	string       fRadiationName;
	Double       fRadiativeFrac;

	virtual void	PrintRHSTemp( TNewtonPtr bt ) = 0;
	virtual void	UpdateSolutionOnePoint( Double * /*y*/, int /*gridPoint*/ ) { cerr << "#error: wrong instance of function UpdateSolution called" << NEWL; exit(2); };

	virtual int		GetOffsetVVelocity( void ) = 0;
	virtual int		GetOffsetUVelocity( void ) = 0;
	virtual int		GetOffsetTemperature( void ) = 0;
	virtual int		GetOffsetMixFrac( void ) = 0;
	virtual int		GetOffsetFirstSpecies( void ) = 0;
	virtual int		GetVariablesWithoutSpecies( void ) = 0;
	virtual ConstStringArray GetVariableNames( void ) = 0;

protected:
	TBVPSolverPtr	fSolver;
	TFlameNodePtr	fFlameNode;
	TFlameNodePtr	fFlameNodeSaved;

	VectorPtr		fSolTemp;
	MatrixPtr		fSolMassFracs;

	VectorPtr		fSavedTemp;
	MatrixPtr		fSavedMassFracs;
	VectorPtr		fSavedGrid;

	VectorPtr		fMmod;

	Flag			fNoDiffCorr;
	VectorPtr		fDiffusivityCorrection; // [m^2/s]	//  this vector contains sum(D_k dY_k/dy)
	//	it has [nGridPoints+2] elements
	VectorPtr		fStrainRate;
	Double			fGeometry;   // 0.0 for planar, 1.0 for axi-symmetric geometry
	Flag			fThermoDiffusion;
	ContinSide		fContinSide;
	Flag			fPrintRHSSpecies;
	Flag				fPrintRHSTemp;
	T1DReaction			fReaction;
	Species				fSpecies;
	T1DPropertiesPtr	fProperties;
	T1DSootPtr			fSoot;

	void			SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void			CheckBC( void );
	void			ReadStartProfiles( TInputDataPtr inp );
	void			CheckInitialGuess( void );
	void			UpdateDimensions( int len );
	void			UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void			UpdateSolution( Double *y, int gridPoint );
	void			SolutionToSolver( void );
	void 			CopyFlameNode( TFlameNodePtr flameNodeSource, TFlameNodePtr flameNodeDest );
	void			CompareFlameNode( TFlameNodePtr flameNodeSource, TFlameNodePtr flameNodeDest );
	Flag			RHSAction( NodeInfoPtr nodeInfo, RHSMode rhsMode );
	Flag			RHSAction(void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
	Double			GetElementMassFraction( Double *Y, const char * const atomName, Double atomMolarMass );

private:
	void			InitT1DFlame( void );
	TFlameNodePtr 	NewTFlameNodeSaved( void );
	void 			DisposeTFlameNodeSaved( TFlameNodePtr flameNodeSaved );
	int				fMixtureSpecificationLeft; // uses enum MixtureSpecification
	int				fMixtureSpecificationRight;
	virtual void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp ) = 0;
	void			SaveGrid( void );
	void			RestoreGrid( void );
	int				GetNOfX( Double theX, int nGridPoints, Double *x );
	int				CheckFlameLocation( Double xMid, Double deltaXVar );
	Double			CheckGradientLeft( void );
	Double			CheckGradientRight( void );
	void			PrintSensitivityCoefficients( TensorPtr sensitivityCoefficients, IntVectorPtr objects );
	void			PrintReactionFluxesReac( Double **fluxes, const char *name, const char *header );
	void 			PrintReactionFluxesSpec( Double **fluxes, const char *name, const char *header );
	void			PrintReduceInfo( Double **fluxes, Double epsilon, const char *header );
	void			PrintOneSensMax( Double **sensCoeff, int gridPoints, int nReactions
			, const char *fileName );
	Double			GetMaxSens( Double **sensCoeff, int gridPoints, int reaction );
	void			PrintSensMax( TensorPtr sensitivityCoefficients, IntVectorPtr objectsPtr );
};

#include "TFlame.hpp"

#endif // TFLAME_H__
