#ifndef __INPUT_H__
#define __INPUT_H__

#include <cstring>   // strcpy
#include <iostream>  // cerr
#include <memory>    // unique_ptr, nullptr

#include "Config.h"
#include "Constants.h"
#include "FlameUtilities.h"
#include "ListTool.h"            // for definition of FatalError
#include "Newton.h"              // enum{kPath, kFileName};
#include "ReadTransportModel.h"  // makeReadSpeciesData factory
#include "ScanManStructs.h"

struct BoundaryInput;
typedef BoundaryInput *BoundaryInputPtr;

enum class TransportModelType { MONO_ATOMIC, MULTI_GEOMETRY, FITTED };
enum class PreFileType { VERSION_3311, VERSION_4000, VERSION_4000_FIT };
enum class FileType	{ kNone, kData, kText, kAppend, kBin };//Krithika

enum FlameType { kNotSpecified, kStartUpDiffusion, kDiffusionPhys, kDiffPhysEigen, kStartUpDiffusion2, 
                kCountPremPhys, kCountPremSim, kCountDiffMix, kCountDiffCont,
                kHomoIsoChor, kHomoIsoBar, kHomoPSR, kTransFlamelet, kTrans1DIsoChor, kUnstrPremPhys, kTransCoalParticle, kEquilibrium };

struct FirstInput {
  FirstInput(int argc, char *argv[])
      : fVVelocityOffset(0),
        fUVelocityOffset(1),
        fTemperatureOffset(2),
        fPressureOffset(3) {
    InitFirstInput(argc, argv);
  };
  //FirstInput(char *fInput)
  //    : fVVelocityOffset(0),
  //      fUVelocityOffset(1),
  //      fTemperatureOffset(2),
  //      fPressureOffset(3) {
  //  InitFirstInput(0, &fInput);
  //};
  ~FirstInput(void);

  void PrintAdditionalData(void);

  const int fVVelocityOffset;
  const int fUVelocityOffset;
  const int fTemperatureOffset;
  const int fPressureOffset;

  TransportModelType TMT = TransportModelType::MONO_ATOMIC;
  PreFileType PFT = PreFileType::VERSION_4000;

  // Exact Backward Reaction Constants
  Flag fExactBackward;

  // 0D-stuff
  Flag fAdditionalOutput;
  Flag fEquidistant;
  int fNOutputs;
  Double fArtificialSource;
  Double fArtSourceTime;
  // For PSR model - krithika
  Flag fIsothermal;
  Double fHeatTransferCoefficient;
  Double fAmbientTemperature;
  // Cai
  Double fTres;
  Double fDpdt;

  // rest
  Flag fScannerProgress;

  char *fDateCreated;
  char *fAuthor;
  int fFlameType;  // uses enum FlameType
  String fContinType;
  String fContinSide;
  Double fContBound;
  Flag fConstantLewisNumber;
  Flag fWithRadiation;
  Flag fArcLengthContin;
  Flag fStrainRateContin;
  Flag fThermoDiffusion;
  Flag fWithSoot;
  int fNSootMoments;
  Flag fCompUnPhysChain;
  int fNUnPhysChain;
  Flag fIsAxiSymmetric;
  Flag fClipNegativeConcs;
  Flag fNoDiffCorr;
  Double fMinStrainRate;
  Double fPressure[1000];
  int fNPressures;
  Double fPressureComm;   // pressure specified on command line
  Double fParameterComm;  // parameter specified on command line
  Double fRadiativeFrac;  // parameter specified on command line
  String fRadiationName;
  String fZlZr;
  int fPremConfiguration;

  // Sensitivity Analysis
  Flag fSensAnal;
  Flag fSensObjAll;    // all species are SensObj
  Flag fSensAnalSpec;  // Sens Anal on Species
  Flag fSensAnalReac;  // Sens Anal on Reactions
  // cai: 24/08/2012
  Flag fSensAnalClas;    // Sens Anal on Classes
  Double fSensAnalFac;   // Factor by which A is multiplied
  Flag fSensMax;         // sensitivity on max
  Flag fSensFinal;       // sensitivity on final values
  char *fSensObj[1000];  // Name of species/variables used as target
  int fNSensObj;
  Flag fReactionFluxes;

  Flag fPrintRHSSpecies;
  Flag fPrintRHSTemp;
  Double fStrainRate[1000];
  int nStrainRates;
  Double fDissRate[1000];
  int nDissRates;
  int fNPhi;
  Double fPhi[1000];
  Flag fKeepMassFracs;
  Flag fLiquidPoolBC;
  Double fDeltaSCont;
  //	Double			fChi;
  char *fGlobalReaction;
  char *fFuelName[1000];
  int fNFuels;
  char *fFromSpeciesName;
  char *fToSpeciesName;
  Double fContInc;
  char *fOxName;
  Flag fPrintMolarFractions;
  Flag fSteadyStatesNumerical;
  int fVariablesWithoutSpecies;

  Flag fNucleation;
  Flag fCondensation;
  Flag fCoagulation;
  Flag fSurfaceGrowth;
  Flag fSurfaceOxidation;
  Flag fThermoPhoresis;
  Flag fOHPAHOxidation;
  Flag fO2PAHOxidation;
  Double fCoagFact;
  Flag fSootRadiation;
  Flag fSootUpdateProdRates;
  Flag fSizeDepDiff;
  Flag fSurfDepCoag;

  BoundaryInputPtr fInitialCond;
  BoundaryInputPtr leftBoundary;
  BoundaryInputPtr rightBoundary;

  Flag fWriteBT;
  Flag fWriteResiduum;
  Flag fWatchGridding;
  Flag fWriteEverySolution;
  char *fOutputPath;
  char *fOutputPathComm;

  Flag fUseModifiedNewton;
  Flag fUseNumericalJac;
  Flag fUseSecOrdJac;
  Flag fUseNumericalDM;
  Flag fWriteFullRes;
  int fInitialEquations;
  int fMaxGridPoints;
  int fInitialGridPoints;
  Flag fDampFlag;
  Flag fTimeDepFlag;
  Flag fContinFlag;
  int fDeltaNewGrid;
  Flag fOneSolOneGrid;
  Flag fAdjustComputationalDomain;
  Double fTolRes;
  Double fTolDy;
  int fMaxIter;
  Double fLeft;
  Double fRight;
  int fGamma;
  Double fKappa;
  Double fTauGrid;
  Double fRadLossPercent; //Rui
  Double fR;
  Double fQ;
  Double fStart;
  Double fEnd;
  Double fAlpha;
  Double fLambdaMin;
  Double fDeltaTStart;
  Double fDeltaTMax;
  int fContSteps;

  Double fCl;
  Double fHv;
  Double fT_B;

  char *fLewisNumberFile;
  char *fExpTempFile;
  char *fReactionFile;
  // THIS FILE NAME WAS GIVEN A LENGTH OF 31 BEFORE.
  // THIS LED TO THE FOLLOWING POINTER BEING OVERWRITTEN.
  // SHOULD CHECK THAT IT DOESN'T OVERWRITE ANYTHING OR DYNAMICALLY ALLOCATE
  char fAdditionalInfoFile[LENOFBUFFER];
  //	char			*fAdditionalInfoFile;
  char *fMechanismFileComm;     // name of mechanism file from the command line
  char *fStartProfileFileComm;  // name of startprofiles file from the command
                                // line
  char *fStartProfileFile;
  char *fOutFileName;
  char *fAddFileNo1;
  FILE *fpR;  // filepointer to fReactionFile
  FILE *fpA;  // filepointer to fAdditionalInfoFile
  FILE *fpS;  // filepointer to fStartProfileFile
  FILE *fpO;  // filepointer to fOutputFile

 private:
  void InitFirstInput(int argc, char *argv[]);
  void ParseCommandLine(int argc, char *argv[]);
  int yylex(void);

  char stringBuffer[LENOFBUFFER];
};

struct BoundaryInput {
  BoundaryInput(const FirstInput *firstInp, int len)
      : fVVelocityOffset(firstInp->fVVelocityOffset),
        fUVelocityOffset(firstInp->fUVelocityOffset),
        fTemperatureOffset(firstInp->fTemperatureOffset),
        fPressureOffset(firstInp->fPressureOffset),
        fNVars(4) {
    InitBoundaryInput(len);
  };
  ~BoundaryInput(void);

  BoundaryInput &operator=(const BoundaryInput &boundary);
  void PrintBoundary(FILE *fp);

  const int fVVelocityOffset;
  const int fUVelocityOffset;
  const int fTemperatureOffset;
  const int fPressureOffset;
  const int fNVars;

  // the following variables contain bc's for the species
  int fSpecifiedSpeciesBCs;
  int fMixtureSpecification;  //	uses enum MixtureSpecification
  int fBcFlagSpecies;         //	uses enum BoundaryCondition
  char **speciesName;
  Double *fValueSpecies;  //	contains value of bc for species

  // the following variables contain bc's for f ( or V ), f' ( or U ), T and p
  int *fBcFlag;    //	specifies kind of bc; uses enum BoundaryCondition
  Double *fValue;  //	contains value of bc

 private:
  void InitBoundaryInput(int len);
  int fLen;
};

class TInputData {
 public:
  TInputData(const FirstInput *firstInp)
      : fVVelocityOffset(firstInp->fVVelocityOffset),
        fUVelocityOffset(firstInp->fUVelocityOffset),
        fPressureOffset(firstInp->fPressureOffset),
        fTemperatureOffset(firstInp->fTemperatureOffset) {
    InitInputData(firstInp);
  };
  ~TInputData(void);

  HeaderPtr GetHeader(void) { return fHeader; };
  CounterPtr GetCounter(void) const { return fCounter; };
  ReactionPtr GetReactions(void) { return fReaction; };
  ReactionPtr GetPAHReactions(void) { return fPAHReaction; };
  ReactionPtr GetSootReactions(void) { return fSootReaction; };
  ThirdBodyPtr GetThirdBody(void) { return fThirdBody; };
  ReadSpeciesData *GetReadSpeciesData(void) const {
    return fReadTransportModel.get();
  };
  DimensionPtr GetDimension(void) { return fDimension; };
  int FindSpecies(ConstString species);
  int FindAtomIndex(ConstString atom);
  int FindSpeciesCoefficient(int speciesIndex, int reactionIndex);
  void PrintAdditionalData(void);

  const int fVVelocityOffset;
  const int fUVelocityOffset;
  const int fTemperatureOffset;
  const int fPressureOffset;

  PreFileType PFT;
  TransportModelType TMT;

  // Exact Backward Reaction Constants
  Flag fExactBackward;

  // 0D-stuff
  Flag fAdditionalOutput;
  Flag fEquidistant;
  int fNOutputs;
  Double fArtificialSource;
  Double fArtSourceTime;

  // For PSR
  Flag fIsothermal;
  Double fHeatTransferCoefficient;
  Double fAmbientTemperature;
  // Cai
  Double fTres;
  Double fDpdt;

  char *fDateCreated;
  char *fAuthor;
  int fFlameType;  // uses enum FlameType
  ContinType fContinType;
  ContinSide fContinSide;
  Double fContBound;
  Flag fConstantLewisNumber;

  // radiation
  Flag fWithRadiation;
  Double fRadiativeFrac;  
  string fRadiationName;
  string fZlZr;

  Flag fArcLengthContin;
  Flag fStrainRateContin;
  Flag fThermoDiffusion;
  Flag fWithSoot;
  int fNSootMoments;
  Flag fCompUnPhysChain;  // also flag for ConstMassFlux in TUnstrechedPremixed
  int fNUnPhysChain;
  Flag fIsAxiSymmetric;
  Flag fClipNegativeConcs;
  Flag fNoDiffCorr;
  Double fMinStrainRate;
  VectorPtr fPressure;
  Double fParameterComm;  // parameter specified on command line
  int fPremConfiguration;

  // Sensitivity Analysis
  Flag fSensAnal;
  Flag fSensObjAll;    // all species are SensObj
  Flag fSensAnalSpec;  // Sens Anal on Species
  Flag fSensAnalReac;  // Sens Anal on Reactions
  // cai: 24/08/2012
  Flag fSensAnalClas;   // Sens Anal on Classes
  Double fSensAnalFac;  // Factor by which A is multiplied
  Flag fSensMax;        // sensitivity on max and location of max
  Flag fSensFinal;      // sensitivity on final values
  char **fSensObj;      // Name of species/variables used as target
  int fNSensObj;
  Flag fReactionFluxes;

  Flag fPrintRHSSpecies;
  Flag fPrintRHSTemp;
  VectorPtr fStrainRate;
  VectorPtr fDissRate;
  VectorPtr fPhi;
  Flag fKeepMassFracs;
  Flag fLiquidPoolBC;
  Double fDeltaSCont;

  char *fGlobalReaction;
  IntVectorPtr fFuelIndex;
  int fFromSpeciesIndex;
  int fToSpeciesIndex;
  Double fContInc;
  int fOxIndex;
  int fH2OIndex;
  int fCO2Index;
  Flag fPrintMolarFractions;
  Flag fSteadyStatesNumerical;
  int fVariablesWithoutSpecies;
  char *fPAHSymbol;
  char *fSootSymbol;

  Flag fNucleation;
  Flag fCondensation;
  Flag fCoagulation;
  Flag fSurfaceGrowth;
  Flag fSurfaceOxidation;
  Flag fThermoPhoresis;
  Flag fOHPAHOxidation;
  Flag fO2PAHOxidation;
  Double fCoagFact;
  Flag fSootRadiation;
  Flag fSootUpdateProdRates;
  Flag fSizeDepDiff;
  Flag fSurfDepCoag;

  BoundaryInputPtr fInitialCond;
  BoundaryInputPtr leftBoundary;
  BoundaryInputPtr rightBoundary;

  Flag fWriteBT;
  Flag fWriteResiduum;
  Flag fWatchGridding;
  Flag fWriteEverySolution;
  char *fOutputPath;

  Flag fWriteFullRes;
  Flag fUseModifiedNewton;
  Flag fUseNumericalJac;
  Flag fUseSecOrdJac;
  Flag fUseNumericalDM;
  int fNVariables;
  int fInitialEquations;
  int fMaxGridPoints;
  int fInitialGridPoints;
  Flag fDampFlag;
  Flag fTimeDepFlag;
  Flag fContinFlag;
  int fDeltaNewGrid;
  Flag fOneSolOneGrid;
  Flag fAdjustComputationalDomain;
  Double fTolRes;
  Double fTolDy;
  int fMaxIter;
  Double fLeft;
  Double fRight;
  int fGamma;
  Double fKappa;
  Double fTauGrid;
  Double fRadLossPercent;  // Rui
  Double fR;
  Double fQ;
  Double fStart;
  Double fEnd;
  Double fAlpha;
  Double fLambdaMin;
  Double fDeltaTStart;
  Double fDeltaTMax;
  int fContSteps;

  Double fCl;
  Double fHv;
  Double fT_B;

  char *fLewisNumberFile;
  char *fExpTempFile;
  char *fStartProfileFile;
  char *fAddFileNo1;
  FILE *fpS;  // filepointer to fStartProfileFile
  char *fReactionFile;
  char fAdditionalInfoFile[LENOFBUFFER];
  char *fOutFileName;

 private:
  void InitInputData(
      const FirstInput *firstInp);  // calls ReadInReactions( void )
  void ReadInAdditionalData(void);
  void ReadInReactions(void);

  HeaderPtr NewHeader(void);
  HeaderPtr ReadHeader(void);
  CounterPtr ReadCounter(void);
  AtomsPtr NewAtomArray(int len);
  AtomsPtr ReadAtoms(void);
  //###SpeciesPtr 		NewSpeciesArray( int len );
  //###SpeciesPtr 		ReadSpecies( void );
  ReactionPtr NewReactionArray(int len);
  ReactionPtr ReadReaction(int len);
  ThirdBodyPtr NewThirdBodyArray(int len);
  ThirdBodyPtr ReadThirdBody(void);
  DimensionPtr NewDimensionArray(int len);
  DimensionPtr ReadDimension(void);
  void FreeAtomsArray(AtomsPtr atoms, int len);
  //###void			FreeSpeciesArray( SpeciesPtr species, int len );
  void FreeReactionArray(ReactionPtr reaction, int len);
  void FreeThirdBodyArray(ThirdBodyPtr thirdBody, int len);
  void FreeDimensionArray(DimensionPtr dimension, int len);
  int GetSpeciesNumber(void);

  CounterPtr fCounter;  // initialized in InitInputData( void )
  AtomsPtr fAtoms;
  std::unique_ptr<ReadSpeciesData> fReadTransportModel = nullptr;
  ReactionPtr fReaction;
  ReactionPtr fPAHReaction;
  ReactionPtr fSootReaction;
  ThirdBodyPtr fThirdBody;
  DimensionPtr fDimension;
  HeaderPtr fHeader;
  String fBuffer;  // initialized in ReadHeader( void )

  FILE *fpR;  // filepointer to fReactionFile
  FILE *fpA;  // filepointer to fAdditionalInfoFile
  FILE *fpO;  // filepointer to fReactionFile
};
typedef TInputData *TInputDataPtr;

#endif  // __INPUT_H__
