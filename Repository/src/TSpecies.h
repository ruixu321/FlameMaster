#ifndef TSPECIES__
#define TSPECIES__

#include <limits>
typedef std::numeric_limits< double > dbl;

#include "Config.h"
#include "Constants.h"
#include "Input.h"
#include "TProperties.h"
#include "TReaction.h"

template <class TransportModel>
class TSpecies {
  // TFlame contains a function UpdateThermoProps for the complete calculation
  // of all species and mixture properties and the calculation of the production
  // rates

  // compute all properties for steady states except for steady state species
  // have no storage for fNOfUsedReactions,fUsedReactions, fNu, fProductionRate

 public:
  TSpecies(TInputDataPtr input, TReactionPtr reaction) : TM(*input), fNOfSpecies(input->GetCounter()->species) {
    InitSpecies(input, reaction);
  };
  ~TSpecies(void);
  TransportModel *GetTransportModel(void) { return &TM; };
  int GetNOfSpecies(void) { return fNOfSpecies; };
  int GetNSteadyStates(void) { return fNSteadyStates; };
  int GetNSpeciesInSystem(void) { return fNOfSpecies - fNSteadyStates; };
  VectorPtr GetMolarMass(void) { return fMolarMass; };
  IntMatrixPtr GetComposition(void) const { return fComposition; };
  VectorPtr GetLewisNumber(void) { return fLewisNumber; };
  const char *GetLewisNumberFile(void) { return fLewisNumberFile; };
  VectorPtr *GetNu(void) { return fNu; };
  IntVectorPtr GetNOfUsedReactions(void) { return fNOfUsedReactions; };
  IntVectorPtr *GetUsedReactions(void) { return fUsedReactions; };
  // DoublePtr		***GetReactionRate( void ) { return fReactionRate; };
  char **GetNames(void) { return fNames; };
  Flag IsConstantLewisNumber(void) { return fConstLewisNumber; };
  Flag IsSteadyState(int i) { return fIsSteadyState[i]; };
  int FindSpecies(const char *species);
  void WriteRoggsSymbolsData(void);
  Double ComputeDCpiDT(int index, Double temp);
  //	Double			CompOnedRhoRhoDkdY( int speciesIndexY, int
  // speciesIndexD, Double *Y, TFlameNodePtr flameNode );
  //	Double			CompOnedRhoDkdY( int speciesIndexY, int
  // speciesIndexD,
  // Double *Y, Double rho, Double M, Double *D, Double **invDij );
  void ComputeProductionRates(Double *productionRate, Double *reactionRate);
  Double GetPosNegProductionRate(int which, Double *reactionRate, Flag pos);
  void ComputeNoProductionRates(Double *productionRate);
  void ComputeProductionRate(int i, Double &productionRate,
                             Double *reactionRate);
  void ComputeDeltaI(Double *deltaI, Double *Y, Double *viscosity);
  void CompProdCons(Double *source, Double *sink, Double *reacRate);
  void CompHeatRelease(Double *heatRel, Double *prodRate, Double *enthalpy);
  int CompDetailedProdRate(int i, Double *prodRate, Double *reactionRate,
                           TReactionPtr reaction);
  void ReadLewisNumbers(const char *file, VectorPtr Le);
  MatrixPtr GetCoeffHigh(void) { return fCoeffHigh; };
  MatrixPtr GetCoeffLow(void) { return fCoeffLow; };
  void ComputeTheProductionRates(Double *productionRate, Double *reactionRate,
                                 Double temp, Double pressure, Double *c,
                                 Double *k, Double *M);
  Double **GetSqrtMjOverMi(void) { return sqrtMjOverMi; };
  Double **GetSqrtMiOverMjPlusOne(void) { return sqrtMiOverMjPlusOne; };
  Flag ComputeCP_Enth(Double temp, Double *cp, Double *enth);

 protected:
  void InitSpecies(TInputDataPtr input, TReactionPtr reaction);
  void CalcUsedReactions(TReactionPtr reaction);
  void FillSpecies(TInputDataPtr input, TReactionPtr reaction);

  void ShowUsedData(const std::size_t i, const std::size_t nOfAtoms);
  void ShowSpecies(TFlameNodePtr flameNode);

  Double CompDeltaI(int i, int nSpecies, Double *M, Double *mu, Double *Y);
  Double CompDeltaIOpt(int i, int nSpeciesInSystem, Double **GijOverWj,
                       Double *M, Double *Y);
  Flag EmptyLine(char *s);

  int fNSteadyStates;              //	number of steady state species
  IntVectorPtr fNOfUsedReactions;  //  array of length nOfSpecies
  IntVectorPtr *fUsedReactions;    //  last element is
  //  fUsedReactions[nOfSpecies]->vec[fNOfUsedReactions]
  VectorPtr *fNu;  //  last element is fNu[nOfSpecies]->vec[fNOfUsedReactions]
                   //	nu of educts are > 0
  Flag fConstLewisNumber;  //
  char *fLewisNumberFile;
  VectorPtr fLewisNumber;  //	contains constant Lewis number for every species
                           //	initialized in TSpecies::FillSpecies
  Flag *fIsSteadyState;    //	signals wether the species is steady state or
                           // not
  //	length is nOfSpecies

  // constants needed for the calculation of the properties
  char **fNames;
  MatrixPtr fCoeffLow;        // last element is fCoeffLow[nOfSpecies][7]
  MatrixPtr fCoeffHigh;       // last element is fCoeffHigh[nOfSpecies][7]
  MatrixPtr fHCoeffLow;       // last element is fHCoeffLow[nOfSpecies][6]
  MatrixPtr fHCoeffHigh;      // last element is fHCoeffHigh[nOfSpecies][6]
  VectorPtr fMolarMass;       //  [kg / kmol]
  IntMatrixPtr fComposition;  // fComposition[species][atoms], last element is
  // fComposition[nOfSpecies][TInputData::fCounter->atoms]

  VectorPtr fSigma;  //	length is nOfSpecies

  MatrixPtr fInvDij;  // last element is fInvDij->mat[nOfSpecies][nOfSpecies]
  Double **sqrtMjOverMi;         //	contains sqrt( Mj/Mi ); last element is
                                 // sqrtMjOverMi[nOfSpecies][nOfSpecies]
  Double **sqrtMiOverMjPlusOne;  //	contains 0.3535534 / sqrt( Mi/Mj + 1 );
                                 // last element is
  // sqrtMiOverMjPlusOne[nOfSpecies][nOfSpecies]
  Double fTThermoLimit;
  const unsigned int fNOfSpecies;
  TransportModel TM;
};

template <class TransportModel>
class T0DSpecies : public TSpecies<TransportModel> {
 public:
  using TSpecies<TransportModel>::fMolarMass;
  using TSpecies<TransportModel>::fCoeffHigh;
  using TSpecies<TransportModel>::fHCoeffHigh;
  using TSpecies<TransportModel>::fCoeffLow;
  using TSpecies<TransportModel>::fHCoeffLow;
  using TSpecies<TransportModel>::fTThermoLimit;
  using TSpecies<TransportModel>::GetNSpeciesInSystem;
  using TSpecies<TransportModel>::GetNOfSpecies;
  using TSpecies<TransportModel>::CompDeltaIOpt;
  using TSpecies<TransportModel>::IsConstantLewisNumber;
  using TSpecies<TransportModel>::fUsedReactions;
  using TSpecies<TransportModel>::fIsSteadyState;
  using TSpecies<TransportModel>::fLewisNumber;
  using TSpecies<TransportModel>::TM;
  using TSpecies<TransportModel>::fNames;
  using TSpecies<TransportModel>::fNOfUsedReactions;
  using TSpecies<TransportModel>::CompProdCons;
  using TSpecies<TransportModel>::FindSpecies;
  using TSpecies<TransportModel>::GetNames;
  using TSpecies<TransportModel>::sqrtMiOverMjPlusOne;
  using TSpecies<TransportModel>::fNu;
  using TSpecies<TransportModel>::fInvDij;
  using TSpecies<TransportModel>::fSigma;
  using TSpecies<TransportModel>::sqrtMjOverMi;
  T0DSpecies(TInputDataPtr input, TReactionPtr reaction)
      : TSpecies<TransportModel>(input, reaction) {
    InitT0DSpecies(input);
  };
  ~T0DSpecies(void);

  VectorPtr GetHeatCapacity(void) { return fHeatCapacity; };
  VectorPtr GetConductivity(void) { return fConductivity; };
  VectorPtr GetViscosity(void) { return fViscosity; };
  VectorPtr GetEnthalpy(void) { return fEnthalpy; };
  VectorPtr GetProductionRate(void) { return fProductionRate; };
  VectorPtr GetDeltaI(void) { return fDeltaI; };
  TensorPtr GetOmegaDOverDCoeff(void) { return fOmegaDOverDCoeff; };
  Double &GetCurrTemp(void) { return fCurrTemp; };

  Flag ComputeSpeciesProperties(Double temp);
  template <class Flame>
  void PrintSpecies(Flame *flame);
  template <class Flame>
  void PrintSpecies(int number, Flame *flame, FILE *fp);
  void ComputeDiffusivityTrans1D(Double *diffusivityTrans1D_k, int GPCurr,
                                 Double temp, Double *Y, Double mixMolMass,
                                 int NOSpecies, Double *molarMass,
                                 Double pressure);

 private:
  void InitT0DSpecies(TInputDataPtr input);
  void ComputeSpeciesProperties(int number, Double temp);

#ifdef TESTING

  // friend declaration of test class and the tests that test the private member
  // void ComputeSpeciesProperties( int number, Double temp );

  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesN2StandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesHStandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesOHStandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesOStandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesHO2StandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCOStandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesHCOStandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCH3StandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCH3OStandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesC2H5StandardTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCH2OHStandardTemp);

  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesN2HighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesHHighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesOHHighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesOHighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesHO2HighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCOHighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesHCOHighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCH3HighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCH3OHighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesC2H5HighTemp);
  template <class TransModel> FRIEND_TEST(T0DSpeciesTest, ComputeSpeciesPropertiesCH2OHHighTemp);
#endif

  VectorPtr fHeatCapacity;    // [J / (kg*K)]	//	last element is
                              // fHeatCapacity[nOfSpeciesIn]
  VectorPtr fEnthalpy;        // [J / kg]			//	last element is
                              // fEnthalpy[nOfSpeciesIn]
  VectorPtr fViscosity;       // [J / (kg*K)]	//	last element is
                              // fViscosity[nOfSpeciesIn]
  VectorPtr fConductivity;    // [J / kg]			//	last element is
                              // fConductivity[nOfSpeciesIn]
  VectorPtr fProductionRate;  //  m_i	//	last element is
                              //  fProductionRate[nOfSpeciesIn]
  VectorPtr fDeltaI;          //  used by TProperties::CompMixtureProps and
                              //  TSpecies::Compute_DTherm
  Double fCurrTemp;
  TensorPtr fOmegaDOverDCoeff;
};

template <class TransportModel>
class T1DSpecies : public TSpecies<TransportModel> {
  // steady state species have no storage for fNOfUsedReactions,fUsedReactions,
  // fNu, fProductionRate
 public:
  using TSpecies<TransportModel>::fMolarMass;
  using TSpecies<TransportModel>::fCoeffHigh;
  using TSpecies<TransportModel>::fHCoeffHigh;
  using TSpecies<TransportModel>::fCoeffLow;
  using TSpecies<TransportModel>::fHCoeffLow;
  using TSpecies<TransportModel>::fTThermoLimit;
  using TSpecies<TransportModel>::GetNSpeciesInSystem;
  using TSpecies<TransportModel>::GetNOfSpecies;
  using TSpecies<TransportModel>::CompDeltaIOpt;
  using TSpecies<TransportModel>::IsConstantLewisNumber;
  using TSpecies<TransportModel>::fUsedReactions;
  using TSpecies<TransportModel>::fIsSteadyState;
  using TSpecies<TransportModel>::fLewisNumber;
  using TSpecies<TransportModel>::TM;
  using TSpecies<TransportModel>::fNames;
  using TSpecies<TransportModel>::fNOfUsedReactions;
  using TSpecies<TransportModel>::CompProdCons;
  using TSpecies<TransportModel>::FindSpecies;
  using TSpecies<TransportModel>::GetNames;
  using TSpecies<TransportModel>::sqrtMiOverMjPlusOne;
  using TSpecies<TransportModel>::fNu;
  using TSpecies<TransportModel>::fInvDij;
  using TSpecies<TransportModel>::fSigma;
  using TSpecies<TransportModel>::sqrtMjOverMi;
  using TSpecies<TransportModel>::ShowSpecies;
  T1DSpecies(TInputDataPtr input, T1DReactionPtr reaction)
      : TSpecies<TransportModel>(input, reaction) {
    InitT1DSpecies(input);
  };
  ~T1DSpecies(void);

  TensorPtr GetBinDij(void) { return fBinDij; };
  TensorPtr GetGijOverWj(void) { return fGijOverWj; };
  TensorPtr GetOmegaDOverDCoeff(void) { return fOmegaDOverDCoeff; };
  TensorPtr GetDThermConst(void) { return fDThermConst; };
  MatrixPtr GetProductionRate(void) { return fProductionRate; };
  MatrixPtr GetViscosity(void) { return fViscosity; };
  MatrixPtr GetHeatCapacity(void) { return fHeatCapacity; };
  MatrixPtr GetConductivity(void) { return fConductivity; };
  MatrixPtr GetEnthalpy(void) { return fEnthalpy; };
  MatrixPtr GetDiffusivity(void) { return fDiffusivity; };
  MatrixPtr GetDiffTherm(void) { return fThermoDiffusivity; };
  MatrixPtr GetDeltaI(void) { return fDeltaI; };
  void PrintSpecies(int k);                                // prints all
  void PrintSpecies(int number, int gridPoint, FILE *fp);  // prints one
  void PrintProductionRate(TNewtonPtr bt);
  template <class Flame>
  void PrintProdCons(TNewtonPtr bt, Flame *flame);
  void PrintDiffusionCoeffs(TNewtonPtr bt);
  template <class Flame>
  void PrintProdRateTerms(const char *name, Flame *flame);
  template <class Flame>
  void PrintProdRateTerms(int i, Flame *flame);
  Flag ComputeSpeciesProperties(
      TFlameNodePtr flameNode, Double pressure,
      Double temp);  // computes properties of all species
  void ComputeSpeciesProperties(
      TFlameNodePtr flameNode, Double temp, Double pressure,
      int number);  // computes properties of one species
  void Compute_D(TFlameNodePtr flameNode, Double temp, Double *Y,
                 Double pressure,
                 Flag newTemp);  // computes diffusion coefficients
  void Compute_D(TFlameNodePtr flameNode);
  void ComputeDeltaIOpt(TFlameNodePtr flameNode, Double *Y, Double **GijOverWj,
                        Flag newTemp);
  void CompBinDiffCoeff(TFlameNodePtr flameNode, Double temp, Double *Y,
                        Double press);
  void Compute_DTherm(TFlameNodePtr flameNode, Flag calcNewConst);
  Double ComputeOneDiffusionCoeff(int number, Double *Y, Double *invDij,
                                  Double mixMolarMass);
  MatrixPtr GetSavedY(void) { return fSavedY; };
  MatrixPtr GetSavedDeltaiY(void) { return fSavedDeltaiY; };
  MatrixPtr GetSumDiff(void) { return fSumDiff; };
  VectorPtr GetTempProp(void) { return fTempProp; };
  VectorPtr GetPressureProp(void) { return fPressureProp; };

 private:
  void InitT1DSpecies(TInputDataPtr input);

  MatrixPtr fViscosity;          // [kg / (m*s)]		//	last element is
                                 // fViscosity[nGridPoints][nOfSpecies]
  MatrixPtr fHeatCapacity;       // [J / (kg*K)]	//	last element is
                                 // fHeatCapacity[nGridPoints][nOfSpecies]
  MatrixPtr fConductivity;       // [W / (m*K)]	//	last element is
                                 // fConductivity[nGridPoints][nOfSpecies]
  MatrixPtr fEnthalpy;           // [J / kg]			//	last element is
                                 // fEnthalpy[nGridPoints][nOfSpecies]
  MatrixPtr fDiffusivity;        // [m^2 / s]		//	last element is
                                 // fDiffusivity[nGridPoints][nOfSpecies]
  MatrixPtr fThermoDiffusivity;  // [m^2 / s]//	last element is
                                 // fDiffusivity[nGridPoints][nOfSpecies]
  MatrixPtr fProductionRate;     //  m_i	//	last element is
                                 //  fViscosity[nGridPoints][nOfSpecies -
                                 //  nSteadyStates]
  TensorPtr fBinDij;             // last element is
  // fBinDij->tensor[nGridPoints][nOfSpeciesIn][nOfSpeciesIn]
  TensorPtr fGijOverWj;  //	contains 0.3535534 / sqrt( Mi/Mj + 1 ); last
                         // element is
  // sqrtMiOverMjPlusOne[nOfSpecies][nOfSpecies]
  TensorPtr fOmegaDOverDCoeff;
  TensorPtr fDThermConst;

  MatrixPtr fSavedY;
  MatrixPtr fSavedDeltaiY;
  MatrixPtr fSumDiff;
  MatrixPtr fDeltaI;  //  used by TProperties::CompMixtureProps and
                      //  TSpecies::Compute_DTherm
  VectorPtr fTempProp;
  VectorPtr fPressureProp;
};

#include "TSpecies.hpp"

#endif  // TSPECIES__
