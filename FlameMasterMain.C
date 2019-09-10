#include "FlameMaster.h"

template <class TransModel>
void chooseFlame(FirstInput& finput) {
    cout << "The solver is " << kStartUpDiffusion << endl;
  if (finput.fFlameType == kStartUpDiffusion) {
    TFlameSheet<T1DSpecies<TransModel> > flame(finput);
    flame.GetInputData()->PrintAdditionalData();
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kDiffusionPhys) {
    TCountDiffFlamePhys<T1DSpecies<TransModel> > flame(finput);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kDiffPhysEigen) {
    TCountDiffPhysEigen<T1DSpecies<TransModel> > flame(finput);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kUnstrPremPhys) {
    TUnstrPremFlamePhys<T1DSpecies<TransModel> > flame(finput);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kCountPremPhys) {
    TCountPremFlamePhys<T1DSpecies<TransModel> > flame(finput);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kCountPremSim) {
    TCountPremFlameSim<T1DSpecies<TransModel> > flame(finput);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kStartUpDiffusion2) {
    TCountDiffFlameSim<T1DSpecies<TransModel> > flame(finput);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kCountDiffCont) {
    // rlanger
    std::cerr << "\n\n\tWARNING: kCountDiffCont is not tested! This was "
                 "previously commented out!\n\n\n";
    //   finput.fFlameType = kStartUpDiffusion2;
    finput.fVariablesWithoutSpecies = 4;
    TCountDiffFlameContPre<T1DSpecies<TransModel> > flamePre(finput);
    flamePre.GetSolver()->Solve(&flamePre);
    Double dTds = flamePre.GetdTds();
    Double dAInvds = flamePre.GetdAInvds();
    Double strainRate = flamePre.GetStrainRate();
    TSolutionPtr sol = flamePre.GetSolution();
    //   finput.fFlameType = kCountDiffCont;
    finput.fVariablesWithoutSpecies = 5;
    TCountDiffFlameCont<T1DSpecies<TransModel> > flame(finput, sol, strainRate);
    delete sol;
    flame.SetdTds(dTds);
    flame.SetdAInvds(dAInvds);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kCountDiffMix) {
    TCountDiffFlameMix<T1DSpecies<TransModel> > flame(finput);
    flame.GetSolver()->Solve(&flame);
  } else if (finput.fFlameType == kHomoIsoChor ||
             finput.fFlameType == kHomoIsoBar) {
    T0DIsoChor<T0DSpecies<TransModel> > flame(finput);
    flame.Solve();
  } else if (finput.fFlameType == kHomoPSR) {
    T0DPSR<T0DSpecies<TransModel> > flame(finput);
    flame.Solve();
  } else if (finput.fFlameType == kTransFlamelet) {
    TTransFlamelet<T0DSpecies<TransModel> > flame(finput);
    flame.Solve();
  }
#ifdef COMPILE_FORTRAN_SRC
  else if (finput.fFlameType == kTrans1DIsoChor) {
    std::cerr << "\n\n\tWARNING: kTrans1DIsoChor is not tested! This was "
                 "previously commented out!\n\n\n";
    TTrans1DIsoChor<T0DSpecies<TransModel> > flame(finput);
    flame.Solve();
  } else if (finput.fFlameType == kEquilibrium) {
    std::cerr << "\n\n\tWARNING: kEquilibrium is not tested! This was "
                 "previously commented out!\n\n\n";
    TEquilibrium<T0DSpecies<TransModel> > flame(finput);
    flame.Solve();
  }
#endif  // COMPILE_FORTRAN_SRC
  else {
    std::cerr << "error: wrong flame type specified: " << finput.fFlameType
              << NEWL;
    exit(2);
  }
}

void run(FirstInput&& finput) {
  // pick appropriate transport model
  switch (finput.TMT) {
    case TransportModelType::MONO_ATOMIC:
      chooseFlame<TMonoAtomicTransportModel>(finput);
      break;
    case TransportModelType::MULTI_GEOMETRY:
      chooseFlame<TMultiGeometryTransportModel>(finput);
      break;
    default:
      std::cerr << "#error: unknown transport model (choose \"MonoAtomic\" "
                   "or \"MultiGeometry\" in your input)";
      exit(2);
  }
}



int main(int argc, char* argv[]) {
  run(FirstInput(argc, argv));
  fputc('\a', stderr);
  return 0;
}
