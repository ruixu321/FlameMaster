#undef FULLDIFFUSION

// following makes use of liquid fuel easier, where the assumption of adiabatic
// droplet evaporation can lead to negative temperatures. Keeping liquid cp
// constant leads to smaller temperature change during evaporation
#define CONSTLOWTCP

// if source code for production rates is generated with ScanMan, it can be
// tested here and compared with the solution of the mechanism itself. For using
// production rate files produced by ScanMan, define PRODRATEFILE. If Fortran
// files should be used, define PRODRATEFILEF77, for C files undefine
// PRODRATEFILEF77
// Needs to be set also in TFlame.C
#undef PRODRATEFILE
#undef PRODRATEFILEF77

// optimization
#define OPTPROPS
#define OPTIMIZEDELTAI
#define OPTOMEGAD
#define OPTDTHERM

#ifdef PRODRATEFILE
void ComputeThermoData(Double *h, Double *cp, Double T);
#endif

#ifdef PRODRATEFILEF77
#define PRODRATES prodrates_
#define COMPTHERMODATA compthermodata_
#define GETMOLARMASS getmolarmass_
extern "C" {
void PRODRATES(Double *cdot, Double *w, Double *k, Double *c, Double *M,
               Double *temp, Double *pressure);
void COMPTHERMODATA(Double *h, Double *cp, Double *T);
void GETMOLARMASS(Double *MM);
}

#endif

template <class TransportModel>
void TSpecies<TransportModel>::InitSpecies(TInputDataPtr input,
                                           TReactionPtr reaction) {
  int i;
  int nOfSpecies = input->GetCounter()->species;
  int nOfReactions = input->GetCounter()->reactions;

  fNSteadyStates = input->GetCounter()->steadyStates;

  fTThermoLimit = 273.15;

  fNOfUsedReactions = NewIntVector(nOfSpecies);
  CalcUsedReactions(reaction);
  fUsedReactions = new IntVectorPtr[nOfSpecies];
  if (!fUsedReactions) FatalError("memory allocation of TSpecies failed");
  fNu = new VectorPtr[nOfSpecies];
  if (!fNu) FatalError("memory allocation of TSpecies failed");
  for (i = 0; i < nOfSpecies; ++i) {
    if (fNOfUsedReactions->vec[i]) {
      fUsedReactions[i] =
          NewIntVector(fNOfUsedReactions->vec[i]);  // last element is
      // fUsedReactions[nOfSpecies]->vec[fNOfUsedReactions->vec[i]]
      fNu[i] = NewVector(fNOfUsedReactions->vec[i]);  // last element is
      // fNu[nOfSpecies]->vec[fNOfUsedReactions->vec[i]]
    } else {
      fUsedReactions[i] = NewIntVector(1);
      fNu[i] = NewVector(1);
      fUsedReactions[i]->len = 0;
      fNu[i]->len = 0;
    }
  }

  fLewisNumberFile = input->fLewisNumberFile;
  fLewisNumber = NewVector(nOfSpecies);

  fNames = new String[nOfSpecies];
  if (!fNames) FatalError("memory allocation of TSpecies failed");
  fCoeffLow =
      NewMatrix(NOFCOEFFS, nOfSpecies,
                kColumnPointers);  // last element is fCoeffLow[nOfSpecies][7]
  fCoeffHigh =
      NewMatrix(NOFCOEFFS, nOfSpecies,
                kColumnPointers);  // last element is fCoeffHigh[nOfSpecies][7]
  fHCoeffLow =
      NewMatrix(NOFCOEFFS - 1, nOfSpecies,
                kColumnPointers);  // last element is fCoeffLow[nOfSpecies][7]
  fHCoeffHigh =
      NewMatrix(NOFCOEFFS - 1, nOfSpecies,
                kColumnPointers);  // last element is fCoeffHigh[nOfSpecies][7]
  fMolarMass = NewVector(nOfSpecies);
  fSigma = NewVector(nOfSpecies);
  fInvDij = NewMatrix(nOfSpecies, nOfSpecies, kColumnPointers);

  fComposition = _NewIntMatrix(nOfSpecies, input->GetCounter()->atoms,
                               kRowPointers);  // ###

  fIsSteadyState = new Flag[nOfSpecies];

  int nSpeciesIn = GetNSpeciesInSystem();

  sqrtMjOverMi = New2DArray(nSpeciesIn, nSpeciesIn);
  sqrtMiOverMjPlusOne = New2DArray(nSpeciesIn, nSpeciesIn);

  FillSpecies(input, reaction);
}

template <class TransportModel>
TSpecies<TransportModel>::~TSpecies(void) {
  int j;
  int nOfSpecies = GetNOfSpecies();

  Free2DArray(sqrtMiOverMjPlusOne);
  Free2DArray(sqrtMjOverMi);

  delete fIsSteadyState;

  for (j = 0; j < nOfSpecies; ++j) {
    delete fNames[j];
  }

  DisposeMatrix(fInvDij);
  DisposeVector(fSigma);
  DisposeVector(fMolarMass);
  DisposeMatrix(fHCoeffHigh);
  DisposeMatrix(fHCoeffLow);
  DisposeMatrix(fCoeffHigh);
  DisposeMatrix(fCoeffLow);
  delete fNames;
  DisposeVector(fLewisNumber);

  for (j = 0; j < nOfSpecies; ++j) {
    DisposeVector(fNu[j]);
    DisposeIntVector(fUsedReactions[j]);
  }
  delete fNu;
  delete fUsedReactions;

  DisposeIntVector(fNOfUsedReactions);

  _DisposeIntMatrix(fComposition);
}

template <class TransportModel>
void TSpecies<TransportModel>::FillSpecies(TInputDataPtr input,
                                           TReactionPtr reaction) {
  int i, j, k, l;
  int nOfAtoms = input->GetCounter()->atoms;
  int nOfSpecies = GetNOfSpecies();
  int nOfReactions = reaction->GetNOfReactions();
  VectorPtr *nuReaction = reaction->GetNu();
  IntVectorPtr *speciesNumber = reaction->GetSpeciesNumber();

  for (i = 0; i < nOfSpecies; ++i) {  // i loops species
    l = 0;
    for (j = 0; j < nOfReactions; ++j) {  // j loops reactions
      for (k = 0; k < speciesNumber[j]->len;
           ++k) {  // k loops reactions speciesvector
        if (i ==
            speciesNumber[j]->vec[k]) {  // then species i appears in reaction j
          fUsedReactions[i]->vec[l] = j;
          fNu[i]->vec[l] = nuReaction[j]->vec[k];
          ++l;
        }
      }
    }
  }

  for (i = 0; i < nOfSpecies; ++i) {  // i loops species
    for (j = 0; j < nOfAtoms; ++j) {  // j loops atoms
      fComposition->mat[i][j] =
          input->GetReadSpeciesData()->GetComposition(i).at(j);
    }
  }

  Double *coeffLow;
  Double *coeffHigh;
  for (i = 0; i < nOfSpecies; ++i) {
    // names
    fNames[i] =
        new char[strlen(input->GetReadSpeciesData()->GetName(i).c_str()) + 1];
    if (!fNames[i]) FatalError("memory allocation of TSpecies failed");
    strcpy(fNames[i], input->GetReadSpeciesData()->GetName(i).c_str());
    coeffLow = fCoeffLow->mat[i];
    coeffHigh = fCoeffHigh->mat[i];
    for (j = 0; j < NOFCOEFFS; ++j) {
      coeffLow[j] = input->GetReadSpeciesData()->GetSpecHeatLow(i).at(j);
      coeffHigh[j] = input->GetReadSpeciesData()->GetSpecHeatHigh(i).at(j);
    }

    coeffLow = fHCoeffLow->mat[i];
    coeffHigh = fHCoeffHigh->mat[i];
    coeffLow[0] = input->GetReadSpeciesData()->GetSpecHeatLow(i).at(0);
    coeffHigh[0] = input->GetReadSpeciesData()->GetSpecHeatHigh(i).at(0);
    coeffLow[1] = input->GetReadSpeciesData()->GetSpecHeatLow(i).at(1) * 0.5;
    coeffHigh[1] = input->GetReadSpeciesData()->GetSpecHeatHigh(i).at(1) * 0.5;
    coeffLow[2] = input->GetReadSpeciesData()->GetSpecHeatLow(i).at(2) / 3.0;
    coeffHigh[2] = input->GetReadSpeciesData()->GetSpecHeatHigh(i).at(2) / 3.0;
    coeffLow[3] = input->GetReadSpeciesData()->GetSpecHeatLow(i).at(3) * 0.25;
    coeffHigh[3] = input->GetReadSpeciesData()->GetSpecHeatHigh(i).at(3) * 0.25;
    coeffLow[4] = input->GetReadSpeciesData()->GetSpecHeatLow(i).at(4) * 0.2;
    coeffHigh[4] = input->GetReadSpeciesData()->GetSpecHeatHigh(i).at(4) * 0.2;
    coeffLow[5] = input->GetReadSpeciesData()->GetSpecHeatLow(i).at(5);
    coeffHigh[5] = input->GetReadSpeciesData()->GetSpecHeatHigh(i).at(5);

#ifndef PRODRATEFILEF77
    fMolarMass->vec[i] = input->GetReadSpeciesData()->GetMolarMass(i);
#endif

    fSigma->vec[i] = input->GetReadSpeciesData()->GetCollisionDiam(i);
    if (input->GetReadSpeciesData()->IsSteadyState(i)) {
      fIsSteadyState[i] = TRUE;
    } else {
      fIsSteadyState[i] = FALSE;
    }
  }

#ifdef PRODRATEFILEF77
  GETMOLARMASS(fMolarMass->vec);
#endif

  Double *M = fMolarMass->vec;
  int nSpeciesIn = GetNSpeciesInSystem();

  for (i = 0; i < nSpeciesIn; ++i) {
    for (j = 0; j < nSpeciesIn; ++j) {
      sqrtMjOverMi[j][i] = sqrt(M[j] / M[i]);
      sqrtMiOverMjPlusOne[i][j] = 0.3535534 / sqrt(M[i] / M[j] + 1.0);
    }
  }

  fConstLewisNumber = input->fConstantLewisNumber;

  //for(std::size_t i = 0; i < nOfSpecies; ++ i) ShowUsedData(i, input->GetCounter()->atoms);
  //for(std::size_t i = 0; i < nOfSpecies; ++ i) TM.ShowUsedData(i);
}

template<class TransportModel>
void TSpecies<TransportModel>::ShowUsedData(const std::size_t i, const std::size_t nOfAtoms){
  std::cout.precision(dbl::max_digits10);
  std::cout << "####### " << fNames[i] << "[" << i << "] #######\n";
  std::cout << "\tComposition: ";
  for(auto j = 0; j < nOfAtoms; ++j) std::cout << fComposition->mat[i][j] << ' ';
  std::cout << '\n';
  std::cout << "\tNASA, cp low: ";
  for(auto j = 0; j<NOFCOEFFS; ++j) std::cout << fCoeffLow->mat[i][j] << ' ';
  std::cout << '\n';
  std::cout << "\tNASA, cp high: ";
  for(auto j = 0; j<NOFCOEFFS; ++j) std::cout << fCoeffHigh->mat[i][j] << ' ';
  std::cout << '\n';
  std::cout << "\tNASA, h low: ";
  for(auto j = 0; j<NOFCOEFFS; ++j) std::cout << fHCoeffLow->mat[i][j] << ' ';
  std::cout << '\n';
  std::cout << "\tNASA, h high: ";
  for(auto j = 0; j<NOFCOEFFS; ++j) std::cout << fHCoeffHigh->mat[i][j] << ' ';
  std::cout << '\n';
  std::cout << "\tMolar Mass: " << fMolarMass->vec[i] << '\n';
  std::cout << "\tCollision diameter: " << fSigma->vec[i] << '\n';
    int nSpeciesIn = GetNSpeciesInSystem();
    std::cout << "\tsqrtMjOverMi[" << i << "] ";
    for (auto j = 0; j < nSpeciesIn; ++j)
    {
      std::cout << " [" << j << "]: " << sqrtMjOverMi[i][j];
    }
    std::cout << '\n';
    std::cout << "\tsqrtMiOverMjPlusOne[" << i << "] ";
    for (auto j = 0; j < nSpeciesIn; ++j)
    {
      std::cout << " [" << j << "]: " << sqrtMiOverMjPlusOne[i][j];
    }
    std::cout << '\n';
}

template <class TransportModel>
void TSpecies<TransportModel>::CalcUsedReactions(TReactionPtr reaction) {
  int i, j, k;

  for (i = 0; i < fNOfUsedReactions->len; ++i) {
    for (j = 0; j < reaction->GetNOfReactions(); ++j) {
      for (k = 0; k < reaction->GetSpeciesNumber()[j]->len; ++k) {
        if (i == reaction->GetSpeciesNumber()[j]->vec[k]) {
          ++fNOfUsedReactions->vec[i];
          continue;
        }
      }
    }
  }
}

template <class TransportModel>
Double TSpecies<TransportModel>::GetPosNegProductionRate(int which,
                                                         Double *reactionRate,
                                                         Flag pos) {
  int j;
  int nSpeciesInSystem = GetNSpeciesInSystem();

  Double productionRate = 0.0;
  Double *nu;
  Double *W = fMolarMass->vec;
  int *usedReactions;
  int nOfUsedReactions;

  nu = fNu[which]->vec;
  usedReactions = fUsedReactions[which]->vec;
  nOfUsedReactions = fNOfUsedReactions->vec[which];
  for (j = 0; j < nOfUsedReactions; ++j) {
    if (pos == TRUE && nu[j] < 0.0) {
      productionRate += nu[j] * reactionRate[usedReactions[j]];
    } else if (pos == FALSE && nu[j] > 0.0) {
      productionRate += nu[j] * reactionRate[usedReactions[j]];
    }
  }
  productionRate *= -W[which];

  return productionRate;
}

template <class TransportModel>
void TSpecies<TransportModel>::ComputeProductionRates(Double *productionRate,
                                                      Double *reactionRate) {
  int i, j;
  int nSpeciesInSystem = GetNSpeciesInSystem();

  Double *nu;
  Double *W = fMolarMass->vec;
  int *usedReactions;
  int nOfUsedReactions;

  for (i = 0; i < nSpeciesInSystem; ++i) {
    //		ComputeProductionRate( i, productionRate[i], reactionRate );
    nu = fNu[i]->vec;
    usedReactions = fUsedReactions[i]->vec;
    nOfUsedReactions = fNOfUsedReactions->vec[i];
    productionRate[i] = 0.0;
    for (j = 0; j < nOfUsedReactions; ++j) {
      productionRate[i] += nu[j] * reactionRate[usedReactions[j]];
    }
    productionRate[i] *= -W[i];
  }
}

#ifdef PRODRATEFILE
template <class TransportModel>
void ComputeProductionRates(Double *cdot, Double *w, Double *k, Double *c,
                            Double *M, Double temp, Double pressure);

template <class TransportModel>
void TSpecies::ComputeTheProductionRates(Double *productionRate,
                                         Double *reactionRate, Double temp,
                                         Double pressure, Double *c, Double *k,
                                         Double *M) {
  ::ComputeProductionRates(productionRate, reactionRate, k, c, M, temp,
                           pressure);
}
#endif

#ifdef PRODRATEFILEF77
template <class TransportModel>
void TSpecies<TransportModel>::ComputeTheProductionRates(
    Double *productionRate, Double *reactionRate, Double temp, Double pressure,
    Double *c, Double *k, Double *M) {
  ::PRODRATES(productionRate, reactionRate, k, c, M, &temp, &pressure);
}
#endif

template <class TransportModel>
void TSpecies<TransportModel>::ComputeNoProductionRates(
    Double *productionRate) {
  int nSpeciesInSystem = GetNSpeciesInSystem();

  for (int i = 0; i < nSpeciesInSystem; ++i) {
    productionRate[i] = 0.0;
  }
}

template <class TransportModel>
void TSpecies<TransportModel>::ComputeProductionRate(int i,
                                                     Double &productionRate,
                                                     Double *reactionRate) {
  productionRate = 0.0;
  Double *nu = fNu[i]->vec;
  int *usedReactions = fUsedReactions[i]->vec;
  int nOfUsedReactions = fNOfUsedReactions->vec[i];

  for (int j = 0; j < nOfUsedReactions; ++j) {
    productionRate += nu[j] * reactionRate[usedReactions[j]];
  }
  productionRate *= -fMolarMass->vec[i];
}

template <class TransportModel>
void TSpecies<TransportModel>::CompHeatRelease(Double *heatRel,
                                               Double *prodRate,
                                               Double *enthalpy) {
  *heatRel = 0.0;

  for (int i = 0; i < GetNSpeciesInSystem(); ++i) {
    *heatRel += prodRate[i] * enthalpy[i];
  }
}

template <class TransportModel>
void TSpecies<TransportModel>::CompProdCons(Double *source, Double *sink,
                                            Double *reacRate) {
  Double productionRate;
  int *nUsedReacs = fNOfUsedReactions->vec;
  int *usedReacs;
  Double *nu;
  //	Double	*molarMass = fMolarMass->vec;

  for (int i = 0; i < GetNSpeciesInSystem(); ++i) {
    usedReacs = fUsedReactions[i]->vec;
    nu = fNu[i]->vec;
    for (int j = 0; j < nUsedReacs[i]; ++j) {
      productionRate = nu[j] * reacRate[usedReacs[j]];
      if (productionRate > 0.0) {
        sink[i] -= productionRate;
      } else {
        source[i] -= productionRate;
      }
    }
    //		sink[i] *= molarMass[i];
    //		source[i] *= molarMass[i];
  }
}

template <class TransportModel>
int TSpecies<TransportModel>::CompDetailedProdRate(int i, Double *prodRate,
                                                   Double *reacRate,
                                                   TReactionPtr reaction) {
  int k;
  int entry = 0;
  int usedReac;
  int *nUsedReacs = fNOfUsedReactions->vec;
  int *forwReacs = reaction->GetForwardReacs()->vec;
  int *backReacs = reaction->GetBackwardReacs()->vec;
  Double stoiCoeffBack;
  char **labels = reaction->GetLabels();

  for (int j = 0; j < nUsedReacs[i]; ++j) {
    usedReac = fUsedReactions[i]->vec[j];
    if (forwReacs[usedReac] < 0) {  // means not a backward reaction
      prodRate[entry] = fNu[i]->vec[j] * reacRate[usedReac];
      if (backReacs[usedReac] >=
          0) {  // means there is a corresponding backward reaction
        for (k = 0; k < nUsedReacs[i]; ++k) {
          if (strcmp(labels[backReacs[usedReac]],
                     labels[fUsedReactions[i]->vec[k]]) == 0) {
            stoiCoeffBack = fNu[i]->vec[k];
            break;
          }
        }
        prodRate[entry] += stoiCoeffBack * reacRate[backReacs[usedReac]];
      }
      prodRate[entry] *= -1.0;  // fMolarMass->vec[i];
      ++entry;
    }
  }

  return entry;
}

template <class TransportModel>
Double TSpecies<TransportModel>::ComputeDCpiDT(int index, Double temp) {
  Double *coeff =
      (temp > 1000.0) ? fCoeffHigh->mat[index] : fCoeffLow->mat[index];

  return (RGAS / fMolarMass->vec[index] *
          (coeff[1] +
           temp * (2.0 * coeff[2] +
                   temp * (3.0 * coeff[3] + temp * 4.0 * coeff[4]))));
}

template <class TransportModel>
int TSpecies<TransportModel>::FindSpecies(const char *species) {
  int nOfSpecies = GetNOfSpecies();

  for (int i = 0; i < nOfSpecies; ++i) {
    if (strcmp(species, fNames[i]) == 0) {
      return i;
    }
  }
  return -1;
}

template <class TransportModel>
void TSpecies<TransportModel>::WriteRoggsSymbolsData(void) {
  // inerte Spezies muessen am Ende der Liste stehen

  FILE *fp = fopen("roggssymbols.data", "w");
  if (!fp) {
    cerr << "#warning: unable to open file 'roggssymbols.data'" << NEWL;
    return;
  }
  fprintf(fp, "LIST OF SYM / CHNO..........................................\n");
  for (int i = 0; i < GetNOfSpecies(); ++i) {
    fprintf(fp, "%-5s%-8s 0000\n", "Y", fNames[i]);
  }
  fprintf(fp, "-END OF SYM   0000\n");
  fprintf(fp, "BOUNDS ON TEMPERATURE.............\n");
  fprintf(fp, "0298.00    2500.00\n");
  fprintf(fp, "-END OF BOUNDS\n");
  fprintf(fp, "CHEMTP: TYPE OF CHEMISTRY.........\n");
  fprintf(fp, "DETAILED CHEMISTRY\n");
  fprintf(fp, "-END OF CHEMTP\n");
  fprintf(fp, "NSTMEC: PARAMETER.................\n");
  fprintf(fp, " -1\n");
  fprintf(fp, "-END OF NSTMEC\n");

  fclose(fp);
}

#ifdef OPTIMIZEDELTAI
template <class TransportModel>
void T1DSpecies<TransportModel>::ComputeDeltaIOpt(TFlameNodePtr flameNode,
                                                  Double *Y, Double **GijOverWj,
                                                  Flag newTemp) {
  int i;
  int nSpeciesIn = GetNSpeciesInSystem();
  Double *deltaI = flameNode->deltaI;
  Double *M = fMolarMass->vec;

#define NEWOPTDELTAI
#ifdef NEWOPTDELTAI
  Double *savedDeltaiY = flameNode->savedDeltaiY;
  int changed = -1;  // -1 means nothing changed, -2 means temperature
  // or more than two Y_i changed, >=0 is changed species, if only
  // one or two changed
  int changed2 = -1;  // -1 means nothing changed, -2 means temperature

  if (newTemp) {
    changed = -2;
  } else {
    for (i = 0; i < nSpeciesIn; ++i) {
      if (Y[i] != savedDeltaiY[i]) {
        if (changed2 < 0) {
          if (changed == -1) {
            changed = i;
          } else {
            changed2 = i;
          }
        } else {
          changed = -2;
          break;
        }
      }
    }
  }

  if (changed >= 0 /*&& changed2 == -1*/) {
    for (i = 0; i < nSpeciesIn; ++i) {
      deltaI[i] =
          deltaI[i] +
          M[i] * GijOverWj[i][changed] * (Y[changed] - savedDeltaiY[changed]);
    }
    if (changed2 >= 0) {
      for (i = 0; i < nSpeciesIn; ++i) {
        deltaI[i] = deltaI[i] +
                    M[i] * GijOverWj[i][changed2] *
                        (Y[changed2] - savedDeltaiY[changed2]);
      }
    }
    savedDeltaiY[changed] = Y[changed];
    if (changed2 >= 0) {
      savedDeltaiY[changed2] = Y[changed2];
    }
  } else if (changed == -2) {
    for (i = 0; i < nSpeciesIn; ++i) {
      savedDeltaiY[i] = Y[i];
      deltaI[i] = CompDeltaIOpt(i, nSpeciesIn, GijOverWj, M, Y);
    }
  } else {
    /*		fprintf( stderr, "changed = %d\tchanged2 = %d\n", changed,
     * changed2
     * );*/
    /*		Double	check = 0.0;*/
    /*		for ( i = 0; i < nSpeciesIn; ++i ) {*/
    /*			check = CompDeltaIOpt( i, nSpeciesIn, GijOverWj, M, Y
     * );*/
    /*		//	fprintf( stderr, "deltaI[%d] = %g\tcheck = %g\n", i,
     * deltaI[i],
     * check );*/
    /*		}*/
  }
#else
  for (i = 0; i < nSpeciesIn; ++i) {
    deltaI[i] = CompDeltaIOpt(i, nSpeciesIn, GijOverWj, M, Y);
  }
#endif
}

template <class TransportModel>
Double TSpecies<TransportModel>::CompDeltaIOpt(int i, int nSpeciesInSystem,
                                               Double **GijOverWj, Double *M,
                                               Double *Y) {
  Double Delta_i = 0;

  /*	Double	*loci = GijOverWj[i];*/
  /*	for ( int j = 0; j < nSpeciesInSystem; ++j ){*/
  /*		Delta_i += loci[j] * Y[j];*/
  /*	}*/
  /*	*/
  /*	return Delta_i * M[i];*/

  Delta_i = DotProd(nSpeciesInSystem, GijOverWj[i], 1, Y, 1);

  return Delta_i * M[i];
}
#endif

template <class TransportModel>
void TSpecies<TransportModel>::ComputeDeltaI(Double *deltaI, Double *Y,
                                             Double *viscosity) {
  int i;
  int nSpeciesIn = GetNSpeciesInSystem();
  Double *M = fMolarMass->vec;

  for (i = 0; i < nSpeciesIn; ++i) {
    deltaI[i] = CompDeltaI(i, nSpeciesIn, M, viscosity, Y);
  }
}

template <class TransportModel>
Double TSpecies<TransportModel>::CompDeltaI(int i, int nSpeciesInSystem,
                                            Double *M, Double *mu, Double *Y) {
  Double G_ij, Delta_i = 0;

  for (int j = 0; j < nSpeciesInSystem; ++j) {
    G_ij = sqrt(sqrtMjOverMi[j][i] * mu[i] / mu[j]) + 1.0;
    //	(old form before bug was fixed)	G_ij = sqrt( sqrtMjOverMi[j][i] * mu[j]
    /// mu[i] ) + 1.0;
    G_ij *= G_ij * sqrtMiOverMjPlusOne[i][j];
    Delta_i += G_ij / M[j] * Y[j];
  }

  return Delta_i * M[i];
}

template <class TransportModel>
void TSpecies<TransportModel>::ReadLewisNumbers(const char *file,
                                                VectorPtr Le) {
  char inSpecies[32];
  Double inLe;
  char buffer[80];
  int i, empty = 0;
  FILE *fp = NULL;
  Double *lewis = Le->vec;

  for (i = 0; i < Le->len; ++i) {
    lewis[i] = 1.0;
  }

  if (!(fp = fopen(file, "r"))) {
    fprintf(stderr,
            "Couldn't open \"%s\"---all Lewis numbers will be set to 1.0\n",
            file);
    return;
  }

  while (fgets(buffer, 79, fp) != NULL) {
    if (EmptyLine(buffer)) {
      ++empty;                     /*	Skip empty lines.				*/
    } else if (buffer[0] == '#') { /*	Skip line comments. */
    } else if (sscanf(buffer, "%s%lf", inSpecies, &inLe) == 2) {
      if ((i = FindSpecies(inSpecies)) < 0) {
      } else
        lewis[i] = inLe;
    }
  }
}

template <class TransportModel>
Flag TSpecies<TransportModel>::EmptyLine(char *s) {
  while (*s) {
    if (!isspace(*s)) return FALSE;
    ++s;
  }

  return TRUE;
}

template <class TransportModel>
Flag TSpecies<TransportModel>::ComputeCP_Enth(Double temp, Double *cp,
                                              Double *enth) {
  int number;
  int nSpeciesInSystem = GetNSpeciesInSystem();
  Double *coeff, *hCoeff; /* pointer to nasa coefficients */
  Double r_over_m;        /* in J / (kg K) */

  for (number = 0; number < nSpeciesInSystem; ++number) {
    r_over_m = RGAS / fMolarMass->vec[number];
    if (temp > 1000.0) {
      coeff = fCoeffHigh->mat[number];
      hCoeff = fHCoeffHigh->mat[number];
    } else {
      coeff = fCoeffLow->mat[number];
      hCoeff = fHCoeffLow->mat[number];
    }
    if (temp < fTThermoLimit) {
      cp[number] =
          r_over_m *
          (coeff[0] +
           fTThermoLimit *
               (coeff[1] +
                fTThermoLimit *
                    (coeff[2] +
                     fTThermoLimit * (coeff[3] + fTThermoLimit * coeff[4]))));

      enth[number] =
          r_over_m *
          (fTThermoLimit *
               (hCoeff[0] +
                fTThermoLimit *
                    (hCoeff[1] +
                     fTThermoLimit *
                         (hCoeff[2] +
                          fTThermoLimit *
                              (hCoeff[3] + fTThermoLimit * hCoeff[4])))) +
           hCoeff[5]);
#ifdef CONSTLOWTCP
      enth[number] += (temp - fTThermoLimit) * cp[number];
#else
      Double B =
          r_over_m *
          (coeff[1] +
           fTThermoLimit * (2.0 * coeff[2] +
                            fTThermoLimit * (3.0 * coeff[3] +
                                             4.0 * fTThermoLimit * coeff[4])));
      Double cpRef = cp[number];
      Double A = cpRef - B * fTThermoLimit;
      cp[number] = A + B * temp;
      enth[number] += 0.5 * (cp[number] + cpRef) * (temp - fTThermoLimit);
#endif
    } else {
      cp[number] =
          r_over_m *
          (coeff[0] +
           temp * (coeff[1] +
                   temp * (coeff[2] + temp * (coeff[3] + temp * coeff[4]))));

      enth[number] =
          r_over_m *
          (temp * (hCoeff[0] +
                   temp * (hCoeff[1] +
                           temp * (hCoeff[2] +
                                   temp * (hCoeff[3] + temp * hCoeff[4])))) +
           hCoeff[5]);
    }
  }

#ifdef PRODRATEFILE
  ComputeThermoData(enth, cp, temp);
#elif defined PRODRATEFILEF77
  COMPTHERMODATA(enth, cp, &temp);
#endif
  return TRUE;
}

template <class TransportModel>
void T0DSpecies<TransportModel>::InitT0DSpecies(TInputDataPtr input) {
  int nOfSpecies = input->GetCounter()->species;
  int nSpeciesInSystem = nOfSpecies - input->GetCounter()->steadyStates;
  int maxGridPoints = input->fMaxGridPoints;
  int initGridPoints = input->fInitialGridPoints;

  fHeatCapacity = NewVector(nSpeciesInSystem);
  fEnthalpy = NewVector(nSpeciesInSystem);
  fViscosity = NewVector(nSpeciesInSystem);
  fConductivity = NewVector(nSpeciesInSystem);
  fProductionRate = NewVector(nSpeciesInSystem);
  fDeltaI = NewVector(nOfSpecies);
  fOmegaDOverDCoeff = NewTensor(initGridPoints + 2, nSpeciesInSystem,
                                nSpeciesInSystem, kColumnPointers);
  fOmegaDOverDCoeff->tensor = &fOmegaDOverDCoeff->tensor[kNext];
}

template <class TransportModel>
T0DSpecies<TransportModel>::~T0DSpecies(void) {
  DisposeVector(fDeltaI);
  DisposeVector(fProductionRate);
  DisposeVector(fConductivity);
  DisposeVector(fViscosity);
  DisposeVector(fEnthalpy);
  DisposeVector(fHeatCapacity);
  fOmegaDOverDCoeff->tensor = &fOmegaDOverDCoeff->tensor[kPrev];
  DisposeTensor(fOmegaDOverDCoeff);
}

template <class TransportModel>
Flag T0DSpecies<TransportModel>::ComputeSpeciesProperties(Double temp) {
  int i;
  int nSpeciesInSystem = GetNSpeciesInSystem();

  if (fabs((fCurrTemp - temp) / temp) < 1.0e-5) {
    //	if ( fCurrTemp == temp ) {\\}
    return FALSE;
  } else {
    fCurrTemp = temp;
    //		cerr << "new thermo props" << NEWL;
  }

  for (i = 0; i < nSpeciesInSystem; ++i) {
    ComputeSpeciesProperties(i, temp);
  }

#ifdef PRODRATEFILE
  ComputeThermoData(fEnthalpy->vec, fHeatCapacity->vec, temp);
#elif defined PRODRATEFILEF77
  COMPTHERMODATA(fEnthalpy->vec, fHeatCapacity->vec, &temp);
#endif
  return TRUE;
}

template <class TransportModel>
void T0DSpecies<TransportModel>::ComputeSpeciesProperties(int number,
                                                          Double temp) {
  //	computes heatCapacity, enthalpy for the
  //	species with index 'number' for the temperature temp

  Double *coeff, *hCoeff; /* pointer to nasa coefficients */
  Double r_over_m = RGAS / fMolarMass->vec[number]; /* in J / (kg K) */
  Double *heatCapacity = fHeatCapacity->vec;
  Double *viscosity = fViscosity->vec;

  if (temp > 1000.0) {
    coeff = fCoeffHigh->mat[number];
    hCoeff = fHCoeffHigh->mat[number];
  } else {
    coeff = fCoeffLow->mat[number];
    hCoeff = fHCoeffLow->mat[number];
  }
  if (temp < fTThermoLimit) {
    heatCapacity[number] =
        r_over_m *
        (coeff[0] +
         fTThermoLimit *
             (coeff[1] +
              fTThermoLimit *
                  (coeff[2] +
                   fTThermoLimit * (coeff[3] + fTThermoLimit * coeff[4]))));

    fEnthalpy->vec[number] =
        r_over_m *
        (fTThermoLimit *
             (hCoeff[0] +
              fTThermoLimit *
                  (hCoeff[1] +
                   fTThermoLimit *
                       (hCoeff[2] +
                        fTThermoLimit *
                            (hCoeff[3] + fTThermoLimit * hCoeff[4])))) +
         hCoeff[5]);
#ifdef CONSTLOWTCP
    fEnthalpy->vec[number] += (temp - fTThermoLimit) * heatCapacity[number];
#else
    Double B =
        r_over_m *
        (coeff[1] +
         fTThermoLimit * (2.0 * coeff[2] +
                          fTThermoLimit * (3.0 * coeff[3] +
                                           4.0 * fTThermoLimit * coeff[4])));
    Double cpRef = heatCapacity[number];
    Double A = cpRef - B * fTThermoLimit;
    heatCapacity[number] = A + B * temp;
    fEnthalpy->vec[number] +=
        0.5 * (heatCapacity[number] + cpRef) * (temp - fTThermoLimit);
#endif
  } else {
    heatCapacity[number] =
        r_over_m *
        (coeff[0] +
         temp * (coeff[1] +
                 temp * (coeff[2] + temp * (coeff[3] + temp * coeff[4]))));

    fEnthalpy->vec[number] =
        r_over_m *
        (temp * (hCoeff[0] +
                 temp * (hCoeff[1] +
                         temp * (hCoeff[2] +
                                 temp * (hCoeff[3] + temp * hCoeff[4])))) +
         hCoeff[5]);
  }

  viscosity[number] = TM.ComputeViscosity(number, temp);
  fConductivity->vec[number] = TM.ComputeConductivity(number, temp,
      viscosity[number], heatCapacity[number], fMolarMass->vec[number]);
}

template <class TransportModel>
template <class Flame>
void T0DSpecies<TransportModel>::PrintSpecies(Flame *flame) {
  FILE *fp = flame->GetOutfile("species", FileType::kText);

  for (int i = 0; i < GetNOfSpecies(); ++i) {
    PrintSpecies(i, flame, fp);
  }

  fclose(fp);
}

template <class TransportModel>
template <class Flame>
void T0DSpecies<TransportModel>::PrintSpecies(int number, Flame *flame,
                                              FILE *fp) {
  int i;

  fprintf(fp, "%s has number %d and appears in %d reactions, which are\n",
          fNames[number], number, fNOfUsedReactions->vec[number]);
  for (i = 0; i < fNOfUsedReactions->vec[number]; ++i) {
    flame->GetReaction()->PrintReactionEquation(fUsedReactions[number]->vec[i],
                                                this, fp);
    fprintf(fp, "\t%g\n", flame->GetReaction()
                              ->GetReactionRate()
                              ->vec[fUsedReactions[number]->vec[i]]);
    //		fprintf( fp, "%d\n", fUsedReactions[number]->vec[i] );
  }

  fprintf(fp, "HeatCapacity = %g\n", fHeatCapacity->vec[number]);
  fprintf(fp, "Enthalpy = %g\n", fEnthalpy->vec[number]);

  fprintf(fp, "fCoeffLow = ");
  for (i = 0; i < 7; ++i) {
    fprintf(fp, "\t%g", fCoeffLow->mat[number][i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "fCoeffHigh = ");
  for (i = 0; i < 7; ++i) {
    fprintf(fp, "\t%g", fCoeffHigh->mat[number][i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "molarMass = %g\n", fMolarMass->vec[number]);
  fprintf(fp, "sigma = %g\n", fSigma->vec[number]);
  TM.printTransportModelData(number, fp);
  if (!fIsSteadyState[number]) {
    fprintf(fp, "productionRate = %g\n", fProductionRate->vec[number]);
  }
  fprintf(fp, "\n\n\n");
}

template <class TransportModel>
void T1DSpecies<TransportModel>::InitT1DSpecies(TInputDataPtr input) {
  int nOfSpecies = input->GetCounter()->species;
  int nSpeciesInSystem = nOfSpecies - input->GetCounter()->steadyStates;
  int maxGridPoints = input->fMaxGridPoints;

  fViscosity = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fViscosity->mat = &fViscosity->mat[kNext];
  fHeatCapacity = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fHeatCapacity->mat = &fHeatCapacity->mat[kNext];
  fConductivity = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fConductivity->mat = &fConductivity->mat[kNext];
  fEnthalpy = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fEnthalpy->mat = &fEnthalpy->mat[kNext];
  fDiffusivity = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fDiffusivity->mat = &fDiffusivity->mat[kNext];
  fThermoDiffusivity =
      NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fThermoDiffusivity->mat = &fThermoDiffusivity->mat[kNext];
  fBinDij = NewTensor(maxGridPoints + 2, nSpeciesInSystem, nSpeciesInSystem,
                      kColumnPointers);
  fBinDij->tensor = &fBinDij->tensor[kNext];
  fGijOverWj = NewTensor(maxGridPoints + 2, nSpeciesInSystem, nSpeciesInSystem,
                         kColumnPointers);
  fGijOverWj->tensor = &fGijOverWj->tensor[kNext];
  fOmegaDOverDCoeff = NewTensor(maxGridPoints + 2, nSpeciesInSystem,
                                nSpeciesInSystem, kColumnPointers);
  fOmegaDOverDCoeff->tensor = &fOmegaDOverDCoeff->tensor[kNext];

  if (input->fThermoDiffusion) {
    fDThermConst = NewTensor(maxGridPoints + 2, nSpeciesInSystem,
                             nSpeciesInSystem, kColumnPointers);
    fDThermConst->tensor = &fDThermConst->tensor[kNext];
  } else {
    fDThermConst = NULL;
  }

  fProductionRate = NewMatrix(nSpeciesInSystem, maxGridPoints, kColumnPointers);

  fSavedDeltaiY = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fSavedDeltaiY->mat = &fSavedDeltaiY->mat[kNext];
  fSavedY = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fSavedY->mat = &fSavedY->mat[kNext];
  fSumDiff = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fSumDiff->mat = &fSumDiff->mat[kNext];
  fDeltaI = NewMatrix(nOfSpecies, maxGridPoints + 2, kColumnPointers);
  fDeltaI->mat = &fDeltaI->mat[kNext];
  fTempProp = NewVector(maxGridPoints + 2);
  fTempProp->vec = &fTempProp->vec[kNext];
  fPressureProp = NewVector(maxGridPoints + 2);
  fPressureProp->vec = &fPressureProp->vec[kNext];
}

template <class TransportModel>
T1DSpecies<TransportModel>::~T1DSpecies(void) {
  fTempProp->vec = &fTempProp->vec[kPrev];
  DisposeVector(fTempProp);
  fPressureProp->vec = &fPressureProp->vec[kPrev];
  DisposeVector(fPressureProp);
  fDeltaI->mat = &fDeltaI->mat[kPrev];
  DisposeMatrix(fDeltaI);
  fSumDiff->mat = &fSumDiff->mat[kPrev];
  DisposeMatrix(fSumDiff);
  fSavedDeltaiY->mat = &fSavedDeltaiY->mat[kPrev];
  DisposeMatrix(fSavedDeltaiY);
  fSavedY->mat = &fSavedY->mat[kPrev];
  DisposeMatrix(fSavedY);

  DisposeMatrix(fProductionRate);

  if (fDThermConst) {
    fDThermConst->tensor = &fDThermConst->tensor[kPrev];
    DisposeTensor(fDThermConst);
  }
  fOmegaDOverDCoeff->tensor = &fOmegaDOverDCoeff->tensor[kPrev];
  DisposeTensor(fOmegaDOverDCoeff);
  fGijOverWj->tensor = &fGijOverWj->tensor[kPrev];
  DisposeTensor(fGijOverWj);
  fBinDij->tensor = &fBinDij->tensor[kPrev];
  DisposeTensor(fBinDij);
  fThermoDiffusivity->mat = &fThermoDiffusivity->mat[kPrev];
  DisposeMatrix(fThermoDiffusivity);
  fDiffusivity->mat = &fDiffusivity->mat[kPrev];
  DisposeMatrix(fDiffusivity);
  fEnthalpy->mat = &fEnthalpy->mat[kPrev];
  DisposeMatrix(fEnthalpy);
  fConductivity->mat = &fConductivity->mat[kPrev];
  DisposeMatrix(fConductivity);
  fHeatCapacity->mat = &fHeatCapacity->mat[kPrev];
  DisposeMatrix(fHeatCapacity);
  fViscosity->mat = &fViscosity->mat[kPrev];
  DisposeMatrix(fViscosity);
}

template <class TransportModel>
void T1DSpecies<TransportModel>::PrintSpecies(int k) {
  FILE *fp;

  if (!(fp = fopen("species.tout", "w"))) {
    cerr << "#warning: unable to open file 'species.tout'" << NEWL;
    return;
  }
  fprintf(fp, "GridPoint no. %d\n", k);
  for (int i = 0; i < GetNOfSpecies(); ++i) {
    PrintSpecies(i, k, fp);
  }
  fclose(fp);
}

template <class TransportModel>
void T1DSpecies<TransportModel>::PrintProductionRate(TNewtonPtr bt) {
  // print productionrate in [mole/cm^3]

  int i, k;
  char fileName[32];
  FILE *fp = NULL;
  Double **mat = fProductionRate->mat;
  Double *molarMass = fMolarMass->vec;
  char **names = GetNames();
  TGridPtr currentGrid = bt->GetGrid()->GetCurrentGrid();
  Double *x = currentGrid->GetX()->vec;
  int gridPoints = currentGrid->GetNGridPoints();
  static int counter = 0;

  if (gridPoints > fProductionRate->cols) {
    fprintf(stderr,
            "#warning: productionrate not printed, because it doesn't "
            "correspond to current grid\n");
    return;
  }

  sprintf(fileName, "ProdRate%d.dout", ++counter);
  if (!(fp = fopen(fileName, "w"))) {
    cerr << "#warning: unable to open file " << fileName << NEWL;
    return;
  }
  if (counter >= 10) counter = 0;

  fprintf(fp, "*\n");
  fprintf(fp, "%-9s\t%-12s", "gridPoint", "eta");
  for (i = 0; i < GetNSpeciesInSystem(); ++i) {
    fprintf(fp, "\t%-12s", names[i]);
  }

  fprintf(fp, "\n");

  for (k = 0; k < gridPoints; ++k) {
    fprintf(fp, "%9d\t%12.4E", k, x[k]);
    for (i = 0; i < GetNSpeciesInSystem(); ++i) {
      fprintf(fp, "\t%12.4E", mat[k][i] / molarMass[i] * 1.0e-3);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

template <class TransportModel>
void T1DSpecies<TransportModel>::PrintSpecies(int number, int gridPoint,
                                              FILE *fp) {
  int i;
  static Flag init = FALSE;

  if (!fp) {
    if (!init) {
      fp = fopen("species.tout", "w");
      fprintf(fp, "GridPoint no. %d\n", gridPoint);
    } else {
      fp = fopen("species.tout", "a");
    }
  }
  if (!fp) {
    cerr << "#warning: unable to open file 'species.tout'" << NEWL;
    return;
  }
  init = TRUE;
  if (!fIsSteadyState[number]) {
    fprintf(fp, "%s has number %d and appears in %d reactions, which are\n",
            fNames[number], number, fNOfUsedReactions->vec[number]);
    for (i = 0; i < fNOfUsedReactions->vec[number]; ++i) {
      fprintf(fp, "%d\n", fUsedReactions[number]->vec[i]);
    }
  } else {
    fprintf(fp, "%s has number %d and is steady state\n", fNames[number],
            number);
  }
  fprintf(fp, "Viscosity = %g\n", fViscosity->mat[gridPoint][number]);
  fprintf(fp, "HeatCapacity = %g\n", fHeatCapacity->mat[gridPoint][number]);
  fprintf(fp, "Conductivity = %g\n", fConductivity->mat[gridPoint][number]);
  fprintf(fp, "Enthalpy = %g\n", fEnthalpy->mat[gridPoint][number]);
  fprintf(fp, "Diffusivity = %g\n", fDiffusivity->mat[gridPoint][number]);
  if (IsConstantLewisNumber()) {
    fprintf(fp, "LewisNumber = %g\n", fLewisNumber->vec[number]);
  }

  fprintf(fp, "fCoeffLow = ");
  for (i = 0; i < 7; ++i) {
    fprintf(fp, "\t%g", fCoeffLow->mat[number][i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "fCoeffHigh = ");
  for (i = 0; i < 7; ++i) {
    fprintf(fp, "\t%g", fCoeffHigh->mat[number][i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "molarMass = %g\n", fMolarMass->vec[number]);
  fprintf(fp, "sigma = %g\n", fSigma->vec[number]);
  TM.printTransportModelData(number, fp);

  if (!fIsSteadyState[number]) {
    fprintf(fp, "productionRate = %g\n",
            fProductionRate->mat[gridPoint][number]);
  }
  fprintf(fp, "\n\n\n");
}

template <class TransportModel>
Flag T1DSpecies<TransportModel>::ComputeSpeciesProperties(
    TFlameNodePtr flameNode, Double temp, Double pressure) {
  int i;
  int nSpeciesInSystem = GetNSpeciesInSystem();

#ifdef OPTPROPS
  if (flameNode->tempProp[kCurr] == temp) {
    return FALSE;
  } else {
    flameNode->tempProp[kCurr] = temp;
  }
#else
  flameNode->tempProp[kCurr] = temp;
#endif  // OPTPROPS

#define INLINESPECIESPROPS
#ifdef INLINESPECIESPROPS
  Double **coMat = fCoeffHigh->mat,
         **hCoMat = fCoeffLow->mat; /* pointer to nasa coefficients */
  Double *coeff, *hCoeff;           /* pointer to nasa coefficients */
  Double *heatCapacity = flameNode->heatCapacity;
  Double *vis = flameNode->viscosity;
  Double *ent = flameNode->enthalpy;
  Double *conduc = flameNode->conductivity;

  Double *W = fMolarMass->vec;
  Double sqTemp = sqrt(temp);
  Double R_gas = RGAS;
  Double r_over_m;

  if (temp < 10.0) {
    temp = 10.0;
  }

  if (temp > 1000.0) {
    coMat = fCoeffHigh->mat;
    hCoMat = fHCoeffHigh->mat;
  } else {
    coMat = fCoeffLow->mat;
    hCoMat = fHCoeffLow->mat;
  }

  /*	Double	t2 = temp * temp;*/
  /*	Double	t3 = t2 * temp;*/
  /*	Double	t4 = t3 * temp;*/
  /*	Double	t5 = t4 * temp;*/
  for (i = 0; i < nSpeciesInSystem; ++i) {
    coeff = coMat[i];
    hCoeff = hCoMat[i];
    r_over_m = R_gas / W[i];
    if (temp < fTThermoLimit) {
      heatCapacity[i] =
          r_over_m *
          (coeff[0] +
           fTThermoLimit *
               (coeff[1] +
                fTThermoLimit *
                    (coeff[2] +
                     fTThermoLimit * (coeff[3] + fTThermoLimit * coeff[4]))));

      ent[i] = r_over_m *
               (fTThermoLimit *
                    (hCoeff[0] +
                     fTThermoLimit *
                         (hCoeff[1] +
                          fTThermoLimit *
                              (hCoeff[2] +
                               fTThermoLimit *
                                   (hCoeff[3] + fTThermoLimit * hCoeff[4])))) +
                hCoeff[5]);
#ifdef CONSTLOWTCP
      ent[i] += (temp - fTThermoLimit) * heatCapacity[i];
#else
      Double B =
          r_over_m *
          (coeff[1] +
           fTThermoLimit * (2.0 * coeff[2] +
                            fTThermoLimit * (3.0 * coeff[3] +
                                             4.0 * fTThermoLimit * coeff[4])));
      Double cpRef = heatCapacity[i];
      Double A = cpRef - B * fTThermoLimit;
      heatCapacity[i] = A + B * temp;
      ent[i] += 0.5 * (heatCapacity[i] + cpRef) * (temp - fTThermoLimit);
#endif  // CONSTLOWTCP
    } else {
      heatCapacity[i] =
          r_over_m *
          (coeff[0] +
           temp * (coeff[1] +
                   temp * (coeff[2] + temp * (coeff[3] + temp * coeff[4]))));

      ent[i] =
          r_over_m *
          (temp * (hCoeff[0] +
                   temp * (hCoeff[1] +
                           temp * (hCoeff[2] +
                                   temp * (hCoeff[3] + temp * hCoeff[4])))) +
           hCoeff[5]);
    }
    /*		heatCapacity[i] = (coeff[0] + coeff[1]*temp + coeff[2]*t2 +
     * coeff[3]*t3 + coeff[4]*t4 ) * r_over_m;*/
    /*		ent[i] = (hCoeff[5] + hCoeff[0]*temp + hCoeff[1]*t2 +
     * hCoeff[2]*t3
     * +
     * hCoeff[3]*t4 + hCoeff[4]*t4 ) * r_over_m;*/

    vis[i] = TM.ComputeViscosity(i, temp);
    conduc[i] = TM.ComputeConductivity(i, temp, vis[i], heatCapacity[i], W[i]);
  }
#else   // NO INLINESPECIESPROPS
  for (i = 0; i < nSpeciesInSystem; ++i) {
    ComputeSpeciesProperties(flameNode, temp, pressure, i);
  }
#endif  // INLINESPECIESPROPS
#if defined(OPTIMIZEDELTAI) || defined(OPTOMEGAD) || defined(OPTDTHERM)
  int j;
  int nSpeciesIn = GetNSpeciesInSystem();
  Double **GijOverWj = flameNode->GijOverWj;
  Double *viscosity = flameNode->viscosity;
  Double *W_i = fMolarMass->vec;
  const Double pressureInvT3 = pressure / (temp * sqrt(temp));
  flameNode->pressureProp[kCurr] = pressure;
#ifdef OPTIMIZEDELTAI
  Double *GijOverWj_i, *sqrtMiOverMjPlusOne_i;
  for (i = 0; i < nSpeciesInSystem; ++i) {
    GijOverWj_i = GijOverWj[i];
    for (j = 0; j <= i; ++j) {
      GijOverWj_i[j] = sqrt(sqrtMjOverMi[j][i] * viscosity[i] / viscosity[j]);
      //	(old form before bug was fixed)	GijOverWj_i[j] = sqrt(
      // sqrtMjOverMi[j][i] * viscosity[j] / viscosity[i] );
    }
  }
  for (i = 0; i < nSpeciesInSystem; ++i) {
    GijOverWj_i = GijOverWj[i];
    for (j = i + 1; j < nSpeciesInSystem; ++j) {
      GijOverWj_i[j] = 1.0 / GijOverWj[j][i];
    }
  }
  for (i = 0; i < nSpeciesInSystem; ++i) {
    GijOverWj_i = GijOverWj[i];
    sqrtMiOverMjPlusOne_i = sqrtMiOverMjPlusOne[i];
    for (j = 0; j < nSpeciesInSystem; ++j) {
      GijOverWj_i[j] += 1.0;
      GijOverWj_i[j] *= GijOverWj_i[j] * sqrtMiOverMjPlusOne_i[j] / W_i[j];
      /*			GijOverWj[i][j] = sqrt( sqrtMjOverMi[j][i] *
       * viscosity[j] / viscosity[i] ) + 1.0;*/
      /*			GijOverWj[i][j] *= GijOverWj[i][j] *
       * sqrtMiOverMjPlusOne[i][j] */
      /*								/
       * W_i[j];*/
    }
  }
#endif  // OPTIMIZEDELTAI
#ifdef OPTOMEGAD
  //### old
  /*	GetFuncOne( flameNode->OmegaDOverDCoeff, nSpeciesInSystem,
   * pressureInvT3, omegaCoeff,*/
  /*				dCoeff );*/
  // new
  /*	GetFuncOne( flameNode->OmegaDOverDCoeff, nSpeciesInSystem,
   * pressureInvT3, k_over_epsjk, redDipoleMomentjk, */
  /*				dCoeff );*/
  Double **omegaDOverDCoeff = flameNode->OmegaDOverDCoeff;
  for (i = 0; i < nSpeciesInSystem; ++i) {
    for (j = 0; j < i; ++j)
      omegaDOverDCoeff[i][j] =
          TM.InvBinaryDiffusionCoefficient(i, j, pressureInvT3, temp);
    omegaDOverDCoeff[i][i] = 0.0;
  }
  for (i = 0; i < nSpeciesInSystem; ++i) {
    for (j = i + 1; j < nSpeciesInSystem; ++j) {
      omegaDOverDCoeff[i][j] = omegaDOverDCoeff[j][i];
    }
  }
#endif  // OPTOMEGAD
#endif  // defined (OPTIMIZEDELTAI) || defined (OPTOMEGAD) || defined
        // (OPTDTHERM)
  return TRUE;
}

/*void TSpecies::GetFuncOne( Double **omegaDOverDCoeff, int nSpeciesInSystem,
 * Double pressureInvT3, double **k_over_eps, double redDipoleMoment, */
/*				Double **dCoeff )*/
/*{*/
/*	int i, j;*/
/*	for( i = 0; i < nSpeciesInSystem; ++i ) {*/
/*		for ( j = 0; j < i; ++j ) {*/
/*	// ### old */
/*	//		omegaDOverDCoeff[i][j] = pressureInvT3 * omega_D( temp *
 * k_over_eps[i][j] ) / dCoeff[i][j];*/
/*	// ### new */
/*			omegaDOverDCoeff[i][j] = pressureInvT3 * omega_D( temp *
 * k_over_eps[i][j], redDipoleMoment[i][j] ) / dCoeff[i][j];*/
/*		}*/
/*		omegaDOverDCoeff[i][i] = 0.0;*/
/*	}*/
/*	for( i = 0; i < nSpeciesInSystem; ++i ) {*/
/*		for ( j = i+1; j < nSpeciesInSystem; ++j ) {*/
/*			omegaDOverDCoeff[i][j] = omegaDOverDCoeff[j][i];*/
/*		}*/
/*	}*/
/*}*/

template <class TransportModel>
void T1DSpecies<TransportModel>::ComputeSpeciesProperties(
    TFlameNodePtr flameNode, Double temp, Double /*pressure*/, int number) {
  // computes heatCapacity, viscosity, enthalpy and thermal conductivity for the
  // species with index 'number' for the temperature temp

  Double *coeff, *hCoeff; /* pointer to nasa coefficients */
  Double R_gas = RGAS;
  Double r_over_m = R_gas / fMolarMass->vec[number]; /* in J / (kg K) */
  Double *heatCapacity = flameNode->heatCapacity;
  Double *viscosity = flameNode->viscosity;

  if (temp < 10.0) {
    temp = 10.0;
  }
  if (temp > 1000.0) {
    coeff = fCoeffHigh->mat[number];
    hCoeff = fHCoeffHigh->mat[number];
  } else {
    coeff = fCoeffLow->mat[number];
    hCoeff = fHCoeffLow->mat[number];
  }

  if (temp < fTThermoLimit) {
    heatCapacity[number] =
        r_over_m *
        (coeff[0] +
         fTThermoLimit *
             (coeff[1] +
              fTThermoLimit *
                  (coeff[2] +
                   fTThermoLimit * (coeff[3] + fTThermoLimit * coeff[4]))));

    flameNode->enthalpy[number] =
        r_over_m *
        (fTThermoLimit *
             (hCoeff[0] +
              fTThermoLimit *
                  (hCoeff[1] +
                   fTThermoLimit *
                       (hCoeff[2] +
                        fTThermoLimit *
                            (hCoeff[3] + fTThermoLimit * hCoeff[4])))) +
         hCoeff[5]);
#ifdef CONSTLOWTCP
    flameNode->enthalpy[number] +=
        (temp - fTThermoLimit) * heatCapacity[number];
#else
    Double B =
        r_over_m *
        (coeff[1] +
         fTThermoLimit * (2.0 * coeff[2] +
                          fTThermoLimit * (3.0 * coeff[3] +
                                           4.0 * fTThermoLimit * coeff[4])));
    Double A = heatCapacity[number] - B * fTThermoLimit;
    Double cpRef = heatCapacity[number];
    heatCapacity[number] = A + B * temp;
    flameNode->enthalpy[number] +=
        0.5 * (heatCapacity[number] + cpRef) * (temp - fTThermoLimit);
#endif
  } else {
    heatCapacity[number] =
        r_over_m *
        (coeff[0] +
         temp * (coeff[1] +
                 temp * (coeff[2] + temp * (coeff[3] + temp * coeff[4]))));

    flameNode->enthalpy[number] =
        r_over_m *
        (temp * (hCoeff[0] +
                 temp * (hCoeff[1] +
                         temp * (hCoeff[2] +
                                 temp * (hCoeff[3] + temp * hCoeff[4])))) +
         hCoeff[5]);
  }
  viscosity[number] = TM.ComputeViscosity(number, temp);
  flameNode->conductivity[number] = TM.ComputeConductivity(number, temp,
      viscosity[number], heatCapacity[number], fMolarMass->vec[number]);
}

template <class TransportModel>
void T1DSpecies<TransportModel>::Compute_D(TFlameNodePtr flameNode) {
  int i;
  int nSpeciesInSystem = GetNSpeciesInSystem();
  Double *diffusivity = flameNode->diffusivity;
  Double num = *flameNode->mixConductivity /
               (*flameNode->mixDensity * *flameNode->mixHeatCapacity);
  Double *lewis = fLewisNumber->vec;

  for (i = 0; i < nSpeciesInSystem; ++i) {
    diffusivity[i] = num / lewis[i];
  }
}

#ifdef FULLDIFFUSION
template <class TransportModel>
void T1DSpecies<TransportModel>::Compute_D(TFlameNodePtr flameNode, Double temp,
                                           Double *Y, Double press,
                                           Flag newTemp) {
  int i, j, k;
  int nSpeciesIn = GetNSpeciesInSystem();
  Double **mat = NULL;
  Double **inv = NULL;
  Double *col = NULL;
  int *index = NULL;

  Double **invDij = fInvDij->mat;
  Double **binDij = flameNode->binDij[kCurr];
  Double *mw = fMolarMass->vec;
  Double pressureInvT3 = press / (temp * sqrt(temp));
  const Double zero = 0.0;

  mat = New2DArray(nSpeciesIn, nSpeciesIn);
  inv = New2DArray(nSpeciesIn, nSpeciesIn);
  index = New1DIntArray(nSpeciesIn);
  col = New1DArray(nSpeciesIn);

  //	compute the reciprocal binary diffusion coefficients
  for (i = 0; i < nSpeciesIn; ++i) {
    invDij[i][i] = zero;
    for (j = 0; j < i; ++j) {
      invDij[i][j] = invDij[j][i] =
#ifdef OPTOMEGAD
          flameNode
              ->OmegaDOverDCoeff[i][j];  // pressureInvT3 is in OmegaDOverDCoeff
#else
          TM.InvBinaryDiffusionCoefficient(i, j, pressureInvT3, temp);
#endif  // OPTOMEGAD
    }
  }

  //	assemble the matrix
  for (i = 0; i < nSpeciesIn; ++i) {
    for (j = 0; j < nSpeciesIn; ++j) {
      if (i == j) continue;

      for (k = 0; k < nSpeciesIn; ++k) {
        if (k == i) continue;

        mat[i][j] += Y[k] / mw[k] * invDij[i][k];
      }
      mat[i][j] *= mw[j];
      mat[i][j] += Y[i] * invDij[i][j];
    }
  }

  InverseMatrix(nSpeciesIn, mat, inv, index, col);

  for (i = 0; i < nSpeciesIn; ++i) {
    for (j = 0; j < nSpeciesIn; ++j) {
      binDij[i][j] = inv[i][j] - mw[i] / mw[j] * inv[i][i];
    }
  }

  Free1DArray(col);
  Free1DIntArray(index);
  Free2DArray(inv);
  Free2DArray(mat);
}
#else
template <class TransportModel>
void T1DSpecies<TransportModel>::Compute_D(TFlameNodePtr flameNode, Double temp,
                                           Double *Y, Double pressure,
                                           Flag newTemp) {
  int j, i;
  int nSpeciesInSystem = GetNSpeciesInSystem();
  //	static int		cont1 = 0, cont2 = 0, contAll = 0, cont0 = 0;
  Double *molarMass = fMolarMass->vec;
  Double *diffusivity = flameNode->diffusivity;
  Double mixtureMolarMass = *flameNode->mixMolarMass;
  Double **invDij = fInvDij->mat;
  Double sumYi = 0.0;
  Double pressureInvT3 = pressure / (temp * sqrt(temp));
  Double *savedY = flameNode->savedY;
  Double *sumDiff = flameNode->sumDiff;
  int changed =
      -1;  // -1 means nothing changed, -2 means temperature or pressure
           // or more than one Y_i changed, >=0 is changed species, if only
           // one changed
  int changed2 =
      -1;  // -1 means nothing changed, -2 means temperature or pressure

  if (newTemp) {
    changed = -2;
  } else {
    for (i = 0; i < nSpeciesInSystem; ++i) {
      if (Y[i] != savedY[i]) {
        if (changed2 < 0) {
          if (changed == -1) {
            changed = i;
          } else {
            changed2 = i;
          }
        } else {
          changed = -2;
          break;
        }
      }
    }
  }

#ifdef OPTOMEGAD
  if (flameNode->pressureProp[kCurr] != pressure) {
    changed = -2;
    for (i = 0; i < nSpeciesInSystem; ++i) {
      invDij[i][i] = 0.0;
      for (j = 0; j < i; ++j) {
        invDij[i][j] = invDij[j][i] =
            flameNode->OmegaDOverDCoeff[i][j] * pressure /
            flameNode
                ->pressureProp[kCurr];  // pressureInvT3 is in OmegaDOverDCoeff
      }
      flameNode->pressureProp[kCurr] = pressure;
    }
  } else {
    //			invDij[i][j] = invDij[j][i]  =
    // flameNode->OmegaDOverDCoeff[i][j];
    invDij = flameNode->OmegaDOverDCoeff;
  }
#else   // NO OPTOMEGAD
  for (i = 0; i < nSpeciesInSystem; ++i) {
    for (j = 0; j < i; ++j)
      invDij[i][j] = invDij[j][i] =
          TM.InvBinaryDiffusionCoefficient(i, j, pressureInvT3, temp);
  }
#endif  // OPTOMEGAD

  if (changed >= 0) {
    for (i = 0; i < nSpeciesInSystem; ++i) {
      sumYi += Y[i];
    }
    for (i = 0; i < nSpeciesInSystem; ++i) {
      sumDiff[i] = sumDiff[i] +
                   invDij[i][changed] / molarMass[changed] *
                       (Y[changed] - savedY[changed]);
      if (changed2 >= 0) {
        sumDiff[i] = sumDiff[i] +
                     invDij[i][changed2] / molarMass[changed2] *
                         (Y[changed2] - savedY[changed2]);
      }

      if (sumDiff[i] > 1.0e-15) {
        diffusivity[i] = (sumYi - Y[i]) / (mixtureMolarMass * sumDiff[i]);
      } else {
        for (j = 1, sumDiff[i] = 0.0; j < nSpeciesInSystem; ++j) {
          sumDiff[i] += invDij[i][j] / molarMass[j];
        }
        diffusivity[i] =
            ((Double)nSpeciesInSystem - 1.0) / (mixtureMolarMass * sumDiff[i]);
      }
      //			if ( isnan( diffusivity[i] ) ) {
      //				fprintf( stderr, "diff[%s] is %g\n",
      // fNames[i], diffusivity[i] );
      //				diffusivity[i] = 0.0;
      //			}
    }
    savedY[changed] = Y[changed];
    if (changed2 >= 0) {
      savedY[changed2] = Y[changed2];
      //			++cont2;
    }
    //		else {
    //			++cont1;
    //		}
  } else if (changed == -2) {
    //		++contAll;
    Double *YiOverWi = New1DArray(nSpeciesInSystem);
    sumYi = 0.0;
    for (i = 0; i < nSpeciesInSystem; ++i) {
      YiOverWi[i] = Y[i] / molarMass[i];
      sumYi += Y[i];
    }
    for (i = 0; i < nSpeciesInSystem; ++i) {
      savedY[i] = Y[i];

      sumDiff[i] = YiOverWi[0] * invDij[i][0];
      for (j = 1; j < nSpeciesInSystem; ++j) {
        sumDiff[i] += YiOverWi[j] * invDij[i][j];
      }

      if (sumDiff[i] > 1.0e-15) {
        diffusivity[i] = (sumYi - Y[i]) / (mixtureMolarMass * sumDiff[i]);
      } else {
        for (j = 1, sumDiff[i] = invDij[i][0] / molarMass[0];
             j < nSpeciesInSystem; ++j) {
          sumDiff[i] += invDij[i][j] / molarMass[j];
        }
        diffusivity[i] =
            ((Double)nSpeciesInSystem - 1.0) / (mixtureMolarMass * sumDiff[i]);
      }
      //			if ( isnan( diffusivity[i] ) ) {
      //				fprintf( stderr, "diff[%s] is %g\n",
      // fNames[i], diffusivity[i] );
      //				fprintf( stderr, "sumDiff[%s] is %g\n",
      // fNames[i], sumDiff[i] );
      //				fprintf( stderr, "sumYi is %g\n", sumYi
      //);
      //				fprintf( stderr, "Y[%s] is %g\n",
      // fNames[i], Y[i] );
      //				fprintf( stderr, "mixtureMolarMass is
      //%g\n", mixtureMolarMass );
      //				fprintf( stderr, "temp is %g\n", temp );
      //				for ( j = 0; j < nSpeciesInSystem; ++j )
      //{
      //					fprintf( stderr, "YiOverWi[%s]
      // is %g\n", fNames[j], YiOverWi[j] );
      //					fprintf( stderr, "invDij[%s][%s]
      // is %g\n", fNames[i], fNames[j], invDij[i][j] );
      //				}
      //				exit(2);
      //			}
      /*		if ( diffusivity[i] < 0.0 ) {*/
      /*			fprintf( stderr, "diff[%s] is %g\n", fNames[i],
       * diffusivity[i] );*/
      /*		}*/
    }
    Free1DArray(YiOverWi);
  } else {
    // do nothing
    //		++cont0;
  }

  //	fprintf( stderr, "cont0 = %d\tcont2 = %d\tcont1 = %d\tcontAll = %d\n",
  // cont0, cont1, cont2, contAll );
}
#endif

template <class TransportModel>
void T0DSpecies<TransportModel>::ComputeDiffusivityTrans1D(
    Double *diffusivityTrans1D_k, int GPCurr, Double temp, Double *Y,
    Double mixMolMass, int NOSpecies, Double *molarMass, Double pressure) {
  int j, i;
  Double *YiOverWi = New1DArray(NOSpecies);
  Double *sumDiff = New1DArray(NOSpecies);
  Double sumYi = 0.0;
  Double pressureInvT3 = pressure / (temp * sqrt(temp));
  Double **invDij = GetOmegaDOverDCoeff()->tensor[GPCurr];

  for (i = 0; i < NOSpecies; ++i) {
    for (j = 0; j < i; ++j)
      invDij[i][j] =
          TM.InvBinaryDiffusionCoefficient(i, j, pressureInvT3, temp);
    invDij[i][i] = 0.0;
  }
  for (i = 0; i < NOSpecies; ++i) {
    for (j = i + 1; j < NOSpecies; ++j) invDij[i][j] = invDij[j][i];
  }

  for (i = 0; i < NOSpecies; ++i) {
    YiOverWi[i] = Y[i] / molarMass[i];
    sumYi += Y[i];
  }
  for (i = 0; i < NOSpecies; ++i) {
    sumDiff[i] = YiOverWi[0] * invDij[i][0];
    for (j = 1; j < NOSpecies; ++j) {
      sumDiff[i] += YiOverWi[j] * invDij[i][j];
    }
    if (sumDiff[i] > 1.0e-15) {
      diffusivityTrans1D_k[i] = (sumYi - Y[i]) / (mixMolMass * sumDiff[i]);
    } else {
      for (j = 1, sumDiff[i] = invDij[i][0] / molarMass[0]; j < NOSpecies;
           ++j) {
        sumDiff[i] += invDij[i][j] / molarMass[j];
      }
      diffusivityTrans1D_k[i] =
          ((Double)NOSpecies - 1.0) / (mixMolMass * sumDiff[i]);
    }
  }
  Free1DArray(sumDiff);
  Free1DArray(YiOverWi);
}

#ifdef OPTDTHERM
template <class TransportModel>
void T1DSpecies<TransportModel>::Compute_DTherm(TFlameNodePtr flameNode,
                                                Flag calcNewConst) {
  int i, j;
  int nSpecIn = GetNSpeciesInSystem();
  Double *M = fMolarMass->vec;
  Double *lambda = flameNode->conductivity;
  Double *Y = flameNode->Y[kCurr];
  Double *diffTherm = flameNode->diffTherm[kCurr];
  Double *deltaI = flameNode->deltaI;
  Double diff;
  Double lambdaOverDeltaI;
  const Double fact = (6.0 / 5.0 * 0.93 - 1.0) / RGAS;
  Double **DThermConst = flameNode->DThermConst;

  if (calcNewConst) {
    for (i = 0; i < nSpecIn; ++i) {
      lambdaOverDeltaI = lambda[i] / deltaI[i];
      for (j = 0, diff = 0.0; j < nSpecIn; ++j) {
        DThermConst[i][j] =
            M[j] / (M[i] + M[j]) * (lambda[j] / deltaI[j] - lambdaOverDeltaI);
        diff += Y[j] * DThermConst[i][j];
      }
      diffTherm[i] = diff * fact * Y[i] * M[i];
    }
  } else {
    for (i = 0; i < nSpecIn; ++i) {
      for (j = 0, diff = 0.0; j < nSpecIn; ++j) {
        diff += Y[j] * DThermConst[i][j];
      }
      diffTherm[i] = diff * fact * Y[i] * M[i];
    }
  }
}

#else

template <class TransportModel>
void T1DSpecies<TransportModel>::Compute_DTherm(TFlameNodePtr flameNode,
                                                Flag /*calcNewConst*/) {
  int i, j;
  int nSpecIn = GetNSpeciesInSystem();
  Double *M = fMolarMass->vec;
  Double *lambda = flameNode->conductivity;
  Double *Y = flameNode->Y[kCurr];
  Double *diffTherm = flameNode->diffTherm[kCurr];
  Double *deltaI = flameNode->deltaI;
  Double diff;
  Double lambdaOverDeltaI;
  const Double fact = (6.0 / 5.0 * 0.93 - 1.0) / RGAS;

  for (i = 0; i < nSpecIn; ++i) {
    lambdaOverDeltaI = lambda[i] / deltaI[i];
    for (j = 0, diff = 0.0; j < nSpecIn; ++j) {
      diff += Y[j] * M[j] / (M[i] + M[j]) *
              (lambda[j] / deltaI[j] - lambdaOverDeltaI);
    }
    diffTherm[i] = diff * fact * Y[i] * M[i];
  }
}
#endif

template <class TransportModel>
Double T1DSpecies<TransportModel>::ComputeOneDiffusionCoeff(
    int number, Double *Y, Double *invDij, Double mixMolarMass) {
  int nSpeciesInSystem = GetNSpeciesInSystem();
  Double *molarMass = fMolarMass->vec;
  Double sum = Y[0] * invDij[0] / molarMass[0]; /* j = 0 */

  for (int j = 1; j < nSpeciesInSystem; ++j)
    sum += Y[j] * invDij[j] / molarMass[j];
  if (sum) {
    return (1.0 - Y[number]) / (sum * mixMolarMass);
  } else {
    cerr << "warning: TSpecies::ComputeOneDiffusionCoeff cannot handle this "
            "case ( Y = 1.0 )"
         << NEWL;
    return 0.0;
  }
}

template <class TransportModel>
template <class Flame>
void T1DSpecies<TransportModel>::PrintProdRateTerms(const char *name,
                                                    Flame* flame) {
  int i;

  if ((i = FindSpecies(name)) >= 0) {
    PrintProdRateTerms(i, flame);
  } else {
    fprintf(stderr, "###warning: no species %s\n", name);
    return;
  }
}

template <class TransportModel>
template <class Flame>
void T1DSpecies<TransportModel>::PrintProdRateTerms(int i,
                                                    Flame* flame) {
  TNewtonPtr bt = flame->GetSolver()->bt;
  int k;
  FILE *fp;
  TGridPtr grid = bt->GetGrid()->GetCurrentGrid();
  Double *x = grid->GetX()->vec;
  int nGridPoints = bt->GetCurrentGridPoints();
  Double dummy;
  int counter;

  int j;
  T1DReaction *reaction = flame->GetReaction();
  Double *nu = fNu[i]->vec;
  int *usedReactions = fUsedReactions[i]->vec;
  int nUsedReactions = fNOfUsedReactions->vec[i];
  char **labels = reaction->GetLabels();
  Double **reactionRate = reaction->GetReactionRate()->mat;
  Double Mi = fMolarMass->vec[i];
  Double prodRate;
  char reac[128];

  counter = (int)(modf(bt->GetNIter() / 10.0, &dummy) * 10.0);
  sprintf(flame->GetOutFileBuff(), "%s%s_ProdReac_%d.dout",
          flame->GetOutputPath(), fNames[i], counter);
  if (!(fp = fopen(flame->GetOutFileBuff(), "w"))) {
    cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
    return;
  }

  //	header
  fprintf(fp, "*\n");
  fprintf(fp, "%-12s", "y");
  for (j = 0; j < nUsedReactions; ++j) {
    //		fprintf( fp, "\t%-12s", labels[usedReactions[j]] );
    //		fprintf( fp, "\t%s: ", labels[usedReactions[j]] );
    reaction->PrintReactionEquation(usedReactions[j], flame->GetSpecies(),
                                    reac);
    //		fprintf( stderr, "%s\n", reac );
    fprintf(fp, "\t%s: %-s", labels[usedReactions[j]], reac);
  }

  for (k = 0; k < nGridPoints; ++k) {
    fprintf(fp, "\n%-12E", x[k]);
    for (j = 0; j < nUsedReactions; ++j) {
      prodRate = -Mi * nu[j] * reactionRate[k][usedReactions[j]];
      fprintf(fp, "\t%-12E", prodRate);
    }
  }

  fclose(fp);
}

template <class TransportModel>
template <class Flame>
void T1DSpecies<TransportModel>::PrintProdCons(TNewtonPtr bt,
                                               Flame* flame) {
  int i, k;
  int nGridPoints = bt->GetCurrentGridPoints();
  int nSpecies = GetNSpeciesInSystem();
  Double **reactionRate = flame->GetReaction()->GetReactionRate()->mat;
  FILE *fp;
  TGridPtr grid = bt->GetGrid()->GetCurrentGrid();
  Double *x = grid->GetX()->vec;
  Double **source = New2DArray(nGridPoints, nSpecies);
  Double **sink = New2DArray(nGridPoints, nSpecies);

  for (k = 0; k < nGridPoints; ++k) {
    CompProdCons(source[k], sink[k], reactionRate[k]);
  }

  sprintf(flame->GetOutFileBuff(), "%sProductionRates.dout",
          flame->GetOutputPath());
  if (!(fp = fopen(flame->GetOutFileBuff(), "w"))) {
    cerr << "#warning: unable to open file " << flame->GetOutFileBuff() << NEWL;
    return;
  }
  fprintf(fp, "*\n");
  fprintf(fp, "%-12s", "eta");
  for (i = 0; i < nSpecies; ++i) {
    fprintf(fp, "\tsource_%-12s\tsink_%-12s\tsum_%-12s", fNames[i], fNames[i],
            fNames[i]);
  }
  for (k = 0; k < nGridPoints; ++k) {
    fprintf(fp, "\n%-9E", x[k]);
    for (i = 0; i < nSpecies; ++i) {
      fprintf(fp, "\t%-9E\t%-9E\t%-9E", source[k][i] * 1.0e-3,
              sink[k][i] * 1.0e-3, (source[k][i] + sink[k][i]) * 1.0e-3);
    }
  }
  fclose(fp);
  Free2DArray(sink);
  Free2DArray(source);
}

template <class TransportModel>
void T1DSpecies<TransportModel>::PrintDiffusionCoeffs(TNewtonPtr bt) {
  int i, k;
  int nSpeciesInSystem = GetNSpeciesInSystem();
  char fileName[32];
  FILE *fp;
  TGridPtr grid = bt->GetGrid()->GetCurrentGrid();
  Double *x = grid->GetX()->vec;
  Double **D = fDiffusivity->mat;
  int nGridPoints = bt->GetCurrentGridPoints();
  Double left = bt->GetLeft();
  Double right = bt->GetRight();
  static int counter = 0;

  sprintf(fileName, "diff%d.dout", ++counter);
  if (counter >= 10) counter = 0;
  if (!(fp = fopen(fileName, "w"))) {
    cerr << "#warning: unable to open file " << fileName << NEWL;
    return;
  }
  fprintf(fp, "*\n");
  fprintf(fp, "%-12s", "eta");
  for (i = 0; i < nSpeciesInSystem; ++i) {
    fprintf(fp, "\tD_%-12s", fNames[i]);
  }
  fprintf(fp, "\n%-9E", left);
  for (i = 0; i < nSpeciesInSystem; ++i) {
    fprintf(fp, "\t%-9E", D[-1][i]);
  }
  for (k = 0; k < nGridPoints; ++k) {
    fprintf(fp, "\n%-9E", x[k]);
    for (i = 0; i < nSpeciesInSystem; ++i) {
      fprintf(fp, "\t%-9E", D[k][i]);
    }
  }
  fprintf(fp, "\n%-9E", right);
  for (i = 0; i < nSpeciesInSystem; ++i) {
    fprintf(fp, "\t%-9E", D[nGridPoints][i]);
  }
  fclose(fp);
}
