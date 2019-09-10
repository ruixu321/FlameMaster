#ifndef __TTRANS_FLAME_SOLVER_HPP__
#define __TTRANS_FLAME_SOLVER_HPP__

#undef DIFFCORR
#undef MOLARDIFF
#undef CONVECTION

// when creating unsteady libraries solution is written for every temperature change of UNSTLIBOUTDELTAT
#define UNSTLIBOUTDELTAT 40.0

// can't remember what this is
#undef ENTHALPYLIB 

#undef READCHI
#undef TIMEDEPCHI
#define ONECOLCHI
#undef READU

#undef LEWISCHANGE
#define LEWISSWITCHTIME 5.6e-3
#define LEWISSWITCHPERIOD 0.5e-3

#undef EQUILIBRIUM

// currently grid movement is central MOVEZR is upwind
#undef NEWCONVECTION
#define INGOOPT
#undef CENTRALZRCONV

#undef CENTRALGRID
#undef MOVEGRID   // Changed by Rui 12/17/2018

#undef MOVEZRIGHT

#undef TEMPFROMENTHALPY

#define NEWCHI
#undef TOTENT
#undef ROSSELAND
#undef NOSOOTRAD
#undef TESTNEWNUCLEATION

#undef SPARSEJAC

#ifdef MOVEGRID
#	undef SPARSEJAC
#endif

#undef LOGNORMALCHI    

#define NEWPROPERTIES
#define ENTFLUX
#define HEATCAPGRAD

#undef UPWIND
#define CENTRALTEMP

#define SOOTCONVECTION
#define SIZEDEPDIFFUSION
#define FRACMOMFACT 0.0

#define UPWINDSOOT
#undef SECORDSOOTUPWIND
#undef NEWSOOTTERMDISC
#define NEWTHERMOPHOR

#define UPWINDTHERMOPHOR
#define PRONE
#undef FIRST
#undef SECOND

#ifdef NEWSOOTTERMDISC
#	undef UPWINDSOOT
#endif

#undef DEBUGRES
#undef DEBUGFTOCMAT
#undef DEBUGGRID
#define DEBUGINITIAL
#undef DEBUGGETSOL
#undef DEBUGBAG
#undef DEBUGSOOTSOURCE
#define DEBUGINITIALDENS

#undef DPDT

#define DELTAPROG 1


#define DELTAINEW

static Double gasdev( int *idum );

template<typename Species>
void TTransFlameSolver<Species>::InitTTransFlameSolver( void )
{
	int	i;
#ifdef ENTFLUX
		fprintf( stderr, "ENTFLUX defined\n" );
#else
		fprintf( stderr, "ENTFLUX undefined\n" );
#endif
#ifdef HEATCAPGRAD
		fprintf( stderr, "HEATCAPGRAD defined\n" );
#else
		fprintf( stderr, "HEATCAPGRAD undefined\n" );
#endif
#ifdef DIFFCORR
		fprintf( stderr, "DIFFCORR defined\n" );
#else
		fprintf( stderr, "DIFFCORR undefined\n" );
#endif
#ifdef MOLARDIFF
		fprintf( stderr, "MOLARDIFF defined\n" );
#else
		fprintf( stderr, "MOLARDIFF undefined\n" );
#endif
#ifdef CONVECTION
		fprintf( stderr, "CONVECTION defined\n" );
#else
		fprintf( stderr, "CONVECTION undefined\n" );
#endif

	if ( fSoot ) {
		if ( fSoot->WithThermoPhoresis() ) {
			fprintf( stderr, "ThermoPhoresis defined\n" );
		}
		else {
			fprintf( stderr, "ThermoPhoresis undefined\n" );
		}

#ifdef FIRST
		fprintf( stderr, "FIRST defined\n" );
#else
		fprintf( stderr, "FIRST undefined\n" );
#endif

#ifdef SECOND
		fprintf( stderr, "SECOND defined\n" );
#else
		fprintf( stderr, "SECOND undefined\n" );
#endif

#ifdef MOVEGRID
		fprintf( stderr, "MOVEGRID defined\n" );
#else
		fprintf( stderr, "MOVEGRID undefined\n" );
#endif

#ifdef CONVECTION
		fprintf( stderr, "CONVECTION defined\n" );
#else
		fprintf( stderr, "CONVECTION undefined\n" );
#endif

#ifdef SOOTCONVECTION
		fprintf( stderr, "SOOTCONVECTION defined\n" );
#else
		fprintf( stderr, "SOOTCONVECTION undefined\n" );
#endif

#ifdef SECORDSOOTUPWIND
		fprintf( stderr, "SECORDSOOTUPWIND defined\n" );
#else
		fprintf( stderr, "SECORDSOOTUPWIND undefined\n" );
#endif

#ifdef SIZEDEPDIFFUSION
		fprintf( stderr, "SIZEDEPDIFFUSION defined\n" );
#else
		fprintf( stderr, "SIZEDEPDIFFUSION undefined\n" );
#endif

#ifdef UPWINDSOOT
		fprintf( stderr, "UPWINDSOOT defined\n" );
#else
		fprintf( stderr, "UPWINDSOOT undefined\n" );
#	ifdef NEWSOOTTERMDISC
		fprintf( stderr, "NEWSOOTTERMDISC defined\n" );
#	endif
#endif
	}
	
// set output file
	if ( !fInputData->fOutFileName || strcmp( fInputData->fOutFileName, "sterr" ) == 0 ) {
		fOutFilePtr = stderr;
		fprintf( stderr, "write output to 'stderr'\n" );
	}
	else if ( strcmp( fInputData->fOutFileName, "stdout" ) == 0 ) {
		fOutFilePtr = stdout;
		fprintf( stderr, "write output to 'stdout'\n" );
	}
	else {
		fOutFilePtr = GetOutfile( fInputData->fOutFileName, FileType::kText );
		fprintf( stderr, "write output to file '%s'\n\n"
			, GetOutfileName( fInputData->fOutFileName, FileType::kText ) );
	}

#ifdef LOGNORMALCHI
	SetRandomNumber();
#else
	fRandomNumber = 1.0;
#endif
	fFirstCall = TRUE;
	fMaxOrd = 5;
	fEquidistant = fInputData->fEquidistant;
	fKappa = fInputData->fKappa;
	fTauGrid = fInputData->fTauGrid;

	if (fInputData->fRadLossPercent > 0.0 && fInputData->fRadLossPercent < 1.0){
	fRadLossPercent = fInputData->fRadLossPercent;
	}
	else{
		fRadLossPercent = 1.0;
	}
	cerr << "fRadLossPercent =  " << fRadLossPercent << NEWL;

	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	if ( fInputData->fDeltaTMax > 0.0 ) {
		fDeltaTMax = fInputData->fDeltaTMax;
	}
	else {
		fDeltaTMax = 1.0e-4;
	}

	fNGridPoints = fInputData->fInitialGridPoints-2;

	fFirstTimeStep = 1.0e-09;

	int	nOfSpecies = fSpecies.GetNOfSpecies();
	int	nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();

	fRTol = New1DArray( 1 );
	fATol = New1DArray( 1 );
	fRTol[0] = ( fInputData->fTolRes > 0.0 ) ? fInputData->fTolRes : 1.0e-4;
	fATol[0] = ( fInputData->fTolDy > 0.0 ) ? fInputData->fTolDy : 1.0e-12;

#	ifdef LOGNORMALCHI
	fprintf( fOutFilePtr, "use gaussian distribution of scalar dissipation\n" );
#	endif
	
	int	NEQ = ( fNGridPoints + 2 ) * fNOfEquations;
	fNOfSolver = 1;
#ifdef SPARSEJAC
	fML = fNOfEquations;
#else
	fML = 2 * fNOfEquations;
#endif


	fNActualStep = new int*[fNOfSolver];
	fNActualOrd = new int*[fNOfSolver];
   	if ( !fNActualStep || !fNActualOrd ) {
		FatalError( "memory allocation of TTransFlameSolver failed" );
	}
	fTime = NewVector( fNOfSolver );
	fSolution = NewMatrix( fNOfEquations, fNGridPoints+2, kColumnPointers );
	fSolution->mat = &fSolution->mat[kNext];
	fSolPrime = NewMatrix( fNOfEquations, fNGridPoints+2, kColumnPointers );
	fSolPrime->mat = &fSolPrime->mat[kNext];
#ifdef CVODE
	fprintf(stderr, "\nUse CVODE solver\n\n");
	fNActualStep[0] = &fNActualStepCV;
	fNActualOrd[0] = &fNActualOrdCV;
	fActualTimeStepSize = &fActualTimeStepSizeCV;
#else
// workspace for ddassl
	fprintf(stderr, "\nUse DASSL solver\n\n");
	fLRW =  40 + ( fMaxOrd + 4 ) * NEQ 
				+ ( 2 * fML + fML + 1 ) * NEQ + 2 * ( NEQ / ( fML + fML + 1 ) + 1 );
	fLIW = 20 + NEQ;

	fDasslNEq = NEQ;
	fInfo = new IntVectorPtr[fNOfSolver];
	fRWork = new VectorPtr[fNOfSolver];
	fIWork = new IntVectorPtr[fNOfSolver];
	for ( i = 0; i < fNOfSolver; ++i ) {
		fInfo[i] = NewIntVector( 16 );
		fRWork[i] = NewVector( fLRW );
		fIWork[i] = NewIntVector( fLIW );
		fNActualStep[i] = &fIWork[i]->vec[kF11];
		fNActualOrd[i] = &fIWork[i]->vec[kF8];
	}
	InitDassl( fNOfSolver );
	if ( fInfo[0]->vec[kF10] == 1 ) {
		fprintf( fOutFilePtr, "###attention: clip negative concentrations\n" );
	}
	fActualTimeStepSize = &fRWork[0]->vec[kF7];
#endif
	
// other workspace
	
	fSolTime = NewVector( fNGridPoints+2 );
	fSolMassFracs = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	if ( fSoot ) {
		fSolSootMoments = NewMatrix( fSoot->GetNSootMoments(), fNGridPoints+2, kColumnPointers );
		fSolOldSootMoments = NewMatrix( fSoot->GetNSootMoments(), fNGridPoints+2, kColumnPointers );
		fSootMomentsWork = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );

		fSolSootMoments->mat = &fSolSootMoments->mat[kNext];
		fSolOldSootMoments->mat = &fSolOldSootMoments->mat[kNext];
		fSootMomentsWork->mat = &fSootMomentsWork->mat[kNext];
	}
	else {
		fSolSootMoments = NULL;
		fSolOldSootMoments = NULL;
		fSootMomentsWork = NULL;
	}
	fSolTemp = NewVector( fNGridPoints+2 );
	fSolOldTime = NewVector( fNGridPoints+2 );
	fSolOldMassFracs = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	fSolOldTemp = NewVector( fNGridPoints+2 );
	fSolGrid = NewVector( fNGridPoints+2 );
	fFDWCurr = New1DArray( fNGridPoints );
	fFDWPlus = New1DArray( fNGridPoints );
	fFDWMinus = New1DArray( fNGridPoints );
	fWCurr = New1DArray( fNGridPoints );
	fWPlus = New1DArray( fNGridPoints );
	fWMinus = New1DArray( fNGridPoints );
    fh = New1DArray( fNGridPoints+2 );
    fhm = New1DArray( fNGridPoints+2 );	
	fMonFct = New1DArray( fNGridPoints+2 );
	fgpDens = NewVector( fNGridPoints+2 );
	fMassFracsWork = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	fTempWork = NewVector( fNGridPoints+2 );
	fOutSolWork = NewVector( fNGridPoints+2 );
	fTotEntStart = NewVector( fNGridPoints+2 );
	fTotEntEnd = NewVector( fNGridPoints+2 );

	fHeatCpMix = NewVector( fNGridPoints+2 );
	fViscosity = NewVector( fNGridPoints+2 );
	fMolarMassMix = NewVector( fNGridPoints+2 );
	fDensity = NewVector( fNGridPoints+2 );
	fLambdaOverCpMix = NewVector( fNGridPoints+2 );

	fSpecHeatCp = NewMatrix( nOfSpeciesIn, fNGridPoints+2, kColumnPointers );
	fSpecEnthalpy = NewMatrix( nOfSpeciesIn, fNGridPoints+2, kColumnPointers );
	fProdRate = NewMatrix( nOfSpeciesIn, fNGridPoints+2, kColumnPointers );
	fSpecSource = NewMatrix( nOfSpeciesIn, fNGridPoints+2, kColumnPointers );  // Rui 

	fDiffTermY = NewVector( nOfSpeciesIn );
	fDiffTermW = NewVector( nOfSpeciesIn );
	fDiffCorrY = NewVector( nOfSpeciesIn );
	fDiffCorrW = NewVector( nOfSpeciesIn );

	fSolTime->vec = &fSolTime->vec[kNext];
	fSolMassFracs->mat = &fSolMassFracs->mat[kNext];
	fSolTemp->vec = &fSolTemp->vec[kNext];
	fSolOldTime->vec = &fSolOldTime->vec[kNext];
	fSolOldMassFracs->mat = &fSolOldMassFracs->mat[kNext];
	fSolOldTemp->vec = &fSolOldTemp->vec[kNext];
	fSolGrid->vec = &fSolGrid->vec[kNext];
    fh = &fh[kNext];
    fhm = &fhm[kNext];
    fgpDens->vec = &fgpDens->vec[kNext];
    fMonFct = &fMonFct[kNext];
	fMassFracsWork->mat = &fMassFracsWork->mat[kNext];
	fTempWork->vec = &fTempWork->vec[kNext];
	fTotEntStart->vec = &fTotEntStart->vec[kNext];
	fTotEntEnd->vec = &fTotEntEnd->vec[kNext];

	fHeatCpMix->vec = &fHeatCpMix->vec[kNext];
	fViscosity->vec = &fViscosity->vec[kNext];
	fMolarMassMix->vec = &fMolarMassMix->vec[kNext];
	fDensity->vec = &fDensity->vec[kNext];
	fLambdaOverCpMix->vec = &fLambdaOverCpMix->vec[kNext];

	fSpecHeatCp->mat = &fSpecHeatCp->mat[kNext];
	fSpecEnthalpy->mat = &fSpecEnthalpy->mat[kNext];
	fProdRate->mat = &fProdRate->mat[kNext];

	fSpecSource->mat = &fSpecSource->mat[kNext];   // Rui


	fZRin = NewVector( 10 );
	for ( i = 0; i < fZRin->len; ++i ) {
		fZRin->vec[i] = 1.0;
	}

	fMaxVals = NewVector( fNOfEquations );

	if ( fSoot ) {
		fSoot->SetMomentsOffset( fSootMoments );
	}

	// set variable names
	fVariableNames = new String[fNOfEquations];

	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );

	fVariableNames[fGrid] = new char[5];
	strcpy( fVariableNames[fGrid], "Grid" );

	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies.GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies.GetNames()[i] );
	}
	if ( fSoot ) {
		int	offset = fSoot->GetOffsetSootMoments();
		for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
			fVariableNames[offset + i] = new char[8];
			sprintf( fVariableNames[offset + i], "M%d", i );
		}
	}
	
	FILE *fpNames = GetOutfile( "VarNames", FileType::kText );
	for ( i = 0; i < fNOfEquations; ++i ) {
		fprintf( fpNames, "%s\n", fVariableNames[i] );
	}
	fclose( fpNames );

	fpNames = GetOutfile( "WanNames", FileType::kText );
	for ( i = 0; i < fNOfEquations; ++i ) {
		fprintf( fpNames, "      VANAMES(%d)='%-20s'\n", i+1, fVariableNames[i] );
	}
	fclose( fpNames );

	CompLewisNumbers( fInputData->fLewisNumberFile );
	fflush( fOutFilePtr );

	char	buff[128];
	char	*theFile;
	int		count = 0;
	FILE	*fpCA = NULL;
	Flag	isNew = FALSE;
	do {
		sprintf( buff, "CA%d.in", count );
		theFile = GetOutfileName( buff, FileType::kNone );
		fpCA = fopen( theFile, "r" );
		if ( fpCA ) {
			fclose( fpCA );
			++count;
		}
		else {
			isNew = TRUE;
		}
	} while( !isNew && count < 50 );
	
	if ( !isNew ) {
		fprintf( stderr, "cannot open CA.in file ( count = %d )\n", count );
		exit(2);
	}
	else {
		fCAinFile = fopen( theFile, "w" );
		fprintf( fCAinFile, "RPM = -1\n" );
		fprintf( fCAinFile, "VarsIn = 6\n" );
		fprintf( fCAinFile, "Time(s)	Pressure(Pa)	TOx(K)	TFuel(K)	Sci(1/s)	ZR\n" );
		fflush( fCAinFile );
	}

#ifdef READCHI
#	ifdef TIMEDEPCHI
	fprintf( stderr, "### use scalar dissipation rate from input file\n" );
	Double	*ZIn = New1DArray( 1000 );
	Double	*TimeIn = New1DArray( 1000 );
	char 	dummy[128];
	int		conv, j;

// read Z
	FILE *fpZIn = fopen( "zi.dat", "r" );
	if ( !fpZIn ) {
		fprintf( stderr, "#error: couldn't open input file 'zi.dat'\n" );
		exit( 2 );
	}
	for ( j = 0; j < 1000; ++j ) {
		conv = fscanf( fpZIn, "%lg", &ZIn[j] );
		if ( !conv || conv == EOF ) {
			break;
		}
	}
	if ( j == 0 ) {
		fprintf( stderr, "### error reading file 'zi.dat'\n" );
		exit( 2 );
	}
	fZIn = NewVector( j );

	for ( j = 0; j < fZIn->len; ++j ) {
		fZIn->vec[j] = ZIn[j];
	}
	fclose(fpZIn);

// read Z
	FILE *fpTimeIn = fopen( "Time.dat", "r" );
	if ( !fpTimeIn ) {
		fprintf( stderr, "#error: couldn't open input file 'Time.dat'\n" );
		exit( 2 );
	}
	for ( j = 0; j < 1000; ++j ) {
		conv = fscanf( fpTimeIn, "%lg", &TimeIn[j] );
		if ( !conv || conv == EOF ) {
			break;
		}
	}
	if ( j == 0 ) {
		fprintf( stderr, "### error reading file 'Time.dat'\n" );
		exit( 2 );
	}
	fTimeIn = NewVector( j );

	for ( j = 0; j < fTimeIn->len; ++j ) {
		fTimeIn->vec[j] = TimeIn[j];
	}
	fclose(fpTimeIn);

// read Chi
	fChiIn = NewMatrix( fZIn->len, fTimeIn->len, kColumnPointers );
	Double	**chiIn = fChiIn->mat;
	FILE *fpChiCount = fopen( "chi.dat", "r" );
	if ( !fpChiCount ) {
		fprintf( stderr, "#error: couldn't open input file 'chi.dat'\n" );
		exit( 2 );
	}

#		ifdef ONECOLCHI
	for ( i = 0; i < fTimeIn->len; ++i ) {
		for ( j = 0; j < fZIn->len; ++j ) {
			conv = fscanf( fpChiCount, "%lg", &chiIn[i][j] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong at i=%d j=%d\n", i, j );
				fprintf( stderr, "iLen actually %d\tjLen actually %d\n", fTimeIn->len, fZIn->len );
				break;
			}
		}
	}
#		else
	for ( j = 0; j < fZIn->len; ++j ) {
		for ( i = 0; i < fTimeIn->len; ++i ) {
			conv = fscanf( fpChiCount, "%lg", &chiIn[i][j] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong at i=%d j=%d\n", i, j );
				break;
			}
		}
	}
#		endif
	fclose( fpChiCount );

/* extrapolate chi to ZMean + 2 sqrt(ZVar)*/
	// first read ZR
	VectorPtr	ZRinVec = NewVector( fTimeIn->len );
	Double		*zrin = ZRinVec->vec;
	FILE 		*fpZRin = fopen( "ZR.in", "r" );
	if ( !fpZRin ) {
		fprintf( stderr, "#error: couldn't open input file 'ZR.in' -> don't extrapolate\n" );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			zrin[i] = 1.0;
		}
	}
	else {
		for ( i = 0; i < fTimeIn->len; ++i ) {
			conv = fscanf( fpZRin, "%lg", &zrin[i] );
			zrin[i] = MIN( 1.0, zrin[i] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong at i=%d j=%d\n", i );
				break;
			}
		}
	}
	
	int zrind, zc;
	for ( i = 0; i < fTimeIn->len; ++i ) {
		// get last non-zero point
		j = fZIn->len-2;
		while ( j > 0 && chiIn[i][j] < 1.0e-8 ) --j;
		// now j points at last nonzero scalar diss rate
		// go to next if zr is not larger than first zero chi mixture fraction
		if ( j == 0 || fZIn->vec[j+1] >= zrin[i]-1.0e-8  ) {
			continue;
		}
		// then get zrpoint 
		zrind = j+1;
		while ( fZIn->vec[zrind] < zrin[i]-1.0e-8 ) ++zrind;
		// now zrind points at first z larger than zr
		for ( zc = j+1; zc < zrind; ++zc ) {
			chiIn[i][zc] = chiIn[i][j] * ( 1.0 - ( fZIn->vec[zc] - fZIn->vec[j] ) / ( zrin[i] - fZIn->vec[j] ) );
		}
	}
	
	

	FILE	*fpZR = GetOutfile( "ZR", FileType::kData );
	fprintf( fpZR, "*\ntime\tZR\n" );

	for ( i = 0; i < fTimeIn->len; ++i ) {
		j = fZIn->len-2;
		while ( chiIn[i][j] < 1.0e-8 ) --j;
		fprintf( fpZR, "%g\t%g\n", fTimeIn->vec[i], fZIn->vec[j+1] );
	}
	fclose( fpZR );

	FILE	*fp = GetOutfile( "ChiIn", FileType::kData );
	fprintf( fp, "*\nZ" );
	for ( i = 0; i < fTimeIn->len; ++i ) {
		fprintf( fp, "\tt%g", fTimeIn->vec[i] );
	}
	fprintf( fp, "\n" );

	for ( j = 0; j < fZIn->len; ++j ) {
		fprintf( fp, "%g", fZIn->vec[j] );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			fprintf( fp, "\t%g", chiIn[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );




	Free1DArray( ZIn );
	Free1DArray( TimeIn );
#	else
	Double	*ZIn = New1DArray( 1000 );
	Double	*chiIn = New1DArray( 1000 );
	char 	dummy[128];
	int		conv, j;

	FILE *fpChiCount = fopen( "Chi.tout", "r" );
	if ( !fpChiCount ) {
		fprintf( stderr, "#error: couldn't open input file 'Chi.tout'\n" );
		exit( 2 );
	}

	fscanf( fpChiCount, "%s%s", dummy, dummy );
	for ( j = 0; j < 1000; ++j ) {
		conv = fscanf( fpChiCount, "%lg%lg", &ZIn[j], &chiIn[j] );
		if ( !conv || conv == EOF ) {
			break;
		}
	}
	fclose( fpChiCount );

	fZCount = NewVector( j );
	fChiCount = NewVector( j );
	for ( j = 0; j < fZCount->len; ++j ) {
		fZCount->vec[j] = ZIn[j];
		fChiCount->vec[j] = chiIn[j];
	}

#	endif
#endif

#ifdef READU
#	ifndef READCHI
	fprintf( stderr, "###error: can't use READU without READCHI\n" );
	exit( 2 );
#	endif
#	ifndef TIMEDEPCHI
	fprintf( stderr, "###error: can't use READU without TIMEDEPCHI\n" );
	exit( 2 );
#	endif
	fprintf( stderr, "### use axial velocity from input file\n" );

// read U
	fUstOverU = NewMatrix( fZIn->len, fTimeIn->len, kColumnPointers );
	Double	**UIn = fUstOverU->mat;
	FILE *fpUCount = fopen( "Ucond.dat", "r" );
	if ( !fpUCount ) {
		fprintf( stderr, "#error: couldn't open input file 'Ucond.dat'\n" );
		exit( 2 );
	}

	for ( i = 0; i < fTimeIn->len; ++i ) {
		for ( j = 0; j < fZIn->len; ++j ) {
			conv = fscanf( fpUCount, "%lg", &UIn[i][j] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong with UIn at i=%d j=%d\n", i, j );
				fprintf( stderr, "iLen actually %d\tjLen actually %d\n", fTimeIn->len, fZIn->len );
				break;
			}
			if ( UIn[i][j] <= 0.0 ) {
				if ( j == 0 ) {
					if ( i > 0 ) {
						UIn[i][j] = UIn[i-1][j];
					}
					else {
						// check next value
						++j;
						conv = fscanf( fpUCount, "%lg", &UIn[i][j] );
						if ( !conv || conv == EOF ) {
							fprintf( stderr, "something wrong with UIn at i=%d j=%d\n", i, j );
							fprintf( stderr, "iLen actually %d\tjLen actually %d\n", fTimeIn->len, fZIn->len );
							break;
						}
						if ( UIn[i][j] <= 0.0 ) {
							fprintf( stderr, "###error: something wrong with Uin 1\n" );
							exit( 2 );
						}
						else {
							UIn[i][j-1] = UIn[i][j];
						}
					}
				}
				else {
					UIn[i][j] = UIn[i][j-1];
				}
			}
		}
	}
	fclose( fpUCount );

	fp = GetOutfile( "UIn", FileType::kData );
	fprintf( fp, "*\nZ" );
	for ( i = 0; i < fTimeIn->len; ++i ) {
		fprintf( fp, "\tt%g", fTimeIn->vec[i] );
	}

	for ( j = 0; j < fZIn->len; ++j ) {
		fprintf( fp, "%g", fZIn->vec[j] );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			fprintf( fp, "\t%g", UIn[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );
	
#endif

#ifdef DELTAINEW
	fTempGSave = New1DArray( fNGridPoints+2 );
    fYGSave = New2DArray( fNGridPoints+2, nOfSpeciesIn );
    fDeltaI = New2DArray( fNGridPoints+2, nOfSpeciesIn );
	fG_ij = New3DArray( fNGridPoints+2, nOfSpeciesIn, nOfSpeciesIn );
	fTempGSave = &fTempGSave[kNext];
    fYGSave = &fYGSave[kNext];
    fDeltaI = &fDeltaI[kNext];
	fG_ij = &fG_ij[kNext];
#endif


}

template<typename Species>
TTransFlameSolver<Species>::~TTransFlameSolver( void )
{
#ifdef DELTAINEW
	fG_ij = &fG_ij[kPrev];
	fDeltaI = &fDeltaI[kPrev];
	fYGSave = &fYGSave[kPrev];
	fTempGSave = &fTempGSave[kPrev];
	Free3DArray( fG_ij );
	Free2DArray( fDeltaI );
	Free2DArray( fYGSave );
	Free1DArray( fTempGSave );
#endif

#ifdef READU
	DisposeMatrix( fUstOverU );
#endif

#ifdef READCHI
#	ifdef TIMEDEPCHI
	DisposeMatrix( fChiIn );
	DisposeVector( fZIn );
	DisposeVector( fTimeIn );
#	else
	DisposeVector( fChiCount );
	DisposeVector( fZCount );
#	endif
#endif

	int nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();

	for ( int i = 0; i < nOfSpeciesIn+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
	
	fMonFct = &fMonFct[kPrev];
	fhm = &fhm[kPrev];
	fh = &fh[kPrev];
	fgpDens->vec = &fgpDens->vec[kPrev];
	fTempWork->vec = &fTempWork->vec[kPrev];
	fMassFracsWork->mat = &fMassFracsWork->mat[kPrev];
	fSolGrid->vec = &fSolGrid->vec[kPrev];
	fSolOldTemp->vec = &fSolOldTemp->vec[kPrev];
	fSolOldMassFracs->mat = &fSolOldMassFracs->mat[kPrev];
	fSolOldTime->vec = &fSolOldTime->vec[kPrev];
	fSolTemp->vec = &fSolTemp->vec[kPrev];
	fSolMassFracs->mat = &fSolMassFracs->mat[kPrev];
	fSolTime->vec = &fSolTime->vec[kPrev];
	fTotEntStart->vec = &fTotEntStart->vec[kPrev];
	fTotEntEnd->vec = &fTotEntEnd->vec[kPrev];

	fHeatCpMix->vec = &fHeatCpMix->vec[kPrev];
	fViscosity->vec = &fViscosity->vec[kPrev];
	fMolarMassMix->vec = &fMolarMassMix->vec[kPrev];
	fDensity->vec = &fDensity->vec[kPrev];
	fLambdaOverCpMix->vec = &fLambdaOverCpMix->vec[kPrev];

	fSpecHeatCp->mat = &fSpecHeatCp->mat[kPrev];
	fSpecEnthalpy->mat = &fSpecEnthalpy->mat[kPrev];
	fProdRate->mat = &fProdRate->mat[kPrev];

	DisposeVector( fgpDens );
	Free1DArray( fMonFct );
	Free1DArray( fhm );
	Free1DArray( fh );
	DisposeVector( fOutSolWork );
	DisposeVector( fTempWork );
	DisposeMatrix( fMassFracsWork );
	Free1DArray( fWMinus );
	Free1DArray( fWPlus );
	Free1DArray( fWCurr );
	Free1DArray( fFDWMinus );
	Free1DArray( fFDWPlus );
	Free1DArray( fFDWCurr );
	DisposeVector( fSolGrid );
	DisposeVector( fSolOldTemp );
	DisposeMatrix( fSolOldMassFracs );
	DisposeVector( fSolOldTime );
	DisposeVector( fSolTemp );
	DisposeVector( fTotEntStart );
	DisposeVector( fTotEntEnd );
	DisposeVector( fZRin );
	DisposeVector( fDiffCorrW );
	DisposeVector( fDiffCorrY );
	DisposeVector( fDiffTermW );
	DisposeVector( fDiffTermY );
	DisposeMatrix( fProdRate );
	DisposeMatrix( fSpecEnthalpy );
	DisposeVector( fLambdaOverCpMix );
	DisposeVector( fDensity );
	DisposeVector( fMolarMassMix );
	DisposeVector( fHeatCpMix );
	DisposeVector( fViscosity );
	if ( fSoot ) {
		fSootMomentsWork->mat = &fSootMomentsWork->mat[kPrev];
		fSolOldSootMoments->mat = &fSolOldSootMoments->mat[kPrev];
		fSolSootMoments->mat = &fSolSootMoments->mat[kPrev];
		DisposeMatrix( fSootMomentsWork );
		DisposeMatrix( fSolOldSootMoments );
		DisposeMatrix( fSolSootMoments );
	}
	DisposeMatrix( fSolMassFracs );
	DisposeVector( fSolTime );

	fSolPrime->mat = &fSolPrime->mat[kPrev];
	fSolution->mat = &fSolution->mat[kPrev];
	DisposeMatrix( fSolPrime );
	DisposeMatrix( fSolution );

#ifndef CVODE
	DisposeIntVector( fIWork[0] );
	DisposeVector( fRWork[0] );
	DisposeIntVector( fInfo[0] );

	delete fIWork;
	delete fRWork;
	delete fInfo;
#else 
  N_VDestroy(fCVY);       /* Free the fCVY vector       */
  CVodeFree(&fMem);       /* Free the integrator memory */
#if SUNDIALS_VERSION_MAJOR > 2
  SUNLinSolFree(LS);      /* Free linear solver memory  */
  SUNMatDestroy(A);       /* Free the matrix memory     */
#endif // SUNDIALS_VERSION_MAJOR > 2
#endif

	delete fNActualOrd;
	delete fNActualStep;

	DisposeVector( fTime );

	Free1DArray( fATol );
	Free1DArray( fRTol );
}

template<typename Species>
void TTransFlameSolver<Species>::SolutionToCVode( void )
{
	copy( ( fNGridPoints + 2 ) * fNOfEquations, fSolution->mat[kPrev], 1, fCYdata, 1 );
}

template<typename Species>
void TTransFlameSolver<Species>::SolutionFromCVode( void )
{
	copy( ( fNGridPoints + 2 ) * fNOfEquations, fCYdata, 1, fSolution->mat[kPrev], 1 );
}

template<typename Species>
int TTransFlameSolver<Species>::InitCVODE( void )
{
	int flag;  // for checking return values from cvode functions
	int	NEQ = ( fNGridPoints + 2 ) * fNOfEquations;

  /* set the solution vector and initialise from solution YInit
     also assign pointer to first memory location in Y
     NOTE: use compiler directive here later for parallel */
	fCVY = N_VMake_Serial(NEQ, fSolution->mat[kPrev]);

	fCYdata = NV_DATA_S(fCVY);

  /* create the memory object for CVODE:
     The system is stiff due to chemistry, thus choose solver appropriately
       CV_BDF:    linear multistep method (Backward Differentiation Formulas) 
       CV_NEWTON: Newton iteration (modified newton iteration for banded) */
	fMem = CVodeCreate(CV_BDF, CV_NEWTON);

  /* allocate internal memory object:
     initialise time:tInit
     RHS function:   cvodeRHS
     CV_SS:          scalar relative and absolute tolerances
     reltol:         relative tolerance
     abstol:         absolute tolerance */

	CVodeInit(fMem, ResTransFlameImpliSolverCV<TTransFlameSolver<Species> >, fTStart, fCVY);
	flag = CVodeSStolerances(fMem, fRTol[0], fATol[0]);

	if (flag != CV_SUCCESS) {
      if (flag == CV_MEM_NULL)
		printf("CVODE ERROR: CVodeMalloc requires CVodeCreate");
		printf("exiting: unable to allocate with CVodeMalloc\n");
		exit(1);
	}

  // set pointer to parameters (structs) needed by the RHS function
  flag = CVodeSetUserData(fMem, this);
  
  // defined dense matrix, set number of equations
#if SUNDIALS_VERSION_MAJOR > 2
    /* Create banded SUNMatrix for use in linear solves -- since this will be factored, 
     set the storage bandwidth to be the sum of upper and lower bandwidths */
  A = SUNBandMatrix(NEQ, fML, fML, 2*fML);
  //if(check_flag((void *)A, "SUNBandMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver object for use by CVode */
  LS = SUNBandLinearSolver(fCVY, A);
  //if(check_flag((void *)LS, "SUNBandLinearSolver", 0)) return(1);
  
  /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
  flag = CVDlsSetLinearSolver(fMem, LS, A);
  //if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);
#else 
  flag = CVBand(fMem, NEQ, fML, fML);
#endif // SUNDIALS_VERSION_MAJOR > 2

  // options
  CVodeSetInitStep(fMem, fFirstTimeStep );
  CVodeSetMaxStep(fMem, fDeltaTMax );
  
  return flag;
}

template<typename Species>
void TTransFlameSolver<Species>::SetInfo( int i )
{
	int		*info = fInfo[i]->vec;

	info[kF1] = 0;	// first call
	info[kF2] = 0;	// RTOL, ATOL are scalars
	info[kF3] = 1;	// solution at TOUT, no intermediates
	info[kF4] = 1;	// integration can be carried out beyond TOUT ( and interpolated )
	info[kF5] = 0;	// numerical jacobian
	info[kF6] = 1;	// full ( not banded ) jacobian
	info[kF7] = 1;	// no choice for maximum stepsize
	info[kF8] = 1;	// no choice for first stepsize
	if ( fMaxOrd == 5 ) {
		info[kF9] = 0;	// choose default for MAXORD ( = 5 )
	}
	else {
		info[kF9] = 1;	// don't choose default for MAXORD ( = 5 )
	}
	info[kF10] = 1;	// solve without non negative constraint // this is clip
	info[kF11] = 1;	// initial values are consistent
}

template<typename Species>
void TTransFlameSolver<Species>::InitDassl( int nGridPoints )
{

	for ( int i = 0; i < nGridPoints; ++i ) {
		SetInfo( i );
		InitDasslWorkSpace( i );
	}
}

template<typename Species>
void TTransFlameSolver<Species>::InitDasslWorkSpace( int i )
{
	int		*info = fInfo[i]->vec;
	int		*iWork = fIWork[i]->vec;
	Double	*rWork = fRWork[i]->vec;

	if ( info[kF6] == 1 ) {
		iWork[kF1] = fML;
		iWork[kF2] = fML;		
	}
	if ( info[kF8] == 1 ) {
		rWork[kF3] = fFirstTimeStep;
	}
	if ( info[kF9] == 1 ) {
		iWork[kF3] = fMaxOrd;
	}
	if ( info[kF7] == 1 ) {
		rWork[kF2] = fDeltaTMax;
	}
}

template<typename Species>
void TTransFlameSolver<Species>::SetMaxTimeStep( int kAct, Double /*t*/ )
{
	fRWork[kAct]->vec[kF2] = fDeltaTMax;
}

template<typename Species>
void TTransFlameSolver<Species>::GetSpeciesData( char **newNames, MatrixPtr highMat
					, MatrixPtr lowMat, Double *molarMass )
{
	int		ind;
	Double	*theMolarMass = fSpecies.GetMolarMass()->vec;
	Double	**high = highMat->mat;
	Double	**low = lowMat->mat;
	Double	**theHigh = fSpecies.GetCoeffHigh()->mat;
	Double	**theLow = fSpecies.GetCoeffLow()->mat;
	
	for ( int i = 0; i < highMat->cols; ++i ) {
		if ( ( ind = fSpecies.FindSpecies( newNames[i] ) ) < 0 ) {
			fprintf( fOutFilePtr, "%s'%s'\n", "#error in function 'GetSpeciesData': no species ", newNames[i]  );
			exit( 2 );
		}
		copy( NOFCOEFFS, theLow[ind], 1, low[i], 1 );
		copy( NOFCOEFFS, theHigh[ind], 1, high[i], 1 );
		molarMass[i] = theMolarMass[ind];
	}
}

template<typename Species>
void TTransFlameSolver<Species>::Initialize( Double timeStart
					, ConstStringArray names
					, Double **startSolution, int vars
					, Double *grid, int gridPointsA
					, Double pressureStart, Double scalarDissRateStart
					, Double firstTimeStep
					, Double ZRef, Double ZRStart )
{
	int 	k;
	Double	*t = fTime->vec;
	
	if ( firstTimeStep > 1.0e-13 ) {
		fprintf( fOutFilePtr, "initial timestep is %g s\n", firstTimeStep );
		fFirstTimeStep = firstTimeStep;
#ifndef CVODE
		for ( k = 0; k < fNOfSolver; ++k ) {
			fRWork[k]->vec[kF3] = firstTimeStep;
		}
#endif
	}
	fTStart = fTEnd = timeStart;
	if ( timeStart > 0.0 ) {
		fprintf( fOutFilePtr, "start at time %g s\n", timeStart );
	}
	for ( k = 0; k < fNOfSolver; ++k ) {
		t[k] = timeStart;
	}
	fPressStart = fPressEnd = pressureStart;
	fChiStart = fChiEnd = scalarDissRateStart;

	fZRStart = fZREnd = ZRStart;
	
	fZl = grid[0];
	fZr = grid[gridPointsA-1];

	if ( gridPointsA != fSolGrid->len ) {
		fprintf( stderr, "make grid, nGIn = %d\tnG = %d\n", gridPointsA, fSolGrid->len );
		MakeGrid( &fSolGrid->vec[kPrev], fSolGrid->len, fZl, fZr, fEquidistant );
	}
	else{
		fprintf( stderr, "use inputgrid, nGIn = %d\tnG = %d\n", gridPointsA, fSolGrid->len );
		Double	*newGrid = &fSolGrid->vec[kPrev];
		for ( k = 0; k < gridPointsA; ++k ) {
			newGrid[k] = grid[k];
		}
	}
	SetWeights();

	SetInitial( names, startSolution, grid, gridPointsA, vars );
	
#ifdef CVODE
	InitCVODE();
#endif

	fTempOxEnd = fSolution->mat[kPrev][fTemperature];
	fTempFuelEnd = fSolution->mat[fNGridPoints][fTemperature];

	fprintf( fCAinFile, "%g\t%g\t%g\t%g\t%g\t%g\n", fTStart, fPressStart
				, fTempOxEnd, fTempFuelEnd, fChiStart, fZRStart );
	fflush( fCAinFile );

	if ( ZRef < 0.0 ) {
		fZRef = GetZStoi();
		fprintf( fOutFilePtr, "ZStoi = %g\n", fZRef );
	}
	else {
		fZRef = ZRef;
	}
	if ( fZRef < grid[0] || fZRef > grid[gridPointsA-1] ) {
		fprintf( fOutFilePtr, "###error: ZRef = %g out of bounds (Zmin = %g Zmax = %g)\n"
				, fZRef, grid[0], grid[gridPointsA-1] );
		exit(2);
	}

#ifdef READU
	/* find Zstoi */

	int		i, j;
	Double	**UIn = fUstOverU->mat;
	for ( j = 0; j < fZIn->len; ++j ) {
		if ( fZIn->vec[j] > GetZStoi() ) break;
	}
	// now zstoi is between j-1 and j
	int zloc = j;
	Double	ustoi;

	for ( i = 0; i < fTimeIn->len; ++i ) {
		ustoi = Interpol( GetZStoi(), UIn[i][zloc-1], fZIn->vec[zloc-1], UIn[i][zloc-1], fZIn->vec[zloc] );
		for ( j = 0; j < fZIn->len; ++j ) {
			UIn[i][j] = ustoi / UIn[i][j];
		}
	}
	
	FILE	*fp = GetOutfile( "UstOverU", FileType::kData );
	fprintf( fp, "*\nZ" );
	for ( i = 0; i < fTimeIn->len; ++i ) {
		fprintf( fp, "\tt%g", fTimeIn->vec[i] );
	}

	for ( j = 0; j < fZIn->len; ++j ) {
		fprintf( fp, "%g", fZIn->vec[j] );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			fprintf( fp, "\t%g", UIn[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );
#endif

#ifdef EQUILIBRIUM
	int		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int		nndd = 7;
	Double	*xmol_i = fSpecies.GetMolarMass()->vec;
	Double	p = GetPressure( GetCurrentTime() ) * 1.0e-5;
	Double	*c0;
	char	*symbol = new char[(nSpeciesIn+4)*20];
	Double	*ponalo = new Double[nSpeciesIn*nndd];
	Double	*ponahi = new Double[nSpeciesIn*nndd];
	Double	totent;

	Double	**lowVals = fSpecies.GetCoeffLow()->mat;
	Double	**highVals = fSpecies.GetCoeffHigh()->mat;
	MatrixPtr	lowMat = FortranToCMat( ponalo, nndd, nndd, nSpeciesIn, nSpeciesIn );
	MatrixPtr	highMat = FortranToCMat( ponahi, nndd, nndd, nSpeciesIn, nSpeciesIn );
	Double	**low = lowMat->mat;
	Double	**high = highMat->mat;
	for ( int i = 0; i < nSpeciesIn; ++i ) {
		copy( nndd, lowVals[i], 1, low[i], 1 );
		copy( nndd, highVals[i], 1, high[i], 1 );
	}

	CToFortranCharArray( symbol, fSpecies.GetNames(), 20, nSpeciesIn );

	for ( k = -1; k <= fNGridPoints; ++k ) {
		c0 = &fSolution->mat[k][fFirstSpecies];
		totent = fTotEntStart->vec[k]*1.0e-3;
		ADIABFLAMETEMP( &fSolGrid->vec[k], &nSpeciesIn, c0, &fSolution->mat[k][fTemperature], &p, symbol, 
                                 xmol_i, &totent,ponalo, ponahi,
                                 &nndd,
                                 &fSolution->mat[k][fTemperature],c0 );
		SaveSolution( k, GetCurrentTime(), fSolution->mat[k] );
		SaveSolution( k, GetCurrentTime(), fSolution->mat[k] );
	}

	DisposeFToCMat( highMat );
	DisposeFToCMat( lowMat );
	fprintf( stderr, "write equilibrium data and exit\n" );
	WriteFlameletFileInitialize( NULL, NULL, "equi" );
	exit( 2 );
#endif
}

template<typename Species>
Double TTransFlameSolver<Species>::GetTotEnt( Double temp, Double *Y )
{
	int		i;
	int		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double	*h = fSpecies.GetEnthalpy()->vec;
	Double	totEnt;
	
	fSpecies.ComputeSpeciesProperties( temp );
	totEnt = 0.0;
	for ( i = 0; i < nSpeciesIn; ++i ) {
		totEnt += Y[i] * h[i];
	}
	
	return totEnt;
}

template<typename Species>
Flag TTransFlameSolver<Species>::Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, int deltaStepsOut )
{
	return Solve( timeEnd, pressureEnd, scalarDissRateEnd, temOxEnd, tempFuelEnd, 1.0, deltaStepsOut );
}

template<typename Species>
Flag TTransFlameSolver<Species>::Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd, int deltaStepsOut )
{
	Flag	error = FALSE;
	Flag	leave = FALSE;
	int		counter = 0;
		
	fprintf( fCAinFile, "%g\t%g\t%g\t%g\t%g\t%g\n", timeEnd, pressureEnd
				, temOxEnd, tempFuelEnd, scalarDissRateEnd, ZREnd );
	fflush( fCAinFile );

	SetEndValues( timeEnd, pressureEnd, scalarDissRateEnd, temOxEnd, tempFuelEnd, ZREnd );

#ifdef LOGNORMALCHI
	SetRandomNumber();
#endif
	
	do {
		leave = OneImpliStep();

		if ( deltaStepsOut > 0 && *fNActualStep[0] % deltaStepsOut == 0 ) {
			char	countstring[128];
			FILE	*fp;
			if ( counter > 0 ) {
				sprintf( countstring, "%d", counter );
				fp = GetOutputFile( fTime->vec[0], NULL, countstring, FileType::kText );
			}
			else {
				fp = GetOutputFile( fTime->vec[0], NULL, NULL, FileType::kText );
			}
			fprintf( fOutFilePtr, "dump output\n" );
			WriteFlameletFile( fp, NULL, NULL );
			fclose( fp );
			counter++;
		}
		error = leave;
		if ( fTime->vec[0] + MIN( ( fTEnd - fTStart ) * 1.0e-7, 1.0e-12 ) >= fTEnd ) {
			leave = TRUE;
			cerr << "dump libout" << NEWL;
			char   tl[128];
			Double         tNow = fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature];
			// sprintf( tl, "Chi%05g_t%05gms_dH%05g_Tmax%05g", GetRefDissRate( fTime->vec[0] ), GetCurrentTime()*1000, DefectMax/1000 , tNow);   // Rui 05142019
			// FILE *fp = GetOutfile( tl, FileType::kText );
			               // FILE *fp = GetOutputFile( fTime->vec[0], tl, NULL, FileType::kText );                           
			WriteFlameletFile( NULL, NULL, NULL );
			// fclose( fp );
		}
#	ifdef HP
		if ( gExit ) {
			DoExit();
		}
#	endif
	} while ( !leave );	

	return error;
}

template<typename Species>
void TTransFlameSolver<Species>::SetEndValues( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd )
{
	int	k;
	fprintf( fOutFilePtr, "ZREnd = %g\n", ZREnd );
	if ( fTEnd >= timeEnd ) {
		fprintf( stderr, "###Error: tStart (%g) >= tEnd (%g)\n", fTEnd, timeEnd );
		exit( 2 );
	}
	fTStart = fTEnd;
	fPressStart = fPressEnd;
	fChiStart = fChiEnd;
	fTempOxStart = fTempOxEnd;
	fTempFuelStart = fTempFuelEnd;
	fZRStart = fZREnd;
	
//set ZREnd
	Double sum = 0.0;
	for ( int i = fZRin->len-1; i > 0 ; --i ) {
		fZRin->vec[i] = fZRin->vec[i-1];
		sum += fZRin->vec[i];
	}
	fZRin->vec[0] = ZREnd;
	fZREnd = MIN( fZRStart, ( sum + fZRin->vec[0] ) / fZRin->len );

	Double	*totEntStart = fTotEntStart->vec;
	Double	*totEntEnd = fTotEntEnd->vec;
	Double	*Z = fSolGrid->vec;
	copy( fNGridPoints+2, &totEntEnd[kPrev], 1, &totEntStart[kPrev], 1 );
		
	fTEnd = timeEnd;
	fPressEnd = pressureEnd;
	fChiEnd = scalarDissRateEnd;

	fTempOxEnd = temOxEnd;
	fTempFuelEnd = tempFuelEnd;
	totEntEnd[kPrev] = GetTotEnt( fTempOxEnd, fSolMassFracs->mat[kPrev] );
	totEntEnd[fNGridPoints] = GetTotEnt( fTempFuelEnd, fSolMassFracs->mat[fNGridPoints] );
	for ( k = 0; k <= fNGridPoints-1; ++k ) {
		totEntEnd[k] = Interpol( Z[k], totEntEnd[kPrev], Z[kPrev], totEntEnd[fNGridPoints], Z[fNGridPoints] );
	}

	fprintf( fOutFilePtr, "solve flamelet\n" );
	fprintf( fOutFilePtr, "from t = %.4g ms with p = %.4g bar  chi = %.4g 1/s  Tox = %.4g K  Tfu = %.4g K  ZR = %.4g\n"
		, fTStart * 1.0e3, fPressStart / 1.0e5, fChiStart, fTempOxStart, fTempFuelStart, fZRStart );
	fprintf( fOutFilePtr, "  to t = %.4g ms with p = %.4g bar  chi = %.4g 1/s  Tox = %.4g K  Tfu = %.4g K  ZR = %.4g\n"
		, fTEnd * 1.0e3, fPressEnd / 1.0e5, fChiEnd, fTempOxEnd, fTempFuelEnd, fZREnd );

	fDTdtOx = ( fTempOxEnd - fTempOxStart ) / ( fTEnd - fTStart );
	fDTdtFu = ( fTempFuelEnd - fTempFuelStart ) / ( fTEnd - fTStart );

	fdDeltaZdt = ( fZREnd - fZRStart ) / ( fTEnd - fTStart );

#ifdef DPDT
	fDPdt = ( fPressEnd - fPressStart ) / ( fTEnd - fTStart );
#else
	fDPdt = 0.0;
#endif
	SetMolarMassOverRInf();

#ifdef CVODE
		  CVodeSetStopTime(fMem, fTEnd);
#else
	for ( k = 0; k < fNOfSolver; ++k ) {
		if ( fInfo[k]->vec[kF4] == 1 ) {
			fRWork[k]->vec[kF1] = fTEnd;
		}
	}
#endif
}

template<typename Species>
void TTransFlameSolver<Species>::SetMolarMassOverRInf( void )
{
	fProperties->ComputeMixtureMolarMass( fWOverRInf, fSolMassFracs->mat[kPrev]
				, fSpecies.GetMolarMass()->vec, fSpecies.GetNSpeciesInSystem() );
	fWOverRInf /= RGAS;
}

template<typename Species>
void TTransFlameSolver<Species>::PrintDissRate( Double t )
{
	static int	counter = 0;		
	char	name[64];
	sprintf( name, "DissRate%d", counter++ );
	FILE	*fp = GetOutfile( name, FileType::kData );

	fprintf( fp, "*\nZ\tChi\n" );
	for ( int k = -1; k <= fNGridPoints; ++k ) 
	{
		fprintf( fp, "%g\t%g\n", fSolGrid->vec[k], GetDissRate( t, fSolGrid->vec[k] ) );
	}
	
	fclose( fp );
}

template<typename Species>
Double TTransFlameSolver<Species>::GetTempOfEnthalpy( Double ent, Double *Y, Double initialGuess )
{
	Double	*h = fSpecies.GetEnthalpy()->vec;
	Double	*cp = fSpecies.GetHeatCapacity()->vec;
	Double	temp = ( initialGuess > 0.0 ) ? initialGuess : 1000;
	Double	deltaT;
	Double	entSum, cpSum;
	int		i;
	int		nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int		count = 0;

	do {
		fSpecies.ComputeSpeciesProperties( temp );
		for ( i = 0, entSum = cpSum = 0; i < nOfSpeciesIn; ++i ) {
			cpSum += cp[i] * Y[i];
			entSum += h[i] * Y[i];
		}
		
		deltaT = -( entSum - ent ) / cpSum;
		temp += deltaT;
		
		count++;
		if ( ++count > 1000 ) {
			fprintf( stderr, "#Error: temperature iteration for h = %g not converged\n"
					, ent );
			exit( 2 );
		}
	} while ( fabs( deltaT / temp ) > 1.0e-3 );

	return temp;
}

template<typename Species>
void TTransFlameSolver<Species>::SetInitial( ConstStringArray names, Double **startSol, Double *grid, int gridPointsA, int vars )
{
	int			i, k, ind;
	int			nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
	MMDataBag	bag( vars );
	
// copy input to fSolution
	bag.Initialize();
	bag.SetOldInpedVar( grid, gridPointsA, 1, "xIn" );
	for ( i = 0; i < vars; ++i ) {
		bag.Insert( &startSol[0][i], gridPointsA, startSol[1] - startSol[0], names[i] );
	}
	bag.SetNewInpedVar( &fSolGrid->vec[kPrev], fSolGrid->len, 1, "xNew" );

	Double		*newY = &fOutSolWork->vec[kNext];
	Double		**sol = fSolution->mat;
	for ( i = 0; i < bag.NumElems(); ++i ) {
		ind = GetVariableIndex( bag[i].Name() );
		if ( ind >= 0 ) {
			bag[i].Map( fOutSolWork );
			for ( k = -1; k < fSolGrid->len-1; ++k ) {
				sol[k][ind] = newY[k];
			}
			if ( fSoot && ind >= fSootMoments 
						&& ind < fSootMoments + fSoot->GetNSootMoments() ) {
			}
		}
		else {
			fprintf( fOutFilePtr, "#warning: no match for variable %s\n", bag[i].Name() );
		}
	}

	for ( k = -1; k < fSolGrid->len-1; ++k ) {
        sol[k][fGrid] = fSolGrid->vec[k];
	}

	Double	*totEntStart = fTotEntStart->vec;
	Double	*h = fSpecies.GetEnthalpy()->vec;
#ifdef TEMPFROMENTHALPY
// set Enthalpy
	Double	*solatk;

	// set total enthalpy
		// at k = -1
	totEntStart[kPrev] = 0.0;
	solatk = fSolution->mat[kPrev];
	fSpecies.ComputeSpeciesProperties( solatk[fTemperature] );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		totEntStart[kPrev] += solatk[fFirstSpecies+i] * h[i];
	}
		// at k = -1
	totEntStart[fNGridPoints] = 0.0;
	solatk = fSolution->mat[fNGridPoints];
	fSpecies.ComputeSpeciesProperties( solatk[fTemperature] );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		totEntStart[fNGridPoints] += solatk[fFirstSpecies+i] * h[i];
	}
		// linearly interpolate enthalpy
	for ( k = 0; k < fNGridPoints; ++k ) {
		solatk = fSolution->mat[k];
		totEntStart[k] = Interpol( fSolGrid->vec[k]
					, totEntStart[kPrev], fSolGrid->vec[kPrev]
					, totEntStart[fNGridPoints], fSolGrid->vec[fNGridPoints] );
		solatk[fTemperature] = GetTempOfEnthalpy( totEntStart[k], &solatk[fFirstSpecies], fSolution->mat[k-1][fTemperature] );
	}	
#endif

#ifdef DEBUGINITIAL
	FILE	*fp = GetOutfile( "InitialSolOld", FileType::kData );
	PrintSolution( fp, fSolGrid->len-2, fNOfEquations, fSolGrid->vec, fSolution->mat, fVariableNames );
	fclose( fp );
#endif

	if ( !fEquidistant ) {
	}

	if ( fSoot ) {
		// solver works with M/rho
		int		nSootMoments = fSoot->GetNSootMoments();
		int		momOff = fSoot->GetOffsetSootMoments();
		Double	mixMolarMass;
		Double	*solk;
		for ( k = -1; k <= fNGridPoints; ++k )
		{
			solk = sol[k];
			fProperties->ComputeMixtureMolarMass( mixMolarMass, &sol[k][fFirstSpecies], fSpecies.GetMolarMass()->vec, nOfSpeciesIn );
			for ( i = 0; i < nSootMoments; ++i ) {
				solk[momOff+i] *= RGAS * solk[fTemperature]
								/ ( GetPressure( fTime->vec[0] ) * mixMolarMass );
			}
		}
	}
	
#ifdef DEBUGINITIAL
	for ( k = -1; k < fSolGrid->len-1; ++k ) {
        fSolGrid->vec[k] = sol[k][fGrid];
	}
	fp = GetOutfile( "InitialSol", FileType::kData );
	PrintSolution( fp, fSolGrid->len-2, fNOfEquations, fSolGrid->vec, fSolution->mat, fVariableNames );
	fclose( fp );
#endif

// save initial solution to fSol
	int		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double	*Y;
	Double	*solVec;
	Double	pressure;
	Double	*totEntEnd = fTotEntEnd->vec;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		solVec = fSolution->mat[k];
		pressure = GetPressure( fTime->vec[0] );
		Y = fSolMassFracs->mat[k];
		// put fSolution to fMassFracs
		SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
		// update steady state concs
#ifdef DELTAINEW
		UpdateThermoProps( k, Y, fSolTemp->vec[k], pressure
					, fProperties->GetDensityRef(), kDensFromPress
					, (fSoot)?&solVec[fSoot->GetOffsetSootMoments()]:NULL );
#else 
		T0DFlame<Species>::UpdateThermoProps( Y, fSolTemp->vec[k], pressure
					, fProperties->GetDensityRef(), kDensFromPress
					, (fSoot)?&solVec[fSoot->GetOffsetSootMoments()]:NULL );
#endif
		fDensity->vec[k] = fProperties->GetDensity();
		// second call puts solution including steady state concs to  down to fMassFracs
		SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
		// third call puts solution down to fSolOld
		SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
		// set total enthalpy
		totEntStart[k] = 0.0;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			totEntStart[k] += Y[i] * h[i];
		}
		totEntEnd[k] = totEntStart[k];
	}

	SetMaxVals();
}

template<typename Species>
void TTransFlameSolver<Species>::SetGrid( void )
{
	int		i, k, countF = 0;
	Double	**sol = fSolution->mat;
	Double	*F = New1DArray( fNGridPoints+2 );
	Double	**solOld = New2DArray( fNGridPoints+2, fNOfEquations );
	Double	*Z = fSolGrid->vec;
	Double	deltaF, FCurr;

	F = &F[kNext];
	solOld = &solOld[kNext];

	for ( k = -1; k <= fNGridPoints; ++k ) {
		for ( i = 0; i < fNOfEquations; ++i ) {
			solOld[k][i] = sol[k][i];
		}
	}

	SetFDWeights( (Double *) *sol );
	ComputeGPDistributionWeights( (Double *) *sol );

// estimate weight fuction at Z=0
	fMonFct[kPrev] = fMonFct[0] - ( fMonFct[1] - fMonFct[0] ) / ( Z[1] - Z[0] ) * ( Z[0] - Z[kPrev] );


	TrapIntegrate( fNGridPoints+2, &fMonFct[kPrev], &Z[kPrev], &F[kPrev] );
	deltaF = F[fNGridPoints] / ( fNGridPoints + 1 );
	
	for ( k = 0; k < fNGridPoints; ++k ) {
		FCurr = ( k + 1 ) * deltaF;
		while ( countF < fNGridPoints && F[countF] < FCurr ) ++countF;
		for ( i = 0; i < fNOfEquations; ++i ) {
			sol[k][i] = Interpol( FCurr, solOld[countF-1][i], F[countF-1], solOld[countF][i], F[countF] );
		}
	}

#ifdef DEBUGINITIAL
	FILE	*fp = GetOutfile( "F", FileType::kData );
	fprintf( fp, "*\nZ\tG\tF\tTOld\tZNew\tTNew\n" );
	for ( k = -1; k <= fNGridPoints; ++k )
	{
		fprintf( fp, "%g\t%g\t%g\t%g\t%g\t%g\n", fSolGrid->vec[k], fMonFct[k], F[k], solOld[k][fTemperature]
								, sol[k][fGrid], sol[k][fTemperature] );
	}
	fclose( fp );
#endif


	Free2DArray(&solOld[kPrev]);
	Free1DArray(&F[kPrev]);
}

template<typename Species>
void TTransFlameSolver<Species>::GetSolution( ConstStringArray names, Double **outSol, Double *grid, int gridPointsA, int vars
									, Double *density )
{
	int			i, k, ind;
	int			nOfSpecies = fSpecies.GetNOfSpecies();
	int			nSootMoments;
	MMDataBag	bag( nOfSpecies+fVariablesWithoutSpecies+( ( density ) ? 1 : 0 ) );
	Double		**theY = fMassFracsWork->mat;
	Double		**theMom;
	int			theYOff = fMassFracsWork->rows;
	int			theMomOff, momOff;
	Double		*theTemp = fTempWork->vec;
	char		**specNames = fSpecies.GetNames();
	VectorPtr	fOutSol = NewVector( gridPointsA );
	VectorPtr	ZGridVec = NewVector( fSolGrid->len );
	Double		*ZGrid = &ZGridVec->vec[kNext];
	Flag		*outSet = new Flag[vars];
	Flag		densitySet = FALSE;
	for ( i = 0; i < vars; ++i ) {
		outSet[i] = FALSE;
	}

	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		momOff = fSoot->GetOffsetSootMoments();
		theMomOff = fSootMomentsWork->rows;
		theMom = fSootMomentsWork->mat;
	}
	
	Double	ZR = Interpol( GetCurrentTime(), fZRStart, fTStart, fZREnd, fTEnd );
	for ( i = -1; i < fSolGrid->len-1; ++i ) {
		ZGrid[i] = fSolGrid->vec[i] * ZR;
	}
	
// copy solution to bag
	bag.Initialize();
	bag.SetOldInpedVar( &ZGrid[kPrev], ZGridVec->len, 1, "xIn" );

	SetOutSolution();
	if ( density ) {
		int		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
		Double	MM;
		Double	*rho = fDensity->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			fProperties->ComputeMixtureMolarMass( MM, theY[k], fSpecies.GetMolarMass()->vec
						, nSpeciesIn );
			rho[k] = GetPressure() * MM / ( RGAS * theTemp[k] );
		}
	}
	bag.Insert( &theTemp[-1], fNGridPoints+2, 1, fVariableNames[fTemperature] );
	for ( i = 0; i < nOfSpecies; ++i ) {
		bag.Insert( &theY[-1][i], fNGridPoints+2, theYOff, specNames[i] );
	}
	if ( fSoot ) {
		for ( i = 0; i < nSootMoments; ++i ) {
			bag.Insert( &theMom[-1][i], fNGridPoints+2, theMomOff, fVariableNames[momOff+i] );
		}
	}
	if ( density ) {
		bag.Insert( &fDensity->vec[-1], fNGridPoints+2, 1, "density" );
	}
	
	bag.SetNewInpedVar( grid, gridPointsA, 1, "xNew" );

#ifdef DEBUGBAG
	cout << bag;
#endif

	Double		**newSol = &outSol[kNext];
	Double		*outWork = &fOutSol->vec[kNext];
	for ( i = 0; i < bag.NumElems(); ++i ) {
		ind = GetVariableIndex( bag[i].Name(), names, vars );
		if ( ind >= 0 ) {
			bag[i].Map( fOutSol );
			for ( k = -1; k < gridPointsA-1; ++k ) {
				newSol[k][ind] = outWork[k];
			}
			outSet[ind] = TRUE;
		}
		else {
			if ( density && strcmp( bag[i].Name(), "density") == 0 ) {
				bag[i].Map( fOutSol );
				for ( k = -1; k < gridPointsA-1; ++k ) {
					density[k] = outWork[k];
				}
				densitySet = TRUE;
			}
		}
	}
	
	for ( i = 0; i < vars; ++i ) {
		if ( outSet[i] == FALSE ) {
			fprintf( fOutFilePtr, "#warning from function 'GetSolution': no match for variable '%s'\n", names[i] );
		}
	}
	if ( density && densitySet == FALSE ) {
		fprintf( fOutFilePtr, "#warning from function 'GetSolution': no match for variable '%s'\n", "density" );
	}
	
#ifdef DEBUGGETSOL
	FILE	*fp = GetOutfile( "NewSol", FileType::kData );
	PrintSolution( fp, gridPointsA-2, vars, &grid[kNext], newSol, fVariableNames );
	fclose( fp );
#endif
	
#ifdef DEBUGINITIALDENS
	if ( density ) {
		FILE	*fpdens = GetOutfile( "DensSol", FileType::kData );
		
		fprintf( fpdens, "*\nZ\tdensity\n" );
		for ( k = -1; k < gridPointsA-1; ++k ) {
			fprintf( fpdens, "%g\t%g\n", grid[k], density[k] );
		}
	
		fclose( fpdens );
	}
#endif

	// clean up
	delete[] outSet;
	DisposeVector( ZGridVec );
	DisposeVector( fOutSol );
}

template<typename Species>
void TTransFlameSolver<Species>::GetSolutionInfo( Double *timeStep )
{
	*timeStep = *fActualTimeStepSize;
}

template<typename Species>
void TTransFlameSolver<Species>::SetOutSolution( void )
{
	int			i, k;
	int			nOfSpecies = fSpecies.GetNOfSpecies();
	int			nSootMoments;
	Double		**theY = fMassFracsWork->mat;
	Double		*theTemp = fTempWork->vec;
	Double		*oTemp = fSolOldTemp->vec;
	Double		*nTemp = fSolTemp->vec;
	Double		**oY = fSolOldMassFracs->mat;
	Double		**nY = fSolMassFracs->mat;
	Double		*oTime = fSolOldTime->vec;
	Double		*nTime = fSolTime->vec;
	Double		currTime = GetCurrentTime();
	Double		**nMom;
	Double		**oMom;
	Double		**theMom;
	Double		oTimek, nTimek, *theYk, *oYk, *nYk;
	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		nMom = fSolSootMoments->mat;
		oMom = fSolOldSootMoments->mat;
		theMom = fSootMomentsWork->mat;
	}

	for ( k = -1; k <= fNGridPoints; ++k )
	{
		oTimek = oTime[k];
		nTimek = nTime[k];
		theYk = theY[k];
		oYk = oY[k];
		nYk = nY[k];
		theTemp[k] = Interpol( currTime, oTemp[k], oTimek, nTemp[k], nTimek );
		for ( i = 0; i < nOfSpecies; ++i ) {
			theYk[i] = Interpol( currTime, oYk[i], oTimek, nYk[i], nTimek );
		}
		if ( fSoot ) {
			for ( i = 0; i < nSootMoments; ++i ) {
				theMom[k][i] = Interpol( currTime, oMom[k][i], oTimek, nMom[k][i], nTimek );
			}
		}
	}
}

template<typename Species>
void TTransFlameSolver<Species>::PrintSolution( FILE *fp, int nGPoints, int nEq, Double *grid
									, Double **sol, ConstStringArray names )
{
	int		i, k;
	char	format[128];
	int 	*len = new int[nEq];
	
	fprintf( fp, "*\n%-12s", "Z" );
	
	for ( i = 0; i < nEq; ++i ) {
		len[i] = maxint( strlen( names[i] ), 12 );
		sprintf( format, "%s%%-%ds", "\t", len[i] );
		fprintf( fp, format, names[i] );
	}
	
	for ( k = -1; k < nGPoints+1; ++k )
	{
		fprintf( fp, "\n%-12E", grid[k] );
		for ( i = 0; i < nEq; ++i ) {
			sprintf( format, "%s%%-%dE", "\t", len[i] );
            fprintf( fp, format, sol[k][i] );
		}
	}
	
	delete[] len;
}

template<typename Species>
void TTransFlameSolver<Species>::MakeGrid( VectorPtr theGrid, Double left, Double right, Flag equidistant )
{
	MakeGrid( theGrid->vec, theGrid->len, left, right, equidistant );
}

template<typename Species>
void TTransFlameSolver<Species>::MakeGrid( Double *grid, int len
							, Double /*left*/, Double /*right*/, Flag equidistant )
{
	Double	deltaZFine;
	
	if ( equidistant ) {
		deltaZFine = 0.0;
	}
	else {
		deltaZFine = 0.2;
	}
	
	int			 	gridPoints = len - 2;
	Double			deltaZCoarse = 1.0 - deltaZFine;
	Double	deltaZ = ( deltaZCoarse + 2.0 * deltaZFine ) 
							/ ( ( Double ) ( gridPoints + 1 ) );
	Double	bound = deltaZFine - 0.5 * deltaZ + 1.0e-10;

	grid[0] = 0.0;
	grid[gridPoints+1] = 1.0;

	for ( int k = 1; k <= gridPoints; ++k ) {
		if ( grid[k-1] <= bound ) {
			grid[k] = grid[k-1] + 0.5 * deltaZ;
		}
		else {
			grid[k] = grid[k-1] + deltaZ;
		}
		if ( grid[k] >= 1.0 ) {
			fprintf( fOutFilePtr, "#error: invalid grid generated: Z[%d] = %g\n", k, grid[k] );
		}
	}
#	ifdef DEBUGGRID
	static int	gridCount = 0;
	char		fName[128];
	sprintf( fName, "grid%2d", gridCount++ );
	FILE		*fp = GetOutfile( fName, FileType::kData );
	
	fprintf( fp, "*\nx\td\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "%g\t1.0\n", grid[k] );
	}
	
	fclose( fp );
#	endif
}


// template<typename Species>                                                // Rui 02/23/2019
// void TTransFlameSolver<Species>::MakeGrid( Double *grid, int len
// 							, Double /*left*/, Double /*right*/, Flag equidistant )
// {
// 	Double	deltaZFine;
	
// 	if ( equidistant ) {
// 		deltaZFine = 0.0;
// 	}
// 	else {
// 		deltaZFine = 0.0012;
// 	}
	
// 	int			 	gridPoints = len - 2;
// 	// Double			deltaZCoarse = 1.0 - deltaZFine;
// 	// Double	deltaZ = ( deltaZCoarse + 2.0 * deltaZFine ) 
// 	// 						/ ( ( Double ) ( gridPoints + 1 ) );
// 	// Double	bound = 0.3;

// 	grid[0] = 0.0;
// 	grid[gridPoints+1] = 1.0;

// 	for ( int k = 1; k <= gridPoints; ++k ) {
// 		if ( grid[k-1] <= 0.3 ) {
// 			grid[k] = grid[k-1] + deltaZFine;
// 		}
// 		else {
// 			grid[k] = grid[k-1] + 1.04*( grid[k-1] - grid[k-2] );
// 		}
// 		if ( grid[k] >= 1.0 ) {
// 			fprintf( fOutFilePtr, "#error: invalid grid generated: Z[%d] = %g\n", k, grid[k] );
// 		}
// 	}

// #	ifdef DEBUGGRID
// 	static int	gridCount = 0;
// 	char		fName[128];
// 	sprintf( fName, "grid%2d", gridCount++ );
// 	FILE		*fp = GetOutfile( fName, FileType::kData );
	
// 	fprintf( fp, "*\nx\td\n" );
// 	for ( k = 0; k < gridPoints+2; ++k ) {
// 		fprintf( fp, "%g\t1.0\n", grid[k] );
// 	}
	
// 	fclose( fp );
// #	endif
// }

template<typename Species>
void TTransFlameSolver<Species>::SetWeights( void )
{
	Double	*grid = fSolGrid->vec;
	Double	h, hm, hnenn;
	
	fgpDens->vec[kPrev] = 1.0 / ( grid[fGrid] - fZl );
	for ( int k = 0; k < fNGridPoints; ++k )
	{
		h = grid[k+1] - grid[k];
		hm = grid[k] - grid[k-1];
		fh[k] = grid[k+1] - grid[k];
		fhm[k] = grid[k] - grid[k-1]; 
		fgpDens->vec[k] = 1.0 / h;
		hnenn = h * hm * ( h + hm );
		fFDWCurr[k] = ( h * h - hm * hm ) / hnenn;
		fFDWMinus[k] = - h * h / hnenn;
		fFDWPlus[k] = hm * hm / hnenn;
		fWCurr[k] = - 2.0 * ( h + hm ) / hnenn;
		fWMinus[k] = 2.0 * h / hnenn;
		fWPlus[k] = 2.0 * hm / hnenn;
	}
}

template<typename Species>
void TTransFlameSolver<Species>::SetFDWeights( Double *y )
{
	int		gpOff;
	double  fhnenn;

	fgpDens->vec[kPrev] = 1.0 / ( y[fGrid] - y[fGrid-fNOfEquations] );
	for ( int k = 0; k < fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		fh[k] = y[gpOff+fGrid+fNOfEquations] - y[gpOff+fGrid];
		fhm[k] = y[gpOff+fGrid] - y[gpOff+fGrid-fNOfEquations];
		fh[kPrev] = fhm[kCurr];
		fhm[fNGridPoints] = fh[fNGridPoints-1];
		fhnenn = fh[k] * fhm[k] * ( fh[k] + fhm[k] );
		fFDWCurr[k] = ( fh[k] * fh[k] - fhm[k] * fhm[k] ) / fhnenn;
		fFDWMinus[k] = - fh[k] * fh[k] / fhnenn;
		fFDWPlus[k] = fhm[k] * fhm[k] / fhnenn;
		fWCurr[k] = - 2.0 * ( fh[k] + fhm[k] ) / fhnenn;
		fWMinus[k] = 2.0 * fh[k] / fhnenn;
		fWPlus[k] = 2.0 * fhm[k] / fhnenn;	
		fgpDens->vec[k] = 1.0 / fh[k];
	}
}

template<typename Species>
void TTransFlameSolver<Species>::ComputeGPDistributionWeights( Double *y )
{
	int			nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double		dTdZ, dYdZ; 
	Double		*maxVals = fMaxVals->vec;
	int			j, k;
	int			gpOff, sootOff, nSootMoments;
	if ( fSoot ) {
		sootOff = fSoot->GetOffsetSootMoments();
		nSootMoments = fSoot->GetNSootMoments();
	}
	
	fMonFct[kPrev] = fMonFct[fNGridPoints] = 1.0;
	for ( k = 0; k < fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		fMonFct[k] = 1.0;

		dTdZ = ( y[gpOff+fTemperature+fNOfEquations] - y[gpOff+fTemperature] ) / fh[k];
		fMonFct[k] += dTdZ * dTdZ / MAX( 1.0e-10, maxVals[fTemperature] * maxVals[fTemperature] );

		for ( j = 0; j < nSpeciesIn; ++j ) {
			dYdZ = ( y[gpOff+fNOfEquations+fFirstSpecies+j] - y[gpOff+fFirstSpecies+j] ) / fh[k];
			fMonFct[k] += dYdZ * dYdZ / MAX( 1.0e-10, maxVals[fFirstSpecies+j] * maxVals[fFirstSpecies+j] ); //Rui 12/31/18
		}

		if ( fSoot ) {
			for ( j = 0; j < nSootMoments; ++j ) {
				dYdZ = ( y[gpOff+fNOfEquations+sootOff+j] - y[gpOff+sootOff+j] ) / fh[k];
					fMonFct[k] += dYdZ * dYdZ / MAX( 1.0e-10, maxVals[sootOff+j] * maxVals[sootOff+j] );
			}
		}
		fMonFct[k] = sqrt( fMonFct[k] );
		
	}

	
	fMonFct[fNGridPoints] = fMonFct[fNGridPoints-1];
}

template<typename Species>
int TTransFlameSolver<Species>::GetActualPoint( Double tEnd )
{
	int		minPoint = 0;
	Double	*t = fSolTime->vec;
	Double	minValue = t[0];
	
	for ( int k = 1; k < fNGridPoints; ++k )
	{
		if ( t[k] < minValue ) {
			minPoint = k;
			minValue = t[k];
		}
	}
	
	if ( minValue < tEnd ) {
		return minPoint;
	}
	else {
		return -1;
	}
}

template<typename Species>
void TTransFlameSolver<Species>::SaveSolution( int k, Double t, Double *y )
{
	int		nOfSpecies = fSpecies.GetNOfSpecies();
	int		nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int		nSootMoments, momOff;
	Double	*YNew = fSolMassFracs->mat[k];
	Double	*YOld = fSolOldMassFracs->mat[k];
	Double	*momNew;
	Double	*momOld;
	
	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		momOff = fSoot->GetOffsetSootMoments();
		momNew = fSolSootMoments->mat[k];
		momOld = fSolOldSootMoments->mat[k];
	}
	
	// save old solution
	fSolOldTime->vec[k] = fSolTime->vec[k];
	copy( nOfSpecies, YNew, 1, YOld, 1 );
	if ( fSoot ) {
		copy( nSootMoments, momNew, 1, momOld, 1 );
	}
	fSolOldTemp->vec[k] = fSolTemp->vec[k];

	// save new solution
	fSolTime->vec[k] = t;
	copy( nOfSpeciesIn, &y[fFirstSpecies], 1, YNew, 1 );

	if ( fSoot ) {
		Double	mixMolarMass;
		fProperties->ComputeMixtureMolarMass( mixMolarMass, YNew, fSpecies.GetMolarMass()->vec, nOfSpeciesIn );
		
		for ( int i = 0; i < nSootMoments; ++i ) {
			momNew[i] = y[momOff+i] * GetPressure( t ) * mixMolarMass 
									/ ( RGAS * y[fTemperature] );
		}
	}
	fSolTemp->vec[k] = y[fTemperature];
	fSolGrid->vec[k] = y[fGrid];
}

template<typename Species>
void TTransFlameSolver<Species>::SetMaxVals( void )
{
	int		nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int		nSootMoments, momOff;
	Double	**Y = fSolMassFracs->mat;
	Double	**mom;
	Double	*maxvals = fMaxVals->vec;
	
	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		momOff = fSoot->GetOffsetSootMoments();
		mom = fSolSootMoments->mat;
	}
	
	// save new solution
	for ( int i = 0; i < nOfSpeciesIn; ++i ) {
		maxvals[fFirstSpecies+i] = Y[LocationOfMax( fNGridPoints+2, &Y[kPrev][i]
									, fSolMassFracs->phys_rows )-1][i];
	}
	
	if ( fSoot ) {		
		int loc;
		for ( int i = 0; i < nSootMoments; ++i ) {
			loc = LocationOfMax( fNGridPoints+2, &mom[kPrev][i]
									, fSolSootMoments->phys_rows )-1;
			maxvals[momOff+i] = mom[loc][i] / fDensity->vec[loc];
		}
	}

	maxvals[fTemperature] = fSolTemp->vec[LocationOfMax( fNGridPoints+2, &fSolTemp->vec[kPrev]
							, 1 )-1];
}

/* template<typename Species> */
/* Flag TTransFlameSolver<Species>::FirstStep( void ) */
/* { */
/* 	Flag	leave; */
/* 	 */
/* 	for ( int k = 0; k < fNGridPoints; ++k ) { */
/* 		fActualPoint = k; */
/* 		if ( leave = OneStep( k ) ) { */
/* 			break; */
/* 		} */
/* 		fInfo[k]->vec[kF11] = 0; */
/* 	} */
/*  */
/* 	fFirstCall = FALSE; */
/* 	return leave; */
/* } */

/*  */
/* template<typename Species> */
/* Flag TTransFlameSolver<Species>::OneStep( int k ) */
/* { */
/* 	SetMaxTimeStep( k, fTime->vec[k] ); */
/* 	DDASSL( ::ResTransFlameSolver<TTransFlameSolver<Species> >, &fDasslNEq, &fTime->vec[k], fSolution->mat[k] */
/* 			, fSolPrime->mat[k], &fTEnd, fInfo[k]->vec, fRTol, fATol, &fIdid */
/* 			, fRWork[k]->vec, &fLRW, fIWork[k]->vec, &fLIW, NULL, ( int * ) this */
/* 			, ::JacTransFlameSolver<TTransFlameSolver<Species> > ); */
/*  */
/* 	fMaxStepTaken[k] = MAX( fMaxStepTaken[k], fRWork[k]->vec[kF7] ); */
/* 	if ( *fNActualStep[k] % DELTAPROG == 0 ) { */
/* 		fprintf( fOutFilePtr, "gp = %3d    max stp = %8.2e    stp = %5d    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n" */
/* 			, k, fMaxStepTaken[k], *fNActualStep[k], *fNActualOrd[k], fTime->vec[k] * 1000.0, fSolution->mat[k][fTemperature] ); */
/* 		fflush( fOutFilePtr ); */
/* 	} */
/*  */
/* 	if ( fIdid < 1 ) { */
/* 		fprintf( fOutFilePtr, "#error: ddassl error no. %d occured at gridpoint no. %d\n" */
/* 			, fIdid, k ); */
/* 		return TRUE; */
/* 	} */
/* 	else { */
/* 		if ( fRWork[k]->vec[kF3] > fTime->vec[k] ) { */
/* 			fprintf( fOutFilePtr, "#warning: integration carried out beyond tEnd\n" ); */
/* 			fprintf( fOutFilePtr, "info[4] = %d\n", fInfo[k]->vec[kF4] ); */
/* 			fprintf( fOutFilePtr, "fIdid = %d\n", fIdid ); */
/* 		} */
/* 		SaveSolution( k, fTime->vec[k], fSolution->mat[k] ); */
/* 		SetMaxVals(); */
/* 		return FALSE; */
/* 	} */
/* } */

template<typename Species>
Flag TTransFlameSolver<Species>::OneImpliStep( void )
{
  static int   initrandomfile=0;
  static Double timrandomchange=0.0;
  const Double deltarandomchange = 1.0e-4;
  static FILE *chiout;
  
    if (initrandomfile==0) {
        ++initrandomfile;
        chiout=GetOutfile( "RandomChi", FileType::kData );
        fprintf(chiout,"*\nChi\trand\tT\n");
	}



#ifdef CVODE
	int flag;
	long int nst;         // current time step number

    flag = CVode(fMem, fTEnd,fCVY, &fTime->vec[0], CV_ONE_STEP);  // ResTransFlameImpliSolver will be called here!!!!!!!!!!!!

	SolutionFromCVode();

    CVodeGetNumSteps(fMem, &nst);
    fNActualStepCV = nst;
    CVodeGetLastStep(fMem, &fActualTimeStepSizeCV);
    CVodeGetLastOrder(fMem, &fNActualOrdCV);

	if ( *fNActualStep[0] % DELTAPROG == 0 ) {
		fprintf( fOutFilePtr, "stpCVode = %5d    dt = %8.2e s    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n"
			, *fNActualStep[0], *fActualTimeStepSize, *fNActualOrd[0], fTime->vec[0] * 1000.0
			, fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
		fflush( fOutFilePtr );
	}

	if (flag < 0) {
		cerr << "#error: cvode error no. " << flag << " occured" << NEWL;
		return TRUE;
	} 
#else
	DDASSL( ::ResTransFlameImpliSolver<TTransFlameSolver<Species> >, &fDasslNEq, &fTime->vec[0], fSolution->mat[kPrev]
			, fSolPrime->mat[kPrev], &fTEnd, fInfo[0]->vec, fRTol, fATol, &fIdid
			, fRWork[0]->vec, &fLRW, fIWork[0]->vec, &fLIW, NULL, ( int * ) this
			, NULL );
	if ( *fNActualStep[0] % DELTAPROG == 0 ) {
		fprintf( fOutFilePtr, "stp = %5d    dt = %8.2e s    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n"
			, *fNActualStep[0], *fActualTimeStepSize, *fNActualOrd[0], fTime->vec[0] * 1000.0
			, fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
		fflush( fOutFilePtr );
	}

	if ( fIdid < 1 ) {
		fprintf( fOutFilePtr, "#error: dassl error no. %d occured at gridpoint no. %d\n"
			, fIdid, *fNActualStep[0] );
		fprintf( fOutFilePtr, "number of error test failures is %d\n", fIWork[0]->vec[kF14] );
		fprintf( fOutFilePtr, "number of convergence test failures is %d\n", fIWork[0]->vec[kF15] );
		fprintf( fOutFilePtr, "dassl suggested order %d and stepsize %g for the next step\n"
			,fIWork[0]->vec[kF7], fRWork[0]->vec[kF3] );
		return TRUE;
	}
	if ( fRWork[0]->vec[kF3] > fTime->vec[0] ) {
		fprintf( fOutFilePtr, "#warning: integration carried out beyond tEnd\n" );
		fprintf( fOutFilePtr, "info[4] = %d\n", fInfo[0]->vec[kF4] );
		fprintf( fOutFilePtr, "fIdid = %d\n", fIdid );
	}
#endif

	for ( int k = -1; k <= fNGridPoints; ++k ) {
		SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
	}

	SetMaxVals();
		
	if ( fSoot ) {
		int maxloc0 = 0;
		int maxloc1 = 0;
		for ( int k = 1; k < fNGridPoints; ++k ) {
			if ( fSolSootMoments->mat[k][0] > fSolSootMoments->mat[maxloc0][0] ) {
				maxloc0 = k;
			}
			if ( fSolSootMoments->mat[k][1] > fSolSootMoments->mat[maxloc1][1] ) {
				maxloc1 = k;
			}
		}
		fprintf( fOutFilePtr, "M0_max at Z=%g is %g\tM1_max at Z=%g is %g\tfv_max = %g\n\n"
				, fSolGrid->vec[maxloc0], fSolSootMoments->mat[maxloc0][0]
				, fSolGrid->vec[maxloc1], fSolSootMoments->mat[maxloc1][1]
				, fSolSootMoments->mat[maxloc1][1] * 24.0 / 1800.0 );
	}
		
#ifdef ENTHALPYLIB 
	Double	tNow = fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature];
	if ( CheckOutput( fTime->vec ) ) {
		cerr << "dump libout" << NEWL; 
		char   tl[128]; 
		sprintf( tl, "Lib_Chi%05g_T%04.0f_t%05.2g", GetRefDissRate( fTime->vec[0] ), tNow, fTime->vec[0]*1000 ); 
		FILE *fp = GetOutfile( tl, FileType::kText );
		WriteFlameletFile( fp, NULL, NULL ); 
		fclose( fp );
	}
#endif

#ifdef LOGNORMALCHI
    if ( fTime->vec[0] > timrandomchange+deltarandomchange ) {
		fprintf( chiout, "%g\t%g\t%g\n", GetRefDissRate( fTime->vec[0] ), fRandomNumber, fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
        fflush( chiout );
		SetRandomNumber();
        timrandomchange=fTime->vec[0];
		fprintf( fOutFilePtr, "new random number: %g\n", fRandomNumber );
	}
#endif

	return FALSE;
}

template<typename Species>
Flag TTransFlameSolver<Species>::CheckOutput( Double *t )
{
	const int		vars = 3;
	static Flag		init = FALSE;
	int				i,k;
	Flag			ord = FALSE;
	Double			max[vars];
	Double			pressure = GetPressure( *t );
	static Double	**valsnew, **valsold;
	Double			**nY = fSolMassFracs->mat;
	Double			*nTemp = fSolTemp->vec;
	Double			*prodRate = fSpecies.GetProductionRate()->vec;

	if ( !init ) {
		valsold = New2DArray( fNGridPoints, vars );
		valsnew = New2DArray( fNGridPoints, vars );
		init = TRUE;
	}

// compute W and C
	int	CO2 = fSpecies.FindSpecies( "CO2" );
	int	CO = fSpecies.FindSpecies( "CO" );
	int	H2O = fSpecies.FindSpecies( "H2O" );
	int	H2 = fSpecies.FindSpecies( "H2" );

	for ( int k = 0; k < fNGridPoints; ++k ) {

#ifdef DELTAINEW
		UpdateThermoProps( k, nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? fSolSootMoments->mat[k] : NULL );
#else
		T0DFlame<Species>::UpdateThermoProps( nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? fSolSootMoments->mat[k] : NULL );
#endif

		valsnew[k][0] = fSolution->mat[k][fTemperature];
		valsnew[k][1] = fSolution->mat[k][fFirstSpecies+CO2] + fSolution->mat[k][fFirstSpecies+CO] 
							+ fSolution->mat[k][fFirstSpecies+H2O] + fSolution->mat[k][fFirstSpecies+H2];
		valsnew[k][2] = prodRate[CO2] + prodRate[CO] + prodRate[H2O] + prodRate[H2];
	}

//	get maxes
	max[0] = fabs( fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
	max[1] = fabs( valsnew[LocationOfMax( fNGridPoints, &valsnew[0][1], vars )][1] );
	max[2] = fabs( valsnew[LocationOfMax( fNGridPoints, &valsnew[0][2], vars )][2] );
		
	for ( k = 0; k < fNGridPoints; ++k ) {
		for ( i = 0; i < vars; ++i ) {
			if ( fabs( ( valsnew[k][i] - valsold[k][i] ) / MAX( 1.0e-10, max[i] ) ) > 0.2 ) {
				fprintf( stderr, "valsnew[k][i] = %g  valsold[k][i] = %g, i = %d, k = %d, max = %g\n", valsnew[k][i], valsold[k][i], i, k, max[i]);
				ord = TRUE;
				break;
			}
		}
		if ( ord == TRUE ) break;
		if ( !ord && fabs( valsnew[k][0] - valsold[k][0] ) > 20.0 ) {
			fprintf( stderr, "Temp 20K up to %g at k = %d\n", valsnew[k][0], k );
			ord = TRUE;
			break;
		}
		if ( fabs( valsnew[k][1] - valsold[k][1] ) > 0.01 ) {
			fprintf( stderr, "ProgVar 0.01 up to %g at k = %d\n", valsnew[k][1], k );
			ord = TRUE;
			break;
		}

	}
	if ( ord == TRUE ) {
		for ( k = 0; k < fNGridPoints; ++k ) {
			for ( i = 0; i < vars; ++i ) {
				valsold[k][i] = valsnew[k][i];
			}
		}
	}
	
	return ord;
}

template<typename Species>
void TTransFlameSolver<Species>::PostIter( Double t )
{
	int		nGm1 = fNGridPoints-1;
	Double	**fSol = fSolution->mat;
	Double	slope = *fActualTimeStepSize * fdDeltaZdt / ( fSolGrid->vec[fNGridPoints] - fSolGrid->vec[nGm1] );
	
	fSolOldTime->vec[fNGridPoints] = fSolTime->vec[fNGridPoints];
	fSolTime->vec[fNGridPoints] = t;

	fSolOldTime->vec[kPrev] = fSolTime->vec[kPrev];
	fSolTime->vec[kPrev] = t;

#ifdef MOVEZRIGHT
	for ( int i = 0; i < fNOfEquations; ++i ) {
		fSol[i][fNGridPoints] = ( fSol[i][fNGridPoints] - slope * fSol[i][nGm1] )
								/ ( 1.0 - slope );
	}
#else
	fSol[fTemperature][fNGridPoints] = Interpol( t, fTempFuelStart, fTStart
									, fTempFuelEnd, fTEnd );
#endif
}

template<typename Flame>
void ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar )
{
	( ( Flame* ) iPar )->ResTransFlameSolver( T, y, yPrime, delta
			, iRes, rPar, iPar );
}

template<typename Species>
void TTransFlameSolver<Species>::ResTransFlameSolver( Double *t, Double *y, Double *yPrime, Double *delta
			, int * /*iRes*/, Double * /*rPar*/, int * /*iPar*/ )
{
	int						nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int						i, ieq, kAct = fActualPoint;
	Double					sum = 0.0;
	Double					tempMinus;
	Double					tempPlus;


	Double					*prodRate = fSpecies.GetProductionRate()->vec;
	Double					*enth = fSpecies.GetEnthalpy()->vec;
	Double					*molarMass = fSpecies.GetMolarMass()->vec;
	Double					*Le = fSpecies.GetLewisNumber()->vec;

	Double					temp = y[fTemperature];
	Double					*YF = &y[fFirstSpecies];

	Double					*oT = &fSolOldTemp->vec[kAct];
	Double					*oTime = &fSolOldTime->vec[kAct];
	Double					*nT = &fSolTemp->vec[kAct];
	Double					*nTime = &fSolTime->vec[kAct];
	Double					**nY = &fSolMassFracs->mat[kAct];

	Double					**nMom;
	Double					*moments;
	int						sootOff;
	int						nSootMoments;
	if ( fSoot ) {
		nMom             =  &fSolSootMoments->mat[kAct];
		moments          =  fSoot->GetMoments()->vec;
		nSootMoments     =  fSoot->GetNSootMoments();
		sootOff          =  fSoot->GetOffsetSootMoments();
	}

	// copy nY to oY and then YF to nY
	copy( nSpeciesIn, YF, 1, nY[kCurr], 1 );

	Double					Z = fSolGrid->vec[kAct];
	Double					chi = GetDissRate( t[kCurr], Z );
	Double					pressure = GetPressure( t[kCurr] );

	
	if ( fSoot ) {
		fSoot->MOverRhoToM( &y[sootOff], moments, nSootMoments
				, nY[kCurr], temp, pressure
				, molarMass, nSpeciesIn, fProperties );
	}

	T0DFlame<Species>::UpdateThermoProps( nY[kCurr], temp, pressure, fProperties->GetDensityRef()
									, kDensFromPress, moments );

	Double				rho = fProperties->GetDensity();
	Double				heatCap = fProperties->GetMixHeatCapacity();
	Double				mixMolarMass = fProperties->GetMixMolarMass();

	for ( i = 0; i < nSpeciesIn; ++i ) {
		ieq = fFirstSpecies + i;
		delta[ieq] = yPrime[ieq] - prodRate[i] / rho - chi / ( 2.0 * Le[i] ) 
			* ( fWMinus[kAct] * nY[kPrev][i]
			+ fWCurr[kAct] * YF[i]
			+ fWPlus[kAct] * nY[kNext][i] );
		sum += enth[i] * prodRate[i];
	}
	
	tempMinus = nT[kPrev];
	tempPlus = nT[kNext];

	if ( kAct == 0 ) {
		tempMinus = Interpol( t[kCurr], oT[kPrev], oTime[kPrev], nT[kPrev], nTime[kPrev] );
	}
	if ( kAct == fNGridPoints-1 ) {
		tempPlus = Interpol( t[kCurr], oT[kNext], oTime[kNext], nT[kNext], nTime[kNext] );
	}

	delta[fTemperature] = yPrime[fTemperature] 
			+ ( sum - fDPdt ) / ( rho * heatCap ) - 0.5 * chi
			* ( fWMinus[kAct] * tempMinus
			+ fWCurr[kAct] * temp
			+ fWPlus[kAct] * tempPlus )
			- GetTempSource( Z * GetZR() );
	if ( fProperties->GetRadiation() ) { 
		delta[fTemperature] -= ( fProperties->GetRadiation()->GetRadiation() ) / ( rho * heatCap );
	}
	
	if ( fSoot ) {
		// ATTENTION (hp)
		// calculation of soot with semi explicit solver currently doesn't
		// work, since y[sootOff] is M/rho and nMom is M
		fprintf( fOutFilePtr, "###error: calculation of soot with semi explicit solver currently doesn't work\n" );
		exit(2);

		Double	fracIndex;
		for ( i = 0; i < nSootMoments; ++i ) {
			fracIndex = i - 2.0 / 3.0;
			ieq = sootOff + i;
			delta[ieq] = yPrime[ieq] 
				- fSoot->NucleationNew( i, temp, fSoot->GetPAHMoments()->vec ) / rho 
				- chi / ( 2.0 * fSoot->GetLewis1() ) 
#ifdef SIZEDEPDIFFUSION
				* ( fWMinus[kAct] * fSoot->FracMom2( fracIndex, nMom[kPrev] )
				+ fWCurr[kAct] * fSoot->FracMom2( fracIndex, &y[sootOff] )
				+ fWPlus[kAct] * fSoot->FracMom2( fracIndex, nMom[kNext] ) );
#else
				* ( fWMinus[kAct] * nMom[kPrev][i]
				+ fWCurr[kAct] * y[sootOff+i]
				+ fWPlus[kAct] * nMom[kNext][i] );
#endif
		}
	}
#ifdef DEBUGRES
	for ( int j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
		fprintf( stderr, "%s\t%g\t%g\n", fVariableNames[j], y[j], delta[j] );
	}
	fprintf( stderr, "\n" );
	for ( j = nSpeciesIn; j < fSpecies.GetNOfSpecies(); ++j ) {
		fprintf( stderr, "SteadyState %s:\t%g\n", fSpecies.GetNames()[j], nY[kCurr][j] );
	}
	fprintf( stderr, "\n" );
#endif
	
}

template<typename Flame>
void ResTransFlameImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar )
{
	( ( Flame* ) iPar )->ResTransFlameImpliSolver( T, y, yPrime, delta
			, iRes, rPar, iPar );
}

template<typename Species>
void TTransFlameSolver<Species>::ResTransFlameImpliSolver( Double *t, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double * /*rPar*/, int * /*iPar*/ )
{
	int						k;
	int						nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int						i, ieq, gpOff, gpSpecOff;
	int						gpGridOff;
	Double					*gpDens = fgpDens->vec;
    Double					kappa_ = fKappa * ( 1.0 + fKappa );
	Double					sum, sumYH;
	Double					*prodRate = fSpecies.GetProductionRate()->vec;
	Double					*enth = fSpecies.GetEnthalpy()->vec;
	Double					*molarMass = fSpecies.GetMolarMass()->vec;
	string                  GasRadiationName = fProperties->GetRadiationName();

#ifdef LEWISCHANGE                                                    
	Double					*LeOrig = fSpecies.GetLewisNumber()->vec;    
	Double					*Le = new Double[nSpeciesIn];
	Double					LeFunc;
	Double					switchTime = LEWISSWITCHTIME;
	Double					switchPeriod = LEWISSWITCHPERIOD;
	const Double			Pi = 4.0 * atan( 1.0 );
	Double					omega = 2.0 * Pi / switchPeriod;
	Double					tStar = *t - ( switchTime - 0.5 * switchPeriod );


	if ( tStar >= 0.0 ) {
		if ( tStar < switchPeriod ) {
			LeFunc = 0.5 * ( cos( omega * tStar ) + 1.0 );
			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
				Le[iLe] = LeFunc * ( LeOrig[iLe] - 1.0 ) + 1.0;
			}
		}
		else {
			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
				Le[iLe] = 1.0;
			}
		}
	}
	else {
		for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
			Le[iLe] = LeOrig[iLe];
		}
	}
#else
	Double					*Le = fSpecies.GetLewisNumber()->vec;    
#endif

	Double					**nY = fSolMassFracs->mat;
	Double					*nTemp = fSolTemp->vec;
	Double					*Z = fSolGrid->vec;
	Double					chi, chiPrev, chiNext, chi_st, MassDiffu;
	Double					rhoChiCurr, rhoChiPrev, rhoChiNext, rhoChiCpOverLambda;
	Double					pressure = GetPressure( *t );
	Double					sumCpidYidZ;
	Double	                *cpiOld = fSpecies.GetHeatCapacity()->vec;
	Double					*cpMix = fHeatCpMix->vec;
	Double					*mu = fViscosity->vec;
	Double					*mixMolarMass = fMolarMassMix->vec;
	Double					*lambdaOverCp = fLambdaOverCpMix->vec;
	Double					**cpi = fSpecHeatCp->mat;
	Double					**hi = fSpecEnthalpy->mat;
	Double					**mDoti = fProdRate->mat;
	Double					*rho = fDensity->vec;
	Double					*diffTermY = fDiffTermY->vec;
	Double					diffCorrY;
	Double					sumdYdZ, sumdYdZOverLe;


#ifdef MOLARDIFF
	Double					diffCorrW;
	Double					*diffTermW = fDiffTermW->vec;
	Double					d2WdZ2, dWdZ;
	Double					diffWFact;
#endif

	Double					halfChi;
	Double					**moments;
	Double					*MOverRhoPrev;
	Double					*MOverRhoCurr;
	Double					*MOverRhoNext;
	int						sootOff;
	int						nSootMoments;

	y = &y[fNOfEquations];                // y is everything, including temp, mass species, grid, 
	yPrime = &yPrime[fNOfEquations];      // zeros here
	delta = &delta[fNOfEquations];        // Enthalpy defects

#ifdef READU			
	int		iU, kU, ieqU;
	Double	UstOverU;
	for ( kU = -1; kU <= fNGridPoints; ++kU ) {
		ieqU = kU * fNOfEquations;
		UstOverU = GetU( *t, Z[kU] );
		for ( iU = 0; iU < fNOfEquations; ++iU ) {
			yPrime[ieqU+iU] /= UstOverU;
		}
	}
#endif

   	    SetFDWeights( y );
		ComputeGPDistributionWeights( y );

	if ( fSoot ) {
		moments = fSolSootMoments->mat;
		nSootMoments = fSoot->GetNSootMoments();
		sootOff = fSoot->GetOffsetSootMoments();
	}
	

	// copy YF to nY
	for ( k = -1; k <= fNGridPoints; ++k ) {
		
		copy( nSpeciesIn, &y[k*fNOfEquations+fFirstSpecies], 1, nY[k], 1 );
		nTemp[k] = y[k*fNOfEquations+fTemperature];                       
		
		if ( fSoot ) {
			if ( k > -1 ) MOverRhoPrev = &y[(k-1)*fNOfEquations+sootOff];
			if ( k < fNGridPoints ) MOverRhoNext = &y[(k+1)*fNOfEquations+sootOff];
			MOverRhoCurr = &y[k*fNOfEquations+sootOff];
			fSoot->MOverRhoToM( &y[k*fNOfEquations+sootOff], moments[k], nSootMoments
					, nY[k], nTemp[k], pressure
					, molarMass, nSpeciesIn, fProperties );
		}

		if ( nTemp[k] < 100.0 || nTemp[k] > 5000.0 ) {
			fprintf( fOutFilePtr, "###Warning: Temp = %g @ gp = %d and Z = %g\n", nTemp[k], k, Z[k] );
			*iRes = -1;
			return;
		}

#ifdef DELTAINEW
		UpdateThermoProps( k, nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? moments[k] : NULL );
#else
		T0DFlame<Species>::UpdateThermoProps( nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? moments[k] : NULL );
#endif

		cpMix[k] = fProperties->GetMixHeatCapacity();
		rho[k] = fProperties->GetDensity();
		mu[k] = fProperties->GetMixViscosity();
		lambdaOverCp[k] = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();
		mixMolarMass[k] = fProperties->GetMixMolarMass();
		copy( nSpeciesIn, fSpecies.GetHeatCapacity()->vec, 1, cpi[k], 1 );
		copy( nSpeciesIn, fSpecies.GetEnthalpy()->vec, 1, hi[k], 1 );
		copy( nSpeciesIn, fSpecies.GetProductionRate()->vec, 1, mDoti[k], 1 );
	}
	
/****************************************************************************************************************
1. Finished updating nTemp[k], nY[k] ... using y calculated from the last step.

2. Now, we should have everything we need to calculate the radiation source term of WSGG/Grey/SNB at each gird point.
******************************************************************************************************************/
	
	if ( GasRadiationName == "WSGG" ||  GasRadiationName == "Grey" ||  GasRadiationName == "SNB" || GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar") {

		vector<double> GasRadSourceTerm(fNGridPoints+2);

	// 1. Compute physical grid from mixture fraction space
		vector<double> physicalGrid(fNGridPoints+2);
		double DeltaX = 0.01;
		vector<double> heatrel(fNGridPoints+2);
		for (int k = 0; k < fNGridPoints+2; ++k ) {
			heatrel[k] = 0.0;
			for ( i = 0; i < nSpeciesIn; ++i ) {
				heatrel[k] += hi[k-1][i] * mDoti[k-1][i];
			}
			heatrel[k] = - heatrel[k];
		}
		double chi_ref = GetRefDissRate(*t);
		physicalGrid[0] = DeltaX;
		for (int k = 1; k < fNGridPoints+1; ++k ) {
			chi = GetDissRate( *t, Z[k-1] ) ;
			MassDiffu = lambdaOverCp[k-1] / rho[k-1] ;
			physicalGrid[k] = physicalGrid[k-1] + sqrt( 2 * MassDiffu / chi ) * (Z[k-1] - Z[k-2]) ; 
		}
		physicalGrid[fNGridPoints+1] = physicalGrid[fNGridPoints] + DeltaX;

	// 2. Compute locTemp, locMoleMass of H2O and CO2
		int fH2OIndex = fInputData->fH2OIndex;
		int fCO2Index = fInputData->fCO2Index;
    	vector <double> locMoleMassH2O (fNGridPoints+2, 0.);
    	vector <double> locMoleMassCO2 (fNGridPoints+2, 0.);
    	vector <double> 	   locTemp (fNGridPoints+2, 0.);
    	for ( int k = 0; k < fNGridPoints+2; ++k) { 
        	locMoleMassH2O[k] = nY[k-1][fH2OIndex] * mixMolarMass[k-1] / molarMass[fH2OIndex];
        	locMoleMassCO2[k] = nY[k-1][fCO2Index] * mixMolarMass[k-1] / molarMass[fCO2Index];
        	locTemp[k]        = nTemp[k-1];
    	}
		
		fProperties->GetRadiation()->ComputeRadiationUnsteady(locTemp, locMoleMassH2O, locMoleMassCO2, nSpeciesIn, physicalGrid, 
    																    DeltaX, chi_ref, heatrel, GasRadiationName);
	}

/****************************************************************************************************************
1. The next for loop is to update y due to combustion and / or radiation. 

2. nTemp[k], nY[k] ... will not be changed for this loop. 

3. Radiation source term can be called depending on gas radiation models.
******************************************************************************************************************/

// **** to be changed ****
// catch first and last point
// *** end of change ***

// **** to be changed ****
#ifdef DEBUGSOOTSOURCE
	fprintf( fOutFilePtr, "\n*\nZ\tNuc\tCond\tSG\tOx\tDiff\n" );
#endif

	// double chi_ref = GetRefDissRate(*t);    // By Rui

	for ( k = 0; k < fNGridPoints; ++k ) {
// *** end of change ***
		gpOff = k * fNOfEquations;
    	gpGridOff = gpOff + fGrid;	
		
		chi = GetDissRate( *t, Z[k] );
		chiPrev = GetDissRate( *t, Z[k-1] );
		chiNext = GetDissRate( *t, Z[k+1] );

		// chi = chi_ref;             // By Rui
		// chiPrev = chi_ref;
		// chiNext = chi_ref;

		// if ( chi_ref < 0.01){     // By Rui
		// 	chi = chi_ref;
		// 	chiPrev = chi_ref;
		// 	chiNext = chi_ref;
		// }

		rhoChiCurr = rho[k] * chi;
		rhoChiPrev = rho[k-1] * chiPrev;
		rhoChiNext = rho[k+1] * chiNext;
		rhoChiCpOverLambda = rhoChiCurr / lambdaOverCp[k];
		halfChi = 0.5 * chi;


		
#ifdef MOLARDIFF
		d2WdZ2 = fWMinus[k] * mixMolarMass[k-1]
				+ fWCurr[k] * mixMolarMass[k]
				+ fWPlus[k] * mixMolarMass[k+1];
		dWdZ = fFDWMinus[k] * mixMolarMass[k-1]
				+ fFDWCurr[k] * mixMolarMass[k]
				+ fFDWPlus[k] * mixMolarMass[k+1];
		diffWFact = halfChi / mixMolarMass[k] * d2WdZ2;
		diffCorrW = 0.0;
#endif

		sumdYdZ = sumdYdZOverLe = 0.0;
		diffCorrY = 0.0;
		Double	dYdZ;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			diffTermY[i] = halfChi / Le[i]
				* ( fWMinus[k] * nY[k-1][i]
				+ fWCurr[k] * nY[k][i]
				+ fWPlus[k] * nY[k+1][i] );
			diffCorrY += diffTermY[i];
			dYdZ = fFDWMinus[k] * nY[k-1][i]
				+ fFDWCurr[k] * nY[k][i]
				+ fFDWPlus[k] * nY[k+1][i];
			sumdYdZ += dYdZ;
			sumdYdZOverLe += dYdZ / Le[i];
#ifdef MOLARDIFF
			diffTermW[i] = diffWFact * nY[k][i] / Le[i];
			diffCorrW += diffTermW[i];
#endif
		}
				
#ifdef NEWPROPERTIES
		// If is true but void. Do nothing here and below
#else

#ifdef DELTAINEW
		UpdateThermoProps( k, nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? moments[k] : NULL );
#else 
		T0DFlame<Species>::UpdateThermoProps( nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? moments[k] : NULL );
#endif
		copy( nSpeciesIn, fSpecies.GetHeatCapacity()->vec, 1, cpi[k], 1 );
		copy( nSpeciesIn, fSpecies.GetEnthalpy()->vec, 1, hi[k], 1 );
		copy( nSpeciesIn, fSpecies.GetProductionRate()->vec, 1, mDoti[k], 1 );
#endif


		gpSpecOff = gpOff + fFirstSpecies;

		Double	constConvVeloY = 0.25 / ( rho[k] * lambdaOverCp[k] )
						* ( fFDWMinus[k] * lambdaOverCp[k-1] * rhoChiPrev 
							+ fFDWCurr[k] * lambdaOverCp[k] * rhoChiCurr
							+ fFDWPlus[k] * lambdaOverCp[k+1] * rhoChiNext );
		sum = 0.0;
		sumYH = 0.0;
		sumCpidYidZ = 0.0;

// Species Mass Fraction Equations
		
		for ( i = 0; i < nSpeciesIn; ++i ) {
			ieq = gpSpecOff + i;
			delta[ieq] = yPrime[ieq] - mDoti[k][i] / rho[k] 
					- diffTermY[i]
#ifdef MOLARDIFF
					- diffTermW[i] 
#endif
#ifdef DIFFCORR
					+ nY[k][i] * ( diffCorrY 
#	ifdef MOLARDIFF
								+ diffCorrW
#	endif
								)
#endif
								;
#ifdef CONVECTION
		// 5. mass fraction convection
			Double	convVeloY = ( 1.0 / Le[i] - 1.0 ) * constConvVeloY;
#	ifdef UPWIND
			if ( convVeloY < 0.0 ) {
				delta[ieq] -= convVeloY * ( nY[k][i] - nY[k-1][i] ) / ( fhm[k] );
			}
			else {
				delta[ieq] -= convVeloY * ( nY[k+1][i] - nY[k][i] ) / ( fh[k] );
			}
#	else
				delta[ieq] -= convVeloY * ( fFDWMinus[k] * nY[k-1][i]
							+ fFDWCurr[k] * nY[k][i]
							+ fFDWPlus[k] * nY[k+1][i] );
#	endif

#	ifdef MOLARDIFF
		// 6. molar mass convection
			Double	convVeloM = 0.25 / ( Le[i] * rho[k] )
				* ( ( fFDWMinus[k] * rhoChiPrev * nY[k-1][i] / mixMolarMass[k-1]
					+ fFDWCurr[k] * rhoChiCurr * nY[k  ][i] / mixMolarMass[k  ]
					+ fFDWPlus[k] * rhoChiNext * nY[k+1][i] / mixMolarMass[k+1] )
				+ rhoChiCpOverLambda 
						* ( fFDWMinus[k] * lambdaOverCp[k-1] * nY[k-1][i] / mixMolarMass[k-1]
							+ fFDWCurr[k] * lambdaOverCp[k  ] * nY[k  ][i] / mixMolarMass[k  ]
							+ fFDWPlus[k] * lambdaOverCp[k+1] * nY[k+1][i] / mixMolarMass[k+1] ) );

#		ifdef UPWIND
			if ( convVeloM > 0.0 ) {
				delta[ieq] -= convVeloM * ( mixMolarMass[k] - mixMolarMass[k-1] ) / ( fhm[k] );
			}
			else {
				delta[ieq] -= convVeloM * ( mixMolarMass[k+1] - mixMolarMass[k] ) / ( fh[k] );
			}
#		else
			delta[ieq] -= convVeloM * ( fFDWMinus[k] * mixMolarMass[k-1]
							+ fFDWCurr[k] * mixMolarMass[k]
							+ fFDWPlus[k] * mixMolarMass[k+1] );
#		endif
#	endif
#	ifdef DIFFCORR
		// 7. mass fraction diffusivity correction convection
			delta[ieq] += 0.25 * sumdYdZOverLe / rho[k] * (  
						( fFDWMinus[k] * rhoChiPrev * nY[k-1][i]
						+ fFDWCurr[k] * rhoChiCurr * nY[k  ][i]
						+ fFDWPlus[k] * rhoChiNext * nY[k+1][i] )
					+ rhoChiCpOverLambda 
						* ( fFDWMinus[k] * lambdaOverCp[k-1] * nY[k-1][i]
							+ fFDWCurr[k] * lambdaOverCp[k] * nY[k  ][i]
							+ fFDWPlus[k] * lambdaOverCp[k+1] * nY[k+1][i] ) );
#		ifdef MOLARDIFF
		// 8. molar mass diffusivity correction convection
			Double	sumDiffCorrMM1 = 0.0;
			Double	sumDiffCorrMM2 = 0.0;
			Double	rhoChiYOverMPrev = rhoChiPrev * nY[k-1][i] / mixMolarMass[k-1];
			Double	rhoChiYOverMCurr = rhoChiCurr * nY[k  ][i] / mixMolarMass[k  ];
			Double	rhoChiYOverMNext = rhoChiNext * nY[k+1][i] / mixMolarMass[k+1];
			Double	lambdaOCpYOMPrev = lambdaOverCp[k-1] * nY[k-1][i] / mixMolarMass[k-1];
			Double	lambdaOCpYOMCurr = lambdaOverCp[k  ] * nY[k  ][i] / mixMolarMass[k  ];
			Double	lambdaOCpYOMNext = lambdaOverCp[k+1] * nY[k+1][i] / mixMolarMass[k+1];
			for ( int j = 0; j < nSpeciesIn; ++j ) {
				sumDiffCorrMM1 += ( fFDWMinus[k] * rhoChiYOverMPrev * nY[k-1][j]
						+ fFDWCurr[k] * rhoChiYOverMCurr * nY[k  ][j]
						+ fFDWPlus[k] * rhoChiYOverMNext * nY[k+1][j] ) / Le[j];
				sumDiffCorrMM2 += ( fFDWMinus[k] * lambdaOCpYOMPrev * nY[k-1][j]
								+ fFDWCurr[k] * lambdaOCpYOMCurr * nY[k  ][j]
								+ fFDWPlus[k] * lambdaOCpYOMNext * nY[k+1][j] ) / Le[j];
			}		

#		ifdef UPWIND
			fprintf( stderr, "upwind not completely implemented\n" );
			if ( convVeloM > 0.0 ) {
				delta[ieq] -= convVeloM * ( mixMolarMass[k] - mixMolarMass[k-1] ) / ( fhm[k] );
			}
			else {
				delta[ieq] -= convVeloM * ( mixMolarMass[k+1] - mixMolarMass[k] ) / ( fh[k] );
			}
#		else
			delta[ieq] += 0.25 / rho[k] * ( sumDiffCorrMM1 + rhoChiCpOverLambda * sumDiffCorrMM2 ) * dWdZ;
#		endif
#		endif
#	endif
#endif
			sum += hi[k][i] * mDoti[k][i];
			sumYH += hi[k][i] * nY[k][i];
            sumCpidYidZ += ( 
#ifdef DIFFCORR
				- cpMix[k] 
#endif
				+ cpi[k][i] ) / Le[i] * ( fFDWMinus[k] * nY[k-1][i]
                                                    + fFDWCurr[k]  * nY[k  ][i]
                                                    + fFDWPlus[k]  * nY[k+1][i]
#ifdef MOLARDIFF
							+ nY[k][i] / mixMolarMass[k] * dWdZ
#endif
				);
// #ifndef NEWCONVECTION
// #	ifdef CENTRALGRID
// 			delta[ieq] += - yPrime[gpGridOff] * ( fFDWMinus[k] * nY[k-1][i] + fFDWCurr[k] * nY[k][i] + fFDWPlus[k] * nY[k+1][i] );
// #	else

// 			if ( ( -yPrime[gpGridOff] ) > 0.0 ) {
// 				delta[ieq] += - ( nY[k][i] - nY[k-1][i] ) / fhm[k] * yPrime[gpGridOff];
// 			}
// 			else {
// 				delta[ieq] += - ( nY[k+1][i] - nY[k][i] ) / fh[k] * yPrime[gpGridOff];
// 			}	
// #	endif
// #endif
		}
		
		ieq = gpOff + fTemperature;
		
// Energy Equation
#ifdef TOTENT
		delta[ieq] = GetTotEnt( k, Z, *t ) - sumYH;
#else
		delta[ieq] = yPrime[ieq] 
				+ ( sum - fDPdt ) / ( rho[k] * cpMix[k] ) 
				- 0.5 * chi
				* ( fWMinus[k] * nTemp[k-1]
				+ fWCurr[k] * nTemp[k]
				+ fWPlus[k] * nTemp[k+1] )
				- GetTempSource( y[gpGridOff] * GetZR() );



#ifdef CENTRALTEMP                                                       
		delta[ieq] -= 0.5 * chi/ cpMix[k] * ( fFDWMinus[k] * nTemp[k-1]
											+fFDWCurr[k] * nTemp[k]
											+fFDWPlus[k] * nTemp[k+1] )
                                         * ( 
#	ifdef ENTFLUX
										 	sumCpidYidZ
#	else
											0.0
#	endif
#	ifdef HEATCAPGRAD
											+ fFDWMinus[k] * cpMix[k-1]
											+ fFDWCurr[k] * cpMix[k]
											+ fFDWPlus[k] * cpMix[k+1]
#	else
											+ 0.0
#	endif
											);
#else
#	ifdef ENTFLUX
		Double	convVeloEntFlux = 0.5 * chi/ cpMix[k] * sumCpidYidZ;
		if ( convVeloEntFlux < 0.0 ) {
			delta[ieq] -= convVeloEntFlux * ( nTemp[k] - nTemp[k-1] ) / ( fhm[k] );
		}
		else {
			delta[ieq] -= convVeloEntFlux * ( nTemp[k+1] - nTemp[k] ) / ( fh[k] );
		}
#	endif
#	ifdef HEATCAPGRAD
		Double	convVelocpdT = 0.5 * chi/ cpMix[k] * ( fFDWMinus[k] * cpMix[k-1]
											+ fFDWCurr[k] * cpMix[k]
											+ fFDWPlus[k] * cpMix[k+1] );
		if ( convVelocpdT < 0.0 ) {
			delta[ieq] -= convVelocpdT * ( nTemp[k] - nTemp[k-1] ) / ( fhm[k] );
		}
		else {
			delta[ieq] -= convVelocpdT * ( nTemp[k+1] - nTemp[k] ) / ( fh[k] );
		}
#	endif
#endif
		

		if ( fProperties->GetRadiation() ) {
				// wird zweimal aufgerufen, koennte aber in TFlame.cp auskommentiert werden : is called twice, but could be commented out in TFlame.cp

			/*********************Rui10112018***********************************************************/
			if (GasRadiationName == "Thin" )  {

					fProperties->GetRadiation()->SetRadiation( nTemp[k], nY[k], molarMass, rho[k] ); 

					delta[ieq] -= ( fProperties->GetRadiation()->GetRadiation() ) / ( rho[k] * cpMix[k] ); // This is the enthalpy defects due to gas radiation

			}	else if (GasRadiationName == "WSGG" || GasRadiationName == "Grey" || GasRadiationName == "SNB" || GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar") {

					delta[ieq] -= fRadLossPercent * ( - fProperties->GetRadiation()->dqrF[k] ) / ( rho[k] * cpMix[k] );   // Note dqrF starts from 0 and doesn't have B.C.

			}	else{
				cerr << "Error: gas radiation model: " << GasRadiationName << " is not found!!!" << NEWL;
			}
		/*********************Rui10112018*********************************************************************/

// 		   	if ( fSoot ) {
// #	ifdef NOSOOTRAD
// #	else
// #		ifdef ROSSELAND
// 			  delta[ieq] += GetRosseRadiation( k, &nTemp[k], &moments[k], rhoChiPrev / lambdaOverCp[k-1]
//                     , rhoChiCpOverLambda, rhoChiNext / lambdaOverCp[k+1] ) / ( rho[k] * cpMix[k] );
// #		else
// 				delta[ieq] += fSoot->GetSootRadiation( nTemp[k], moments[k] ) / ( rho[k] * cpMix[k] );
// #		endif
// #	endif
// 			}
		}
#endif


// #ifndef NEWCONVECTION
// #	ifdef CENTRALGRID
// 			delta[ieq] += - yPrime[gpGridOff] * ( fFDWMinus[k] * nTemp[k-1] + fFDWCurr[k] * nTemp[k] + fFDWPlus[k] * nTemp[k+1] );
// #	else
// 			if ( ( -yPrime[gpGridOff] ) > 0.0 ) {
// 				delta[ieq] += - ( nTemp[k] - nTemp[k-1] ) / fhm[k] * yPrime[gpGridOff];
// 			}
// 			else {
// 				delta[ieq] += - ( nTemp[k+1] - nTemp[k] ) / fh[k] * yPrime[gpGridOff];
// 			}	
// #	endif
// #endif

// Soot Equations                                     // commented out by Rui 01/24/2018
// 		if ( fSoot ) {
// 			/*if ( k > -1 )*/ 
// 			MOverRhoPrev = &y[(k-1)*fNOfEquations+sootOff];
// 			/*if ( k < fNGridPoints )*/ 
// 			MOverRhoNext = &y[(k+1)*fNOfEquations+sootOff];
// 			MOverRhoCurr = &y[k*fNOfEquations+sootOff];
// 			// UpdateSoot wird zwei mal aufgerufen
// 			// soot muss noch veraendert werden
// #ifdef NEWPROPERTIES
// 			fSoot->ComputePolymereConcs( nY[k], nTemp[k], rho[k], fSpecies.GetMolarMass()->vec
// 					, fSoot->GetPij()->mat, fSoot->GetSumPi()->vec, fSoot->GetPAHMoments()->vec
// 					, moments[k], &fReaction );
// 			fSoot->UpdateSoot( &fReaction, &fSpecies, moments[k], nTemp[k], nY[k], rho[k]
// 								, mixMolarMass[k] );
// #endif
// 			Double	fracIndex;
// 			gpSpecOff = gpOff + sootOff;
// #ifdef PRONE
// 			Double	PrCurr = 0.7;
// #else
// 			Double	PrCurr = mu[k] / lambdaOverCp[k];
// #endif
// 			for ( i = 0; i < nSootMoments; ++i ) {
// 				fracIndex = i - 2.0 / 3.0;
// 				ieq = gpSpecOff + i;
				
// 				delta[ieq] = yPrime[ieq] 
					 
// #ifdef SIZEDEPDIFFUSION
// 					- chi / ( 2.0 * fSoot->GetLewis1() ) * FRACMOMFACT * ( fWMinus[k] * fSoot->FracMom2( fracIndex, MOverRhoPrev )
// 					+ fWCurr[k] * fSoot->FracMom2( fracIndex, MOverRhoCurr )
// 					+ fWPlus[k] * fSoot->FracMom2( fracIndex, MOverRhoNext ) );
// #else
// 					- 0.5 * chi * 
// 					( fWMinus[k] * MOverRhoPrev[i]
// 					+ fWCurr[k] * MOverRhoCurr[i]
// 					+ fWPlus[k] * MOverRhoNext[i] );
// #endif
// 				if ( fSoot->WithThermoPhoresis() ) {
// #ifdef NEWTHERMOPHOR
// 					Double	hnenn = fhm[k] * fh[k] * ( fhm[k] + fh[k] );
// 					delta[ieq] -= 0.275 * PrCurr / rho[k] 
// 						* ( fhm[k] * ( rhoChiCurr * MOverRhoCurr[i] / nTemp[k]  + rhoChiNext * MOverRhoNext[i] / nTemp[k+1]  ) 
// 								* ( nTemp[k+1] - nTemp[k] )
// 							+ fh[k] * ( rhoChiCurr * MOverRhoCurr[i] / nTemp[k]  + rhoChiPrev * MOverRhoPrev[i] / nTemp[k-1]
// 								* ( nTemp[k-1] - nTemp[k] ) ) / hnenn );


					
// #ifdef CONVECTION
// 					Double	convVeloThPhor = 0.1375 * PrCurr / rho[k] * MOverRhoCurr[i] / nTemp[k] * constConvVeloY;
// 					if ( convVeloThPhor > 0.0 ) {
// 						delta[ieq] -= convVeloThPhor * ( nTemp[k+1] - nTemp[k] ) / ( fh[k] );
// 					}
// 					else {
// 						delta[ieq] -= convVeloThPhor * ( nTemp[k] - nTemp[k-1] ) / ( fhm[k] );
// 					}
// #endif

// #else
// 					delta[ieq] -= 0.275 * MOverRhoCurr[i] / ( nTemp[k] ) * PrCurr * chi
// 										 * ( fWMinus[k] * nTemp[k-1]
// 											+fWCurr[k] * nTemp[k]
// 											+fWPlus[k] * nTemp[k+1] );
// #endif
// 				}
// 				if ( fSoot->WithNucleation() ) {
// #ifdef TESTNEWNUCLEATION
// 						delta[ieq] -= fSoot->NucleationNow( i, nTemp[k], nY[k]
// 							, rho[k], molarMass ) / rho[k];
// #else
// 						delta[ieq] -= fSoot->NucleationNew( i, nTemp[k]
// 							, fSoot->GetPAHMoments()->vec ) / rho[k];
// #endif
// 				}
// 				if ( fSoot->WithCoagulation() ) {
// 						delta[ieq] -= fSoot->SourceCoagulationNew( i, nTemp[k]
// 										, MOverRhoCurr );
// 				}
// 				if ( fSoot->WithCondensation() ) {
// 					delta[ieq] -= fSoot->SourceCondensationNew( i, nTemp[k]
// 									, fSoot->GetPAHMoments()->vec, MOverRhoCurr,  nY[k], rho[k], molarMass );
// 				}
// 				if ( fSoot->WithSurfaceGrowth() ) {
// 					delta[ieq] -= fSoot->SourceSurfGrowthNew( i, MOverRhoCurr, nY[k], rho[k], molarMass );
// 				}
// 				if ( fSoot->WithSurfaceOxidation() ) {
// 					delta[ieq] -= fSoot->SourceSootOxidationNew( i, MOverRhoCurr, nY[k], rho[k], molarMass );
// 				}
// #ifdef SOOTCONVECTION
// 					Double	convVeloSoot = constConvVeloY;
// #	ifdef UPWINDSOOT
// #		ifdef SECORDSOOTUPWIND
// 					Double	h1, h2;
// 					if ( convVeloSoot > 0.0 ) {
// 						if ( k == 0 ) {
// 							delta[ieq] += convVeloSoot * ( ( MOverRhoCurr[i] - MOverRhoPrev[i] ) 
// 									/ ( fhm[k] )
// 									- 1.0 / fSoot->GetLewis1() 
// 									* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, MOverRhoNext )
// 									- fSoot->FracMom2( fracIndex, MOverRhoCurr ) ) 
// 										/ ( fh[k] ) );
// 						}
// 						else {
// 							h1 = fhm[k];
// 							h2 = fhm[k]+fhm[k-1];
// 							delta[ieq] += convVeloSoot
// 								* ( ( h2 * h2 - h1 * h1 ) * MOverRhoCurr[i]
// 									- h2 * h2 * MOverRhoPrev[i] 
// 									+ h1 * h1 * MOverRhoPrev[i-fNOfEquations] ) 
// 									/ ( h1 * h2 * ( h2 - h1 ) );
// 						}
// 					}
// 					else {
// 						if ( k == (fNGridPoints-1) ) {
// 							delta[ieq] += convVeloSoot * ( ( MOverRhoNext[i] - MOverRhoCurr[i] ) 
// 									/ ( fh[k] )
// 									- 1.0 / fSoot->GetLewis1() 
// 									* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, MOverRhoCurr )
// 									- fSoot->FracMom2( fracIndex, MOverRhoPrev ) ) 
// 										/ ( fhm[k] ) );
// 						}
// 						else {
// 							h1 = fh[k];
// 							h2 = fh[k]+fh[k+1];
// 							delta[ieq] += convVeloSoot
// 								* ( -( h2 * h2 - h1 * h1 ) * MOverRhoCurr[i] 
// 									+ h2 * h2 * MOverRhoNext[i]
// 									- h1 * h1 * MOverRhoNext[i+fNOfEquations] ) 
// 									/ ( h1 * h2 * ( h2 - h1 ) );
// 						}
// 					}
// #		else
// 					if ( convVeloSoot > 0.0 ) {
// 						delta[ieq] += convVeloSoot * ( ( MOverRhoCurr[i] - MOverRhoPrev[i] ) 
// 								/ ( fhm[k] )
// 								- 1.0 / fSoot->GetLewis1() 
// 								* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, MOverRhoNext )
// 								- fSoot->FracMom2( fracIndex, MOverRhoCurr ) ) 
// 									/ ( fh[k] ) );
// 					}
// 					else {
// 						delta[ieq] += convVeloSoot * ( ( MOverRhoNext[i] - MOverRhoCurr[i] ) 
// 								/ ( fh[k] )
// 								- 1.0 / fSoot->GetLewis1() 
// 								* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, MOverRhoCurr )
// 								- fSoot->FracMom2( fracIndex, MOverRhoPrev ) ) 
// 									/ ( fhm[k] ) );
// 					}
// #		endif
// #	else
// 					// central differences
// #		ifdef NEWSOOTTERMDISC
// 					Double	hnenn = fhm[k] * fh[k] * ( fhm[k] + fh[k] );
// 					delta[ieq] += 0.25 / rho[k] * ( 
// 							( fhm[k] * ( rhoChiCurr + rhoChiNext ) * ( MOverRhoNext[i] - MOverRhoCurr[i] )
// 							+ fh[k] * ( rhoChiCurr + rhoChiPrev ) * ( MOverRhoPrev[i] - MOverRhoCurr[i] ) ) / hnenn
// 							+ rhoChiCurr / lambdaOverCp[k] 
// 							* ( fhm[k] * ( lambdaOverCp[k] + lambdaOverCp[k+1] ) * ( MOverRhoNext[i] - MOverRhoCurr[i] )
// 							+ fh[k] * ( lambdaOverCp[k] + lambdaOverCp[k-1] ) * ( MOverRhoPrev[i] - MOverRhoCurr[i] ) ) / hnenn
// 							- 2.0 * rhoChiCurr * ( fWMinus[k] * MOverRhoPrev[i]
// 											+ fWCurr[k] * MOverRhoCurr[i]
// 											+ fWPlus[k] * MOverRhoNext[i] )
// 							);
// #		else
// 					delta[ieq] -= convVeloSoot * ( 
// 							1.0 / fSoot->GetLewis1() 
// 							* FRACMOMFACT * ( fFDWMinus[k] * fSoot->FracMom2( fracIndex, MOverRhoPrev )
// 							+ fFDWCurr[k] * fSoot->FracMom2( fracIndex, MOverRhoCurr )
// 							+ fFDWPlus[k] * fSoot->FracMom2( fracIndex, MOverRhoNext ) )
							
// 							- ( fFDWMinus[k] * MOverRhoPrev[i]
// 							+ fFDWCurr[k] * MOverRhoCurr[i]
// 							+ fFDWPlus[k] * MOverRhoNext[i] )
// 							);
// #		endif
// #	endif
// //				}
// #endif

				
// #ifdef DEBUGSOOTSOURCE
// 				if ( fSoot->WithSurfaceOxidation() &&  fSoot->WithCondensation() && i == 1 ) {
// 					fprintf( fOutFilePtr, "%g\t%g\t%g\t%g\t%g\t%g\n", Z[k]
// 									, fSoot->NucleationNew( i, nTemp[k], fSoot->GetPAHMoments()->vec )
// 									, fSoot->SourceCondensationNew( i, nTemp[k], fSoot->GetPAHMoments()->vec
// 									, moments[k],  nY[k], rho[k], molarMass ) 
// 									,
// 									fSoot->SourceSurfGrowthNew( i, moments[k], nY[k], rho[k], molarMass ) 
// 									,
// 									fSoot->SourceSootOxidationNew( i, moments[k], nY[k], rho[k], molarMass )
// 									,
// 									chi / ( 2.0 * fSoot->GetLewis1() ) 
// #	ifdef SIZEDEPDIFFUSION
// 					* ( fWMinus[k] * fSoot->FracMom2( fracIndex, &y[gpOff-fNOfEquations+sootOff] )
// 					+ fWCurr[k] * fSoot->FracMom2( fracIndex, &y[gpOff+sootOff] )
// 					+ fWPlus[k] * fSoot->FracMom2( fracIndex, &y[gpOff+fNOfEquations+sootOff] ) );
// #	else
// 					* ( fWMinus[k] * y[ieq-fNOfEquations]
// 					+ fWCurr[k] * y[ieq]
// 					+ fWPlus[k] * y[ieq+fNOfEquations] )
// #	endif
// 						);
// 				}
// #endif
// #ifndef NEWCONVECTION
// #	ifdef CENTRALGRID
// 				delta[ieq] += - yPrime[gpGridOff] * ( fFDWMinus[k] * y[ieq-fNOfEquations] 
// 								+ fFDWCurr[k] * y[ieq] + fFDWPlus[k] * y[ieq+fNOfEquations] );
// #	else
// 				if ( ( -yPrime[gpGridOff] ) > 0.0 ) {
// 					delta[ieq] += - ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k] * yPrime[gpGridOff];
// 				}
// 				else {
// 					delta[ieq] += - ( y[ieq+fNOfEquations] - y[ieq] ) / fh[k] * yPrime[gpGridOff];
// 				}	
// #	endif
// #endif
// 			}
// 		}

// Grid                 // Commented out by Rui 01/24/2018
// #ifdef MOVEGRID
// 		if ( k == 0 ) {
// 			delta[gpGridOff] = - 2.0 * yPrime[gpGridOff] + yPrime[gpGridOff+fNOfEquations]; 
// 		}
// 		else if ( k == 1 ) {

// 			delta[gpGridOff] = fTauGrid * (
// 				  ( - kappa_ * gpDens[k-1] * gpDens[k-1] / fMonFct[k]
// 			  	+ ( - kappa_ * gpDens[k-2] * gpDens[k-2] 
// 				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff-fNOfEquations]
// 			  + ( ( kappa_ * gpDens[k-1] * gpDens[k-1]
//     			  + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] 
// 				+ ( kappa_ * gpDens[k] * gpDens[k] 
// 	  			+ ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff]
// 			  + ( - kappa_ * gpDens[k] * gpDens[k] / fMonFct[k-1]
// 			  	+ ( - kappa_ * gpDens[k+1] * gpDens[k+1] 
// 				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] ) * yPrime[gpGridOff+fNOfEquations]
// 			  + ( kappa_ * gpDens[k+1] * gpDens[k+1] / fMonFct[k] ) * yPrime[gpGridOff+2*fNOfEquations]
// 			  )
// 			  + ( - kappa_ * gpDens[k-1] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] - kappa_ * gpDens[k+1] ) / fMonFct[k]
// 			  - ( - kappa_ * gpDens[k-2] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] - kappa_ * gpDens[k] ) / fMonFct[k-1];	
// 		}	
// 		else if ( k == fNGridPoints-1 ) {
// 			delta[gpGridOff] = - 2.0 * yPrime[gpGridOff] + yPrime[gpGridOff-fNOfEquations];
// 		}
// 	 	else if ( k == fNGridPoints-2 ) {

// 			delta[gpGridOff] = fTauGrid * (
// 				( kappa_ * gpDens[k-2] * gpDens[k-2] / fMonFct[k-1] ) * yPrime[gpGridOff-2*fNOfEquations]
// 			  + ( - kappa_ * gpDens[k-1] * gpDens[k-1] / fMonFct[k]
// 			  	+ ( - kappa_ * gpDens[k-2] * gpDens[k-2] 
// 				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff-fNOfEquations]
// 			  + ( ( kappa_ * gpDens[k-1] * gpDens[k-1]
// 			      + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] 
// 				+ ( kappa_ * gpDens[k] * gpDens[k] 
// 				  + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff]
// 			  + ( - kappa_ * gpDens[k] * gpDens[k] / fMonFct[k-1]
// 			  	+ ( - kappa_ * gpDens[k+1] * gpDens[k+1] 
// 				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] ) * yPrime[gpGridOff+fNOfEquations]
// 			  )
// 			  + ( - kappa_ * gpDens[k-1] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] - kappa_ * gpDens[k+1] ) / fMonFct[k]
// 			  - ( - kappa_ * gpDens[k-2] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] - kappa_ * gpDens[k] ) / fMonFct[k-1];	
// 		}	
// 		else if ( k > 1 && k < fNGridPoints-2 ) {

// 			delta[gpGridOff] = fTauGrid * (
// 				( kappa_ * gpDens[k-2] * gpDens[k-2] / fMonFct[k-1] ) * yPrime[gpGridOff-2*fNOfEquations]
// 			  + ( - kappa_ * gpDens[k-1] * gpDens[k-1] / fMonFct[k]
// 			  	+ ( - kappa_ * gpDens[k-2] * gpDens[k-2] 
// 	    		  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff-fNOfEquations]
//   			+ ( ( kappa_ * gpDens[k-1] * gpDens[k-1]
//   			    + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] 
// 				+ ( kappa_ * gpDens[k] * gpDens[k] 
// 			  + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff]
// 			  + ( - kappa_ * gpDens[k] * gpDens[k] / fMonFct[k-1]
// 			  	+ ( - kappa_ * gpDens[k+1] * gpDens[k+1] 
// 				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] ) * yPrime[gpGridOff+fNOfEquations]
// 			  + ( kappa_ * gpDens[k+1] * gpDens[k+1] / fMonFct[k] ) * yPrime[gpGridOff+2*fNOfEquations]
//  			 )
//  			 + ( - kappa_ * gpDens[k-1] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] - kappa_ * gpDens[k+1] ) / fMonFct[k]
// 			 - ( - kappa_ * gpDens[k-2] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] - kappa_ * gpDens[k] ) / fMonFct[k-1];	
//   		}
// #else
// 		delta[gpGridOff] = yPrime[gpGridOff];  // Grid will not be changed
// #endif

#ifdef DEBUGRES
		fprintf( stderr, "gp = %d\n", k );
		for ( int j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
			fprintf( stderr, "%s\t%g\t%g\n", fVariableNames[j]
									, y[gpOff + j]
									, delta[gpOff + j] );
		}
		fprintf( stderr, "\n" );
#endif
	} 


// boundary conditions
	for ( k = -1; k <= fNGridPoints; k+=fNGridPoints+1 ) {
		gpOff = k*fNOfEquations;
		gpSpecOff = gpOff + fFirstSpecies;
    	gpGridOff = gpOff + fGrid;	

//	species
		sum = 0.0;
		sumYH = 0.0;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			delta[gpSpecOff+i] = yPrime[gpSpecOff+i];
			delta[gpSpecOff+i] -= mDoti[k][i] / rho[k];
			sum += hi[k][i] * mDoti[k][i];
#	ifdef TOTENT
			sumYH += hi[k][i] * nY[k][i];
#	endif
		}

//	temperature
#	ifdef TOTENT
		delta[gpOff+fTemperature] = GetTotEnt( k, Z, *t ) - sumYH;
#	else

		delta[gpOff+fTemperature] = yPrime[gpOff+fTemperature];
		delta[gpOff+fTemperature] += ( sum - fDPdt ) / ( rho[k] * cpMix[k] ) 
				- GetTempSource( y[gpOff+fGrid] * GetZR() );
		
		if ( fProperties->GetRadiation() ) {
			// wird zweimal aufgerufen, koennte aber in TFlame.cp auskommentiert werden : is called twice, but could be commented out in TFlame.cp
		
		/*********************Rui10112018**********************************************************************************/
			if (GasRadiationName == "Thin" )  {

				fProperties->GetRadiation()->SetRadiation( nTemp[k], nY[k], molarMass, rho[k] ); 

				delta[gpOff+fTemperature] -= ( fProperties->GetRadiation()->GetRadiation() ) / ( rho[k] * cpMix[k] );   // Updated boundary radiation for thin

			} else if (GasRadiationName == "WSGG" || GasRadiationName == "Grey" || GasRadiationName == "SNB" || GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar") {

				// delta[gpOff+fTemperature] -= ( - fProperties->GetRadiation()->dqrF[k+1] ) / ( rho[k] * cpMix[k] ); // No radiation value at Bounday
			} else{
				cerr << "Error: gas radiation model: " << GasRadiationName << " is not found!!!" << NEWL;
			}

		/*********************Rui10112018***********************************************************************************/

// 			if ( fSoot ) {    //commented out by Rui 01/24/2018
// #		ifdef NOSOOTRAD
// #		else
// 				delta[gpOff+fTemperature] += ( fSoot->GetSootRadiation( nTemp[k], moments[k] ) ) 
// 												/ ( rho[k] * cpMix[k] );
// #		endif
// 			}
		}
#	endif

//	soot moments
		// if ( fSoot ) {             /commented out by Rui 01/24/2018 
		// 	gpSpecOff = gpOff + sootOff;
		// 	for ( i = 0; i < nSootMoments; ++i ) {
		// 		ieq = gpSpecOff + i;
		// 		delta[ieq] = yPrime[ieq];
		// 		fSoot->ComputePolymereConcs( nY[k], nTemp[k], rho[k], fSpecies.GetMolarMass()->vec
		// 				, fSoot->GetPij()->mat, fSoot->GetSumPi()->vec, fSoot->GetPAHMoments()->vec
		// 				, moments[k], &fReaction );
		// 		fSoot->UpdateSoot( &fReaction, &fSpecies, moments[k], nTemp[k], nY[k], rho[k]
		// 							, mixMolarMass[k] );
		// 		if ( fSoot->WithNucleation() ) {
		// 			delta[ieq] -= fSoot->NucleationNew( i, nTemp[k]
		// 					, fSoot->GetPAHMoments()->vec ) / rho[k];
		// 		}
		// 		if ( fSoot->WithCoagulation() ) {
		// 			delta[ieq] -= fSoot->SourceCoagulationNew( i, nTemp[k]
		// 								, moments[k] ) / rho[k];
		// 		}
		// 		if ( fSoot->WithCondensation() ) {
		// 			delta[ieq] -= fSoot->SourceCondensationNew( i, nTemp[k]
		// 							, fSoot->GetPAHMoments()->vec, moments[k],  nY[k], rho[k], molarMass ) 
		// 							/ rho[k];
		// 		}
		// 		if ( fSoot->WithSurfaceGrowth() ) {
		// 			delta[ieq] -= fSoot->SourceSurfGrowthNew( i, moments[k], nY[k], rho[k], molarMass ) 
		// 							/ rho[k];
		// 		}
		// 		if ( fSoot->WithSurfaceOxidation() ) {
		// 			delta[ieq] -= fSoot->SourceSootOxidationNew( i, moments[k], nY[k], rho[k], molarMass )
		// 							/ rho[k];
		// 		}
		// 	}
		// }
	
//	grid equation
		delta[gpOff+fGrid] = yPrime[gpOff+fGrid];
	}  // Finish B.C. 	


#ifdef MOVEZRIGHT
	Double	ZR = GetZR();
#	ifdef CENTRALZRCONV
// start of central part
	for ( k = 0; k < fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		for ( i = 0; i < fNOfEquations; ++i ) {
			if ( i != fGrid ) {
				ieq = gpOff+i;
#		ifdef NEWCONVECTION
				Double	theConvVelo = /*-y[gpOff+fGrid] / ZR * fdDeltaZdt*/ -yPrime[gpOff + fGrid];
				delta[ieq] += theConvVelo * ( 
						fFDWMinus[k] * y[ieq-fNOfEquations]
						+ fFDWCurr[k] * y[ieq]
						+ fFDWPlus[k] * y[ieq+fNOfEquations] );
#		else
				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
							* ( fFDWMinus[k] * y[ieq-fNOfEquations]
								+fFDWCurr[k] * y[ieq]
								+fFDWPlus[k] * y[ieq+fNOfEquations] );
#		endif
			}
		}
	}
	gpOff = fNGridPoints * fNOfEquations;
	for ( i = 0; i < fNOfEquations; ++i ) {
		if ( i != fGrid ) {
			ieq = gpOff+i;
#		ifdef NEWCONVECTION
			Double	theConvVelo = /*-y[gpOff+fGrid] / ZR * fdDeltaZdt */-yPrime[gpOff + fGrid];
			delta[ieq] += theConvVelo * ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k];
#		else
			delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
						* ( y[ieq] - y[ieq-fNOfEquations] ) / ( fhm[k] );
#		endif
		}
	}
// end of central part
// start of upwind part
	for ( k = 0; k <= fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		for ( i = 0; i < fNOfEquations; ++i ) {
			if ( i != fGrid ) {
				ieq = gpOff+i;
#ifdef NEWCONVECTION
				Double	theConvVelo = -y[gpOff+fGrid] / ZR * fdDeltaZdt/* -yPrime[gpOff + fGrid]*/;
				if ( theConvVelo > 0.0 || k == fNGridPoints ) {
					delta[ieq] += theConvVelo * ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k];
				}
				else {
					delta[ieq] += theConvVelo * ( y[ieq+fNOfEquations] - y[ieq] ) / fh[k];
				}	
#else
				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
							* ( y[ieq] - y[ieq-fNOfEquations] ) / ( fhm[k] );
#endif
			}
		}
	}
// end of upwind part
#	else
	for ( k = 0; k <= fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		for ( i = 0; i < fNOfEquations; ++i ) {
			if ( i != fGrid ) {
				ieq = gpOff+i;
#ifdef NEWCONVECTION
				Double	theConvVelo = -y[gpOff+fGrid] / ZR * fdDeltaZdt -yPrime[gpOff + fGrid];
				if ( theConvVelo > 0.0 || k == fNGridPoints ) {
					delta[ieq] += theConvVelo * ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k];
				}
				else {
					delta[ieq] += theConvVelo * ( y[ieq+fNOfEquations] - y[ieq] ) / fh[k];
				}	
#else
				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
							* ( y[ieq] - y[ieq-fNOfEquations] ) / ( fhm[k] );
#endif
			}
		}
	}
#	endif
#endif

#	ifdef READU			
	for ( kU = -1; kU <= fNGridPoints; ++kU ) {
		ieqU = kU * fNOfEquations;
		UstOverU = GetU( *t, Z[kU] );
		for ( iU = 0; iU < fNOfEquations; ++iU ) {
			yPrime[ieqU+iU] *= UstOverU;
			delta[ieqU+iU] *= UstOverU;
		}
	}
#	endif

#	ifdef LEWISCHANGE
	delete Le;
#	endif
}

template<typename Flame>
int ResTransFlameImpliSolverCV( realtype T, N_Vector u, N_Vector udot, void *data )
{
	Flame* flame = (Flame*) data;
	int		idum, iRes = 0;
	int		i, nEq = ( flame->fNGridPoints + 2 ) * flame->fNOfEquations;
	Double	ddum;
	
	Double	*y = NV_DATA_S(u);
	Double	*delta = NV_DATA_S(udot);

	Double	*yPrime = flame->fSolPrime->mat[kPrev];

	for ( i = 0; i < nEq; ++i ) {
		yPrime[i] = 0.0;
	}

	flame->ResTransFlameImpliSolver( &T, y, yPrime, delta, &idum, &ddum, &idum );

	for ( i = 0; i < nEq; ++i ) {
		delta[i] = -delta[i];
	}

	if ( iRes < 0 ) {
		return 1; // means recoverable error
	}
	else {
		return 0; // successful
	}
}

template<typename Flame>
void JacTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *pd
			, Double *cj, Double *rPar, int *iPar )
{
	Flame*	flame = ( Flame* ) iPar;
	int						nEq = flame->fNOfEquations;
	MatrixPtr	pdNew = FortranToCMat( pd, nEq, nEq, nEq, nEq );
	( ( Flame* ) iPar )->JacTransFlameSolver( T, y, yPrime, pdNew->mat
			, *cj, rPar, iPar );
	DisposeFToCMat( pdNew );
}

template<typename Species>
void TTransFlameSolver<Species>::JacTransFlameSolver( Double * /*T*/, Double *y, Double * /*yPrime*/     //Not used for me
			, Double **pd, Double cj, Double * /*rPar*/, int * /*iPar*/ )
{
	int						nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int						i, j, ieq, jeq, kAct = fActualPoint;
	Double					sum = 0.0;	

	Double					*prodRate = fSpecies.GetProductionRate()->vec;
	Double					*enth = fSpecies.GetEnthalpy()->vec;
	Double					*cp = fSpecies.GetHeatCapacity()->vec;
	Double					*molarMass = fSpecies.GetMolarMass()->vec;
	Double					*Le = fSpecies.GetLewisNumber()->vec;

	Double					temp = y[fTemperature];
	Double					*YF = &y[fFirstSpecies];

	Double					*nTime = &fSolTime->vec[kAct];
	Double					Z = fSolGrid->vec[kAct];
	Double					chi = GetDissRate( nTime[kCurr], Z );
	Double					pressure = GetPressure( nTime[kCurr] );
	Double					**dMdY = fDmDy;
	Double					*dMdT = fDmDy[nSpeciesIn];
							// dMdY[i][j] = dM_j/dY_i

	if ( fSoot ) {
		fprintf( fOutFilePtr, "###error: analytical evaluation of jacobian for soot not yet implemented\n" );
		exit( 2 );
	}

	T0DFlame<Species>::UpdateThermoProps( YF, temp, pressure, fProperties->GetDensityRef()   
									, kDensFromPress, NULL );

	Double				rho = fProperties->GetDensity();
	Double				heatCap = fProperties->GetMixHeatCapacity();
	Double				mixMolarMass = fProperties->GetMixMolarMass();

	fReaction.FilldMdYOnePointAnal( dMdY, YF, fReaction.GetReactionRate()->vec
						, mixMolarMass, molarMass, fReaction.GetTBConc()->vec
						, prodRate, rho );
	fReaction.FilldMdTOnePointAnal( dMdT, temp, fReaction.GetReactionRate()->vec
						, fReaction.GetRateCoefficients()->vec, pressure, molarMass
						, fReaction.GetTBConc()->vec );

//	fill dG_i/dY_j
	for ( i = 0; i < nSpeciesIn; ++i ) {
		ieq = fFirstSpecies + i;
		pd[ieq][ieq] += cj;
		pd[ieq][ieq] -= chi / ( 2.0 * Le[i] ) * fWCurr[kAct];
		pd[ieq][fTemperature] -= ( dMdT[i] + prodRate[i] / temp ) / rho;
		for ( j = 0; j < nSpeciesIn; ++j ) {
			jeq = fFirstSpecies + j;
			pd[ieq][jeq] -= ( dMdY[j][i] + prodRate[i] * mixMolarMass / molarMass[j] ) / rho;
		}
		sum += enth[i] * prodRate[i];
	}

	Double	oneOverRhoCp = 1.0 / ( rho * heatCap );
	
	pd[fTemperature][fTemperature] += cj;
	pd[fTemperature][fTemperature] -= 0.5 * chi * fWCurr[kAct];
	pd[fTemperature][fTemperature] +=
			sum / temp * oneOverRhoCp;
	for ( i = 0; i < nSpeciesIn; ++i ) {
		pd[fTemperature][fTemperature] +=
				( enth[i] * dMdT[i] + cp[i] * prodRate[i] ) * oneOverRhoCp;
		pd[fTemperature][fFirstSpecies+i] +=
				sum * mixMolarMass / molarMass[i] * oneOverRhoCp;
		for ( j = 0; j < nSpeciesIn; ++j ) {
			jeq = fFirstSpecies + j;
			pd[fTemperature][jeq] +=
					enth[i] * dMdY[j][i] * oneOverRhoCp;
		}
	}

#ifdef DEBUGRES
	for ( j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
		fprintf( stderr, "%s\t%g\n", fVariableNames[j], pd[fTemperature][j] );
	}
	fprintf( stderr, "\n" );
#endif
	
}


template<typename Species>
Double TTransFlameSolver<Species>::GetPressure( Double t )
{
	return fPressStart + ( fPressEnd - fPressStart ) / MAX( 1e-30, fTEnd - fTStart ) 
					* ( t - fTStart );
}

template<typename Species>
Double TTransFlameSolver<Species>::GetTempSource( Double Z )
{
#ifdef DPDT
	return 0.0;
#else
	return fDTdtOx + Z * ( fDTdtFu - fDTdtOx );
#endif
}

template<typename Species>
Double TTransFlameSolver<Species>::Interpol( Double t, Double valOld, Double tOld, Double valNew, Double tNew )
{
	if ( valNew == valOld ) {
		return valNew;
	}
	return valOld + ( valNew - valOld ) / ( tNew - tOld ) * ( t - tOld );
}

#ifdef READCHI
#	ifdef TIMEDEPCHI
template<typename Species>
Double TTransFlameSolver<Species>::GetDissRate( Double t, Double Z )
{
	int	i = 1;
	int	j = 1;
	int	zlen = fZIn->len;
	int	tlen = fTimeIn->len;
	Double	*zIn = fZIn->vec;
	Double	*tIn = fTimeIn->vec;
	Double	**chiIn = fChiIn->mat;
	Double	chiAtZt1, chiAtZt2;


	while ( j < zlen-1 && Z > zIn[j] ) ++j;
	if ( j == zlen-1 && Z > zIn[j]+1.0e-8 ) {
		fprintf( stderr, "Z = %g out of range[%g,%g]\n", Z, zIn[0]
								, zIn[zlen-1] );
	}
	
	while ( i < tlen-1 && t > tIn[i] ) ++i;
	if ( i == tlen-1 && t > tIn[i] ) {
		fprintf( stderr, "t = %g out of range[%g,%g]: linear extrapolation\n", Z, tIn[0]
								, tIn[tlen-1] );
	}
	
	if ( zIn[j] == zIn[j-1] ) {
		fprintf( stderr, "something's wrong\n" );
	}
	else {
		chiAtZt1 = Interpol( Z, chiIn[i-1][j-1], zIn[j-1], chiIn[i-1][j], zIn[j] );
		chiAtZt2 = Interpol( Z, chiIn[i][j-1], zIn[j-1], chiIn[i][j], zIn[j] );
	}
	
	if ( tIn[i] == tIn[i-1] ) {
		fprintf( stderr, "something else's wrong\n" );
		return 0.0;
	}
	else {
		return Interpol( t, chiAtZt1, tIn[i-1], chiAtZt2, tIn[i] );
	}
}

#	else // !defined TIMEDEPCHI
template<typename Species>
Double TTransFlameSolver<Species>::GetDissRate( Double /*t*/, Double z )
{
	int	i = 1;
	int	len = fZCount->len;
	Double	*zCount = fZCount->vec;
	Double	*chiCount = fChiCount->vec;

	while ( i < len && z > zCount[i] ) ++i;
	if ( i == len && z > zCount[len-1] ) {
		fprintf( stderr, "z = %g out of range[%g,%g]\n", z, zCount[0]
								, zCount[len-1] );
	}
	
	if ( zCount[i] == zCount[i-1] ) {
		return 0.0;
	}
	else {
		return chiCount[i-1] + ( chiCount[i] - chiCount[i-1] ) / ( zCount[i] - zCount[i-1] ) 
											* ( z - zCount[i-1] );
	}
}
#	endif // !defined TIMEDEPCHI
#else // !defined READCHI

#ifdef MOVEZRIGHT
template<typename Species>
Double TTransFlameSolver<Species>::GetDissRate( Double t, Double z )
{
	if ( z < 1.0e-10 || z >= 1.0 ) {
		return 0.0;
	}

	Double	zetaRef = fZRef / Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );

	return GetRefDissRate( t ) * z * z / ( zetaRef * zetaRef ) * log( z ) / log( zetaRef );
}
#else // !defined MOVEZRIGHT
template<typename Species>
Double TTransFlameSolver<Species>::GetDissRate( Double t, Double z )
{
	return GetRefDissRate( t ) * ExactChi( z ) / ExactChi( fZRef );           // This one is being used!!!!!!!!!!!!!
}
#		endif // defined MOVEZRIGHT
#endif // READCHI

template<typename Species>
Double TTransFlameSolver<Species>::GetRefDissRate( Double t )
{
	return fRandomNumber * Interpol( t, fChiStart, fTStart, fChiEnd, fTEnd );
}

template<typename Species>
Double TTransFlameSolver<Species>::GetDissRate( Double t, Double z, Double rho )
{
	Double	fRhoRef = 0.0;
	fprintf( stderr, "error in GetDissRate\n" );
	exit(2);
	return DissRateFact( t, rho ) / DissRateFact( t, fRhoRef ) * GetDissRate( t, z );
}

template<typename Species>
Double TTransFlameSolver<Species>::ExactChi( Double Z )
{
	int				i;
	const int		nSteps = 100;
	Double			twoZ = 2.0 * Z;
	Double			deltax = 4.0 / ( ( Double ) nSteps ), twoErfc;
	static Double	erFunc[4*nSteps];
	static Flag		init = FALSE;
	
	if ( !init ) {
		init = TRUE;
		for ( i = 0; i < 4 * nSteps; ++i ) {
			erFunc[i] = erfc( i*deltax );
		}
	}
	
	if ( Z > 0.5 ) {
		twoZ = 2.0 * ( 1.0 - Z );
	}
	if ( Z < 1.0e-7 ) {
		return 0.0;
	}
	
	i = 0;
	while ( erFunc[i] > twoZ ) ++i;
	
	twoErfc = Interpol( twoZ, (i-1)*deltax, erFunc[i-1], i*deltax, erFunc[i] );
	
	return exp( -2.0 * twoErfc * twoErfc );
}

template<typename Species>
Double TTransFlameSolver<Species>::DissRateFact( Double t, Double rho )
{
	Double rhoInf = fWOverRInf * GetPressure( t ) / Interpol( t, fSolOldTemp->vec[kPrev], fSolOldTime->vec[kPrev]
									, fSolTemp->vec[kPrev], fSolTime->vec[kPrev] );
	Double	densRatio = sqrt( rhoInf / rho );

	return 1.5 * ( densRatio + 1 ) * ( densRatio + 1 ) / ( 2.0 * densRatio + 1 );
}

template<typename Species>
Double TTransFlameSolver<Species>::GetZStoi( void )
{ 
	int				i;
	int				nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	int				*indexFuel = GetFuelIndexVec()->vec;
	int				indexOx;
	TInputDataPtr	input = GetInputData();
	Double			zStoi;
	Double			nuOx = fInputData->fOxIndex;
	Double			nuFuel;
	Double			*molarMass = fSpecies.GetMolarMass()->vec;
	Double			nu = 0.0;
	Double			sumYFuel = 0.0;
	Double			*Y2 = &fSolution->mat[kPrev][fFirstSpecies];
	Double			*Y1 = &fSolution->mat[fNGridPoints][fFirstSpecies];
	char			**names = fSpecies.GetNames();
	
	indexOx = fInputData->fOxIndex;
	
	nuOx = GetNu( fInputData->fGlobalReaction, names[indexOx] );

	if ( nuOx == -1 ) {
		fprintf( fOutFilePtr, "###error: oxidizer '%s' is not in the global reaction\n", names[indexOx] );
		exit(2);
	}
	
	for ( i = 0; i < GetNFuels(); ++i ) {
		nuFuel = GetNu( input->fGlobalReaction, names[indexFuel[i]] );
		if ( nuFuel == -1 ) {
			fprintf( fOutFilePtr, "###error: fuel '%s' is not in the global reaction\n", names[indexFuel[i]] );
			exit(2);
		}
		nu += nuFuel * molarMass[indexFuel[i]];
		sumYFuel += Y1[indexFuel[i]];
	}

	nu = nuOx * molarMass[indexOx] / nu;

	zStoi = 1.0 / ( 1.0 + nu * sumYFuel / Y2[indexOx] );

	if ( zStoi < 1.0e-10 || zStoi > 0.9 ) {
		fprintf( fOutFilePtr, "#warning: the value of ZStoi = %g indicates an error in the inital conditons\n"
				, zStoi );
	}
	else {
		fprintf( fOutFilePtr, "ZStoi = %g\n"
				, zStoi );
	}
	
	return zStoi;
}

template<typename Species>
Double TTransFlameSolver<Species>::GetDelQ( Double t, Double Z )
{
	int	i = 1;
	int	j = 1;
	int	zlen = fZIn->len;
	int	tlen = fTimeIn->len;
	Double	*zIn = fZIn->vec;
	Double	*tIn = fTimeIn->vec;
	Double	**delQIn = fDelQ->mat;
	Double	delQAtZt1, delQAtZt2;
	while ( j < zlen-1 && Z > zIn[j] ) ++j;
	if ( j == zlen-1 && Z > zIn[j]+1.0e-8 ) {
		fprintf( stderr, "Z = %g out of range[%g,%g]\n", Z, zIn[0]
								, zIn[zlen-1] );
	}
	
	while ( i < tlen-1 && t > tIn[i] ) ++i;
	if ( i == tlen-1 && t > tIn[i] ) {
		fprintf( stderr, "t = %g out of range[%g,%g]: linear extrapolation\n", Z, tIn[0]
								, tIn[tlen-1] );
	}
	
	if ( zIn[j] == zIn[j-1] ) {
		fprintf( stderr, "something's wrong\n" );
	}
	else {
		delQAtZt1 = Interpol( Z, delQIn[i-1][j-1], zIn[j-1], delQIn[i-1][j], zIn[j] );
		delQAtZt2 = Interpol( Z, delQIn[i][j-1], zIn[j-1], delQIn[i][j], zIn[j] );
	}
	
	if ( tIn[i] == tIn[i-1] ) {
		fprintf( stderr, "something else's wrong\n" );
		return 0.0;
	}
	else {
		return Interpol( t, delQAtZt1, tIn[i-1], delQAtZt2, tIn[i] );
	}
}

template<typename Species>
int	TTransFlameSolver<Species>::GetVariableIndex( const char *name )
{
	for ( int i = 0; i < fNOfEquations; ++ i ) {
		if ( strcmp( name, fVariableNames[i] ) == 0 ) {
			return i;
		}
	}

	return -1;
}

template<typename Species>
int	TTransFlameSolver<Species>::GetVariableIndex( const char *name, ConstStringArray array, int len )
{
	for ( int i = 0; i < len; ++ i ) {
		if ( strcmp( name, array[i] ) == 0 ) {
			return i;
		}
	}

	return -1;
}

/************************Rui101218***********************************************************************/


#undef WRITEGLOBALPRODRATE

template<typename Species>
void TTransFlameSolver<Species>::WriteFlameletFileInitialize( FILE *fp, const char *head, const char *tail )
{
	T0DPropertiesPtr	props = GetProperties();
	Double				*molarMass = GetSpecies()->GetMolarMass()->vec;
	int					i, k;
	int					nOfSpecies = GetSpecies()->GetNOfSpecies();
	int					nSpeciesIn = GetSpecies()->GetNSpeciesInSystem();
	int					nSootMoments;
	time_t				theDate;
	char				buffer[80];
	char				**names = GetSpecies()->GetNames();
	Flag				fpOpen = FALSE;
	Double				*x = fSolGrid->vec;
	Double				**Y = fMassFracsWork->mat;
	Double				**theMom;
	Double				*temp = fTempWork->vec;
	Double				theTime = GetCurrentTime();
	char				tmpName[127];
	VectorPtr 			ZBilgerVec = NewVector( fNGridPoints + 2 );
	Double				*ZBilger = &ZBilgerVec->vec[kNext];
	string              GasRadiationName = fProperties->GetRadiationName();

	static const Double	molarMassC = 12.01, 
						molarMassO = 16.0,
						molarMassH = 1.008;
	static Double		elemMassLastCInit = 10000.0;
	static Double		elemMassLastHInit = 10000.0;
	static Double		elemMassLastOInit = 10000.0;

	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		theMom = fSootMomentsWork->mat;
	}

	SetOutSolution();

	if ( fZREnd >= 0.9999 ) {
		elemMassLastCInit = GetElementMassFraction( Y[fNGridPoints], "C", molarMassC );
		elemMassLastHInit = GetElementMassFraction( Y[fNGridPoints], "H", molarMassH );
		elemMassLastOInit = GetElementMassFraction( Y[fNGridPoints], "O", molarMassO );
	}

	// if ( !fp ) {
	// 	fpOpen = TRUE;
	// 	fp = GetOutputFile( theTime, head, tail, FileType::kText );
	// }

/*************************************Replaced the above by Rui**********************/
//(1) Find the C_st
	Double	zR = Interpol( theTime, fZRStart, fTStart, fZREnd, fTEnd );
	Double	*mixfrac = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		mixfrac[k+1] = x[k] * zR;
		if ( mixfrac[k+1] < 1.0e-10 ) {
			mixfrac[k+1] = 0.0;
		}
	}
	int	CO2 = fSpecies.FindSpecies( "CO2" );
	int	CO = fSpecies.FindSpecies( "CO" );
	int	H2O = fSpecies.FindSpecies( "H2O" );
	int	H2 = fSpecies.FindSpecies( "H2" );
	int C2H2 = fSpecies.FindSpecies( "C2H2" );
	Double	*prog = New1DArray( fNGridPoints+2 );

	Double C_st = 0.0; 
	Double dH_st = 0.0;


	for ( k = 0; k <= fNGridPoints+1; ++k ) {
		prog[k] = Y[k-1][CO2] + Y[k-1][CO] + Y[k-1][H2O] + Y[k-1][H2];
	}
	for ( k = -1; k <= fNGridPoints; ++k ) {
		if ( mixfrac[k+1] > fZRef ) {
			C_st = prog[k]+(prog[k+1]-prog[k])*(fZRef-mixfrac[k])/(mixfrac[k+1]-mixfrac[k]);
			dH_st = abs( GetTotEnt(k-1,x,theTime) - GetTotEnt(temp[k-1], Y[k-1]) ) + 
						( abs( GetTotEnt(k,x,theTime) - GetTotEnt(temp[k],Y[k]) ) - abs( GetTotEnt(k-1,x,theTime) - GetTotEnt(temp[k-1],Y[k-1]) ) ) * 
							(fZRef-mixfrac[k])/(mixfrac[k+1]-mixfrac[k]);
			break;
		}
	}
//(2) Find the deltaH_max
	// char   tl[128];
	// if ( !fp ) {
	// 	fpOpen = TRUE;
	// 	sprintf( tl, "T%05g_dHmax%05g_Cst%05g_Chi%05g_t%05g", temp[LocationOfMax( fNGridPoints+2, &temp[kPrev] )-1], 
	// 		0.0, C_st, GetRefDissRate( theTime ), theTime*1000 );
	// 	fp = GetOutfile( tl, FileType::kText );
	// }

//(3) write the deltaH_st into title
		char   tl[128];
	if ( !fp ) {
		fpOpen = TRUE;
		sprintf( tl, "T%05g_dHst%05g_Cst%05g_Chi%05g_t%05g", temp[LocationOfMax( fNGridPoints+2, &temp[kPrev] )-1], 
			dH_st/1000.0, C_st, GetRefDissRate( theTime ), theTime*1000 );
		fp = GetOutfile( tl, FileType::kText );
	}
 /*************************************End**********************/

	// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"transient flamelet\"\n" );
	fprintf( fp, "mechanism = \"%s\"\n", fInputData->fReactionFile );
	fprintf( fp, "author = \"%s\"\n", GetAuthor() );
	time( &theDate );
	strcpy( buffer, ctime( &theDate ) );
	if ( buffer[strlen(buffer)-1] == '\n' )
		buffer[strlen(buffer)-1] = '\0';
	fprintf( fp, "date = \"%s\"\n\n", buffer );
	for ( i = 0; i < GetNFuels(); ++i ) {
		fprintf( fp, "fuel = \"%s\"\n", fVariableNames[fFirstSpecies+GetFuelIndex( i )] );
	}
	fprintf( fp, "time = %g [ms]\n", theTime * 1.0e3 );
	fprintf( fp, "pressure = %g [bar]\n", GetPressure( theTime ) / 1.0e5 );
	fprintf( fp, "chi_ref = %g [1/s]\n", GetRefDissRate( theTime ) );

	fprintf( fp, "Zst = %g [1/s]\n", fZRef );

	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( fNGridPoints+2
												, &temp[kPrev] )-1] );
	fprintf( fp, "ZofTmax = %g [K]\n", x[LocationOfMax( fNGridPoints+2
												, &temp[kPrev] )-1] );

	fprintf( fp, "DefectMax = %e [J/kg]\n",  0.0);

	// Double	zR = Interpol( theTime, fZRStart, fTStart, fZREnd, fTEnd );
	fprintf( fp, "ZR = %g\n", zR );

	int gp = fNGridPoints;
	while ( gp >= 0 && temp[gp-1] > temp[gp] ) --gp;
	fprintf( fp, "LocTmax = %g [K]\n", temp[gp] );

	Double	EIFuel = 0.0;
	for ( i = 0; i < GetNFuels(); ++i ) {
		EIFuel += ComputeEmissionIndex( theTime, GetFuelIndex( i ), NULL, NULL );
	}
	fprintf( fp, "EmissionIndexFuel = %g [kg/m^3s]\n", EIFuel );
	int	index = GetSpecies()->FindSpecies( "NO" );
	if ( index > -1 ) {
		fprintf( fp, "EmissionIndexNO = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( theTime, index, NULL, NULL ) / EIFuel );
	}
	
	index = GetSpecies()->FindSpecies( "NO2" );
	if ( index > -1 && !GetSpecies()->IsSteadyState(index)) {
		fprintf( fp, "EmissionIndexNO2 = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( theTime, index, NULL, NULL ) / EIFuel );
	}
	
	index = GetSpecies()->FindSpecies( "N2O" );
	if ( index > -1 && !GetSpecies()->IsSteadyState(index)) {
		fprintf( fp, "EmissionIndexN2O = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( theTime, index, NULL, NULL ) );
	}

	if ( fSoot ) {
		fprintf( fp, "EmissionIndexSoot = %g [kg/sm^3]\n"
				, ComputeEmissionIndexSoot( theTime, 1, NULL, NULL ) * 24 );
	}

	fprintf( fp, "CurrTimeStep = %g [s]\n", *fActualTimeStepSize );

	fprintf( fp, "\nFuelSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempFuelEnd );
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( fabs( Y[fNGridPoints][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[fNGridPoints][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "OxidizerSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempOxEnd );
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( fabs( Y[kPrev][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[kPrev][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
	fprintf( fp, "gridPoints = %d\n\n", fNGridPoints+2 );

	fprintf( fp, "body\n" );

	
	ZBilger[kPrev] = x[kPrev];
	ZBilger[fNGridPoints] = x[fNGridPoints];
	for ( k = 0; k < fNGridPoints; ++k ) {
		ZBilger[k] = ComputeZBilger( Y[k], Y[fNGridPoints], Y[kPrev] );
	}

// write solution
	// Double	*mixfrac = New1DArray( fNGridPoints+2 );

	// for ( k = -1; k <= fNGridPoints; ++k ) {
	// 	mixfrac[k+1] = x[k] * zR;
	// 	if ( mixfrac[k+1] < 1.0e-10 ) {
	// 		mixfrac[k+1] = 0.0;
	// 	}
	// }

	PrintFlameletVector( fNGridPoints+2, mixfrac, "Z", fp );
	PrintFlameletVector( fNGridPoints+2, &x[kPrev], "zeta", fp );

	Free1DArray( mixfrac );

	PrintFlameletVector( fNGridPoints+2, &temp[kPrev], "temperature [K]", fp );

	// write massfractions of species
	for ( i = 0; i < nOfSpecies; ++i ) {
		sprintf( tmpName, "massfraction-%s", names[i] );
		PrintFlameletVector( fNGridPoints+2, &Y[kPrev][i], tmpName
				, fp, fMassFracsWork->phys_rows );
	}

	if ( fPrintMolarFractions ) {

		MatrixPtr	moleMat = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
		Double		**molefracs = &moleMat->mat[kNext];
		for ( k = -1; k <= fNGridPoints; ++k )
		{
			fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
					, Y[k], molarMass, nSpeciesIn );
			for ( i = 0; i < nOfSpecies; ++i ) {
				molefracs[k][i] = Y[k][i] * props->GetMixMolarMass() / molarMass[i];
			}
		}
		for ( i = 0; i < nOfSpecies; ++i ) {
			sprintf( tmpName, "molarfraction-%s", names[i] );
			PrintFlameletVector( fNGridPoints+2, &molefracs[kPrev][i], tmpName
					, fp, moleMat->phys_rows );
		}
		DisposeMatrix( moleMat );
	}
	
	// write progvar 
	// int	CO2 = fSpecies.FindSpecies( "CO2" );
	// int	CO = fSpecies.FindSpecies( "CO" );
	// int	H2O = fSpecies.FindSpecies( "H2O" );
	// int	H2 = fSpecies.FindSpecies( "H2" );
	// Double	*prog = New1DArray( fNGridPoints+2 );
	sprintf( tmpName, "ProgVar" );
	// for ( k = 0; k <= fNGridPoints+1; ++k ) {
	// 	prog[k] = Y[k-1][CO2] + Y[k-1][CO] + Y[k-1][H2O] + Y[k-1][H2];
	// }
	PrintFlameletVector( fNGridPoints+2, prog, "ProgVar", fp );
	Free1DArray( prog );

	PrintFlameletVector( fNGridPoints+2, &ZBilger[kPrev], "ZBilger", fp );

// element mixture fraction
	Double	*ZElem = New1DArray( fNGridPoints+2 );
	Double	elemMassFirst;
	
	
// C
	elemMassFirst = GetElementMassFraction( Y[kPrev], "C", molarMassC );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZElem[k+1] = ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirst ) / MAX( elemMassLastCInit - elemMassFirst, 1.0e-30 );
	}
	PrintFlameletVector( fNGridPoints+2, ZElem, "CMassFrac", fp );
	
// H
	elemMassFirst = GetElementMassFraction( Y[kPrev], "H", molarMassH );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZElem[k+1] = ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirst ) / MAX( elemMassLastHInit - elemMassFirst, 1.0e-30 );
	}
	PrintFlameletVector( fNGridPoints+2, ZElem, "HMassFrac", fp );

// O
	elemMassFirst = GetElementMassFraction( Y[kPrev], "O", molarMassO );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZElem[k+1] = ( GetElementMassFraction( Y[k], "O", molarMassO ) - elemMassFirst ) / MAX( elemMassLastOInit - elemMassFirst, 1.0e-30 );
	}
	PrintFlameletVector( fNGridPoints+2, ZElem, "OMassFrac", fp );

	Free1DArray( ZElem );

// Mixture fraction following Barlows definition
	Double	*ZBarlow = New1DArray( fNGridPoints+2 );
	Double	elemMassFirstH = GetElementMassFraction( Y[kPrev], "H", molarMassH );
	Double	elemMassFirstC = GetElementMassFraction( Y[kPrev], "C", molarMassC );
	Double	denomBarlow = 0.5 * ( elemMassLastHInit - elemMassFirstH ) / molarMassH
				+ 2.0 * ( elemMassLastCInit - elemMassFirstC ) / molarMassC;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZBarlow[k+1] = ( 0.5 * ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirstH ) / molarMassH
						+ 2.0 * ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirstC ) / molarMassC )
						/ denomBarlow;
	}

	PrintFlameletVector( fNGridPoints+2, ZBarlow, "ZBarlow", fp );

	Free1DArray( ZBarlow );


// density
	Double	*density = New1DArray( fNGridPoints+2 );
	Double	*locchi = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
				, Y[k], molarMass, nSpeciesIn );
		density[k+1] = GetPressure( theTime ) * props->GetMixMolarMass() / ( RGAS * temp[k] );
		locchi[k+1] = GetDissRate( theTime, x[k] );
	}

	PrintFlameletVector( fNGridPoints+2, density, "density [kg/m^3]", fp );
	PrintFlameletVector( fNGridPoints+2, locchi, "chi [1/s]", fp );

	Free1DArray( locchi );

	Double	*dummy = New1DArray( fNGridPoints+2 );

	for ( k = -1; k <= fNGridPoints; ++k ) {
		dummy[k+1] = k+1;
	}

	Free1DArray( dummy );

		Double	*sootrad;
		Double	*gasrad;
		if ( fProperties->GetRadiation() ) {
			sootrad = New1DArray( fNGridPoints+2 );
			gasrad = New1DArray( fNGridPoints+2 );
		}
		Double	*summh = New1DArray( fNGridPoints+2 );
		Double	*prodRate = fSpecies.GetProductionRate()->vec;
		Double	**prodRateGrid = fProdRate->mat;
		Double	*enth = fSpecies.GetEnthalpy()->vec;
		Double	press = GetPressure( theTime );
		Double	*cpMix = fHeatCpMix->vec;
		Double	*mu = fViscosity->vec;
		Double	*mixMolarMass = fMolarMassMix->vec;
		Double	*lambdaOverCp = fLambdaOverCpMix->vec;
		
		for ( k = -1; k <= fNGridPoints; ++k ) {
#ifdef DELTAINEW
			UpdateThermoProps( k, Y[k], temp[k], press, fProperties->GetDensityRef()
										, kDensFromPress, (fSoot) ? theMom[k] : NULL );
#else 
			T0DFlame<Species>::UpdateThermoProps( Y[k], temp[k], press, fProperties->GetDensityRef()
										, kDensFromPress, (fSoot) ? theMom[k] : NULL );
#endif
		
			summh[k+1] = 0.0;
			for ( i = 0; i < nSpeciesIn; ++i ) {
				summh[k+1] += enth[i] * prodRate[i];
				prodRateGrid[k][i] = prodRate[i];
			}
			cpMix[k] = fProperties->GetMixHeatCapacity();
			mu[k] = fProperties->GetMixViscosity();
			lambdaOverCp[k] = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();
			mixMolarMass[k] = fProperties->GetMixMolarMass();

			if ( fProperties->GetRadiation() ) {
				// if (GasRadiationName == "Thin" ) {
				// 	gasrad[k+1] = fProperties->GetRadiation()->GetRadiation();
				// }
				// else if (GasRadiationName == "WSGG" || GasRadiationName == "Grey" || GasRadiationName == "SNB" ) {
				// 	if ( k == -1 || k == fNGridPoints ) {
				// 		gasrad[k+1] = 0.0;
				// 	}
				// 	else {
				// 		gasrad[k+1] = - fProperties->GetRadiation()->dqrF[k];   // Note dqrF starts from 0 and doesn't have B.C.
				// 	}
				// }
				// else{
				// 	cerr << "Error: gas radiation model: " << GasRadiationName << " is not found!!!" << NEWL;
				// }

				if ( fSoot ) {
#ifdef ROSSELAND
				  if (k==-1||k==fNGridPoints) {
					sootrad[k+1] = 0.0;
				  }
				  else {
					Double rfM = density[k] * GetDissRate( theTime, x[k-1] ) * lambdaOverCp[k-1];
					Double rfC = density[k+1] * GetDissRate( theTime, x[k] ) * lambdaOverCp[k];
					Double rfP = density[k+2] * GetDissRate( theTime, x[k+1] ) * lambdaOverCp[k+1];

					sootrad[k+1] = -GetRosseRadiation( k, &temp[k], &theMom[k], rfM, rfC, rfP );
				  }
#else

					sootrad[k+1] = -fSoot->GetSootRadiation( temp[k], theMom[k] );
#endif
				}
			}
		}
	
		if ( fProperties->GetRadiation() ) {
			if ( fSoot ) {
				PrintFlameletVector( fNGridPoints+2, sootrad, "SootRadiation [J/m^3s]", fp );
			}
			PrintFlameletVector( fNGridPoints+2, gasrad, "GasRadiation [J/m^3s]", fp );
		}
		PrintFlameletVector( fNGridPoints+2, summh, "SumMH [J/m^3s]", fp );
		PrintFlameletVector( fNGridPoints+2, &cpMix[kPrev], "Cp", fp );
		PrintFlameletVector( fNGridPoints+2, &mu[kPrev], "mu", fp );
		PrintFlameletVector( fNGridPoints+2, &lambdaOverCp[kPrev], "lambdaOverCp", fp );
		PrintFlameletVector( fNGridPoints+2, &mixMolarMass[kPrev], "W", fp );
	
	
//  write source term ProdVar
	fprintf( fp, "ProdRateProgVar [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][CO2] + prodRateGrid[k][CO] + prodRateGrid[k][H2O] + prodRateGrid[k][H2] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );
	
//  write source term CO2
	fprintf( fp, "ProdRateCO2 [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][CO2] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term H2O
	fprintf( fp, "ProdRateH2O [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][H2O] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term H2
	fprintf( fp, "ProdRateH2 [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][H2] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term CO
	fprintf( fp, "ProdRateCO [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][CO] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );


		Free1DArray( summh );
		if ( fProperties->GetRadiation() ) {
			Free1DArray( gasrad );
			Free1DArray( sootrad );
		}
		Free1DArray( density );
	
	if ( fSoot ) {
		for ( i = 0; i < nSootMoments; ++i ) {
			sprintf( tmpName, "M%d", i );
			PrintFlameletVector( fNGridPoints+2, &theMom[kPrev][i], tmpName
					, fp, fSootMomentsWork->phys_rows );
		}
		
		WriteSootInfo( theTime, temp, Y, GetPressure( theTime ), theMom, fp );		
	}

	Double	*totEnt = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k )
	{
		totEnt[k+1] = GetTotEnt( k, x, theTime );
	}
	PrintFlameletVector( fNGridPoints+2, totEnt, "TotalEnthalpy [J/kg]", fp );

	for ( k = -1; k <= fNGridPoints; ++k )
	{
		totEnt[k+1] = GetTotEnt( temp[k], Y[k] );
	}
	PrintFlameletVector( fNGridPoints+2, totEnt, "TotalEnthalpy2 [J/kg]", fp );

	Free1DArray( totEnt );

	PrintFlameletVector( fNGridPoints+2, &fMonFct[kPrev], "MonitorFunc", fp );

	// Production Rate of NO
	index = GetSpecies()->FindSpecies( "NO" );
	if ( index > -1 ) {
		PrintProdRate( theTime, index, fp );
	}
	// Production Rate of NO
	index = GetSpecies()->FindSpecies( "CH4" );
	if ( index > -1 ) {
		PrintProdRate( theTime, index, fp );
	}
	
#ifdef WRITEGLOBALPRODRATE
	PrintProdRateGlobalReac( theTime );	
#endif
	
	fprintf( fp, "trailer\n" );
	
#ifdef LEWISCHANGE
	Double					*LeOrig = fSpecies.GetLewisNumber()->vec;
	Double					*Le = new Double[nSpeciesIn];
	Double					LeFunc;
	Double					switchTime = LEWISSWITCHTIME;
	Double					switchPeriod = LEWISSWITCHPERIOD;
	const Double			Pi = 4.0 * atan( 1.0 );
	Double					omega = 2.0 * Pi / switchPeriod;
	Double					tStar = theTime - ( switchTime - 0.5 * switchPeriod );

	if ( tStar >= 0.0 ) {
		if ( tStar < switchPeriod ) {
			LeFunc = 0.5 * ( cos( omega * tStar ) + 1.0 );
			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
				Le[iLe] = LeFunc * ( LeOrig[iLe] - 1.0 ) + 1.0;
			}
		}
		else {
			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
				Le[iLe] = 1.0;
			}
		}
	}
	else {
		for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
			Le[iLe] = LeOrig[iLe];
		}
	}
#else
	Double					*Le = fSpecies.GetLewisNumber()->vec;
#endif
	for ( i = 0; i < nSpeciesIn; ++i ) {
		fprintf( fp, "%s\t%g\n", names[i], Le[i] );
	}
	if ( fpOpen ) {
		fclose( fp );
	}
#ifdef LEWISCHANGE
	delete Le;
#endif
	DisposeVector( ZBilgerVec );
}
/************************Rui101218***********************************************************************/


#undef WRITEGLOBALPRODRATE

template<typename Species>
void TTransFlameSolver<Species>::WriteFlameletFile( FILE *fp, const char *head, const char *tail )
{
	T0DPropertiesPtr	props = GetProperties();
	Double				*molarMass = GetSpecies()->GetMolarMass()->vec;
	int					i, k;
	int					nOfSpecies = GetSpecies()->GetNOfSpecies();
	int					nSpeciesIn = GetSpecies()->GetNSpeciesInSystem();
	int					nSootMoments;
	time_t				theDate;
	char				buffer[80];
	char				**names = GetSpecies()->GetNames();
	Flag				fpOpen = FALSE;
	Double				*x = fSolGrid->vec;
	Double				**Y = fMassFracsWork->mat;
	Double				**theMom;
	Double				*temp = fTempWork->vec;
	Double				theTime = GetCurrentTime();
	char				tmpName[127];
	VectorPtr 			ZBilgerVec = NewVector( fNGridPoints + 2 );
	Double				*ZBilger = &ZBilgerVec->vec[kNext];
	string              GasRadiationName = fProperties->GetRadiationName();

	static const Double	molarMassC = 12.01, 
						molarMassO = 16.0,
						molarMassH = 1.008;
	static Double		elemMassLastCInit = 10000.0;
	static Double		elemMassLastHInit = 10000.0;
	static Double		elemMassLastOInit = 10000.0;

	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		theMom = fSootMomentsWork->mat;
	}

	SetOutSolution();

	if ( fZREnd >= 0.9999 ) {
		elemMassLastCInit = GetElementMassFraction( Y[fNGridPoints], "C", molarMassC );
		elemMassLastHInit = GetElementMassFraction( Y[fNGridPoints], "H", molarMassH );
		elemMassLastOInit = GetElementMassFraction( Y[fNGridPoints], "O", molarMassO );
	}

	// if ( !fp ) {
	// 	fpOpen = TRUE;
	// 	fp = GetOutputFile( theTime, head, tail, FileType::kText );
	// }
/*************************************Replaced the above by Rui**********************/
//(1) Find the C_st and dH_st
	Double	zR = Interpol( theTime, fZRStart, fTStart, fZREnd, fTEnd );
	Double	*mixfrac = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		mixfrac[k+1] = x[k] * zR;
		if ( mixfrac[k+1] < 1.0e-10 ) {
			mixfrac[k+1] = 0.0;
		}
	}
	int	CO2 = fSpecies.FindSpecies( "CO2" );
	int	CO = fSpecies.FindSpecies( "CO" );
	int	H2O = fSpecies.FindSpecies( "H2O" );
	int	H2 = fSpecies.FindSpecies( "H2" );
	int C2H2 = fSpecies.FindSpecies( "C2H2" );
	Double	*prog = New1DArray( fNGridPoints+2 );

	// sprintf( tmpName, "ProgVar" );
	Double C_st = 0.0; 
	Double dH_st = 0.0;

	for ( k = 0; k <= fNGridPoints+1; ++k ) {
		prog[k] = Y[k-1][CO2] + Y[k-1][CO] + Y[k-1][H2O] + Y[k-1][H2];
	}

	for ( k = -1; k <= fNGridPoints; ++k ) {
		if ( mixfrac[k+1] > fZRef ) {

			C_st = prog[k]+(prog[k+1]-prog[k])*(fZRef-mixfrac[k])/(mixfrac[k+1]-mixfrac[k]);

			dH_st = abs( GetTotEnt(k-1,x,theTime) - GetTotEnt(temp[k-1], Y[k-1]) ) + 
						( abs( GetTotEnt(k,x,theTime) - GetTotEnt(temp[k],Y[k]) ) - abs( GetTotEnt(k-1,x,theTime) - GetTotEnt(temp[k-1],Y[k-1]) ) ) * 
							(fZRef-mixfrac[k])/(mixfrac[k+1]-mixfrac[k]);

			break;
		}
	}

//(2) Find the deltaH_max
	// Double DefectMax = 0.0;        
	// Double DefectCurr = 0.0;
	// char   tl[128];
	// for ( k = -1; k <= fNGridPoints; ++k ){
	// 	DefectCurr = abs( GetTotEnt( k, x, theTime ) - GetTotEnt( temp[k], Y[k] )) ;
	// 	if (  DefectCurr > DefectMax ){
	// 		DefectMax = DefectCurr;
	// 	}
	// }
	// if ( !fp ) {
	// 	fpOpen = TRUE;
	// 	sprintf( tl, "T%05g_dHmax%05g_Cst%05g_Chi%05g_t%05g", temp[LocationOfMax( fNGridPoints+2, &temp[kPrev] )-1], 
	// 		DefectMax/1000, C_st, GetRefDissRate( theTime ), theTime*1000 );
	// 	fp = GetOutfile( tl, FileType::kText );
	// }

//(3) Write deltaH_st in file title , note the unit is kJ/kg
	char   tl[128];

	if ( !fp ) {
		fpOpen = TRUE;
		sprintf( tl, "T%05g_dHst%05g_Cst%05g_Chi%05g_t%05g", temp[LocationOfMax( fNGridPoints+2, &temp[kPrev] )-1], 
			dH_st/1000, C_st, GetRefDissRate( theTime ), theTime*1000 );
		fp = GetOutfile( tl, FileType::kText );
	}
 /*************************************End**********************/

	// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"transient flamelet\"\n" );
	fprintf( fp, "mechanism = \"%s\"\n", fInputData->fReactionFile );
	fprintf( fp, "author = \"%s\"\n", GetAuthor() );
	time( &theDate );
	strcpy( buffer, ctime( &theDate ) );
	if ( buffer[strlen(buffer)-1] == '\n' )
		buffer[strlen(buffer)-1] = '\0';
	fprintf( fp, "date = \"%s\"\n\n", buffer );
	for ( i = 0; i < GetNFuels(); ++i ) {
		fprintf( fp, "fuel = \"%s\"\n", fVariableNames[fFirstSpecies+GetFuelIndex( i )] );
	}
	fprintf( fp, "time = %g [ms]\n", theTime * 1.0e3 );
	fprintf( fp, "pressure = %g [bar]\n", GetPressure( theTime ) / 1.0e5 );
	fprintf( fp, "chi_ref = %g [1/s]\n", GetRefDissRate( theTime ) );

	fprintf( fp, "Zst = %g [1/s]\n", fZRef );

	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( fNGridPoints+2
												, &temp[kPrev] )-1] );

	fprintf( fp, "ZofTmax = %g [K]\n", x[LocationOfMax( fNGridPoints+2
												, &temp[kPrev] )-1] );

	fprintf( fp, "Cst = %e [-]\n",  C_st);
	fprintf( fp, "dHst = %e [J/kg]\n",  dH_st);

	// Double	zR = Interpol( theTime, fZRStart, fTStart, fZREnd, fTEnd );
	fprintf( fp, "ZR = %g\n", zR );

	int gp = fNGridPoints;
	while ( gp >= 0 && temp[gp-1] > temp[gp] ) --gp;
	fprintf( fp, "LocTmax = %g [K]\n", temp[gp] );


/*********************Commmented out by Rui 12/18/2018*********************************************/
	// Double	EIFuel = 0.0;
	// for ( i = 0; i < GetNFuels(); ++i ) {
	// 	EIFuel += ComputeEmissionIndex( theTime, GetFuelIndex( i ), NULL, NULL );
	// }
	// fprintf( fp, "EmissionIndexFuel = %g [kg/m^3s]\n", EIFuel );
	// int	index = GetSpecies()->FindSpecies( "NO" );
	// if ( index > -1 ) {
	// 	fprintf( fp, "EmissionIndexNO = %g [kg/sm^3]\n"
	// 			, ComputeEmissionIndex( theTime, index, NULL, NULL ) / EIFuel );
	// }
	
	// index = GetSpecies()->FindSpecies( "NO2" );
	// if ( index > -1 && !GetSpecies()->IsSteadyState(index)) {
	// 	fprintf( fp, "EmissionIndexNO2 = %g [kg/sm^3]\n"
	// 			, ComputeEmissionIndex( theTime, index, NULL, NULL ) / EIFuel );
	// }
	
	// index = GetSpecies()->FindSpecies( "N2O" );
	// if ( index > -1 && !GetSpecies()->IsSteadyState(index)) {
	// 	fprintf( fp, "EmissionIndexN2O = %g [kg/sm^3]\n"
	// 			, ComputeEmissionIndex( theTime, index, NULL, NULL ) );
	// }

	// if ( fSoot ) {
	// 	fprintf( fp, "EmissionIndexSoot = %g [kg/sm^3]\n"
	// 			, ComputeEmissionIndexSoot( theTime, 1, NULL, NULL ) * 24 );
	// }
/*********************Commmented out by Rui 12/18/2018*********************************************/



	fprintf( fp, "CurrTimeStep = %g [s]\n", *fActualTimeStepSize );

	fprintf( fp, "\nFuelSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempFuelEnd );

	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( fabs( Y[fNGridPoints][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[fNGridPoints][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "OxidizerSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempOxEnd );
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( fabs( Y[kPrev][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[kPrev][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
	fprintf( fp, "gridPoints = %d\n\n", fNGridPoints+2 );

	fprintf( fp, "body\n" );

	
	ZBilger[kPrev] = x[kPrev];
	ZBilger[fNGridPoints] = x[fNGridPoints];
	for ( k = 0; k < fNGridPoints; ++k ) {
		ZBilger[k] = ComputeZBilger( Y[k], Y[fNGridPoints], Y[kPrev] );
	}

// write solution
	// Double	*mixfrac = New1DArray( fNGridPoints+2 );

	// for ( k = -1; k <= fNGridPoints; ++k ) {
	// 	mixfrac[k+1] = x[k] * zR;
	// 	if ( mixfrac[k+1] < 1.0e-10 ) {
	// 		mixfrac[k+1] = 0.0;
	// 	}
	// }

	PrintFlameletVector( fNGridPoints+2, mixfrac, "Z", fp );
	PrintFlameletVector( fNGridPoints+2, &x[kPrev], "zeta", fp );
	Free1DArray( mixfrac );
	PrintFlameletVector( fNGridPoints+2, &temp[kPrev], "temperature [K]", fp );

// write massfractions of species
	for ( i = 0; i < nOfSpecies; ++i ) {
		sprintf( tmpName, "massfraction-%s", names[i] );
		// if (names[i] == "CO2" || names[i] == "H2O"){
		PrintFlameletVector( fNGridPoints+2, &Y[kPrev][i], tmpName, fp, fMassFracsWork->phys_rows );
		// }
	}

// write molar fractions of species
	if ( fPrintMolarFractions ) {

		MatrixPtr	moleMat = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
		Double		**molefracs = &moleMat->mat[kNext];
		for ( k = -1; k <= fNGridPoints; ++k )
		{
			fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
					, Y[k], molarMass, nSpeciesIn );
			for ( i = 0; i < nOfSpecies; ++i ) {
				molefracs[k][i] = Y[k][i] * props->GetMixMolarMass() / molarMass[i];
			}
		}
		for ( i = 0; i < nOfSpecies; ++i ) {
			sprintf( tmpName, "molarfraction-%s", names[i] );
			PrintFlameletVector( fNGridPoints+2, &molefracs[kPrev][i], tmpName
					, fp, moleMat->phys_rows );
		}
		DisposeMatrix( moleMat );
	}
	
// write progvar 
	// int	CO2 = fSpecies.FindSpecies( "CO2" );
	// int	CO = fSpecies.FindSpecies( "CO" );
	// int	H2O = fSpecies.FindSpecies( "H2O" );
	// int	H2 = fSpecies.FindSpecies( "H2" );
	// int C2H2 = fSpecies.FindSpecies( "C2H2" );


	// Double	*prog = New1DArray( fNGridPoints+2 );
	sprintf( tmpName, "ProgVar" );
	// for ( k = 0; k <= fNGridPoints+1; ++k ) {
	// 	prog[k] = Y[k-1][CO2] + Y[k-1][CO] + Y[k-1][H2O] + Y[k-1][H2];
	// }

	PrintFlameletVector( fNGridPoints+2, prog, "ProgVar", fp );
	Free1DArray( prog );

	PrintFlameletVector( fNGridPoints+2, &ZBilger[kPrev], "ZBilger", fp );

// element mixture fraction
	// Double	*ZElem = New1DArray( fNGridPoints+2 );
	// Double	elemMassFirst;
	
/*********************Commmented out by Rui 12/18/2018*********************************************/

// // C
// 	elemMassFirst = GetElementMassFraction( Y[kPrev], "C", molarMassC );
// 	for ( k = -1; k <= fNGridPoints; ++k ) {
// 		ZElem[k+1] = ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirst ) / MAX( elemMassLastCInit - elemMassFirst, 1.0e-30 );
// 	}
// 	PrintFlameletVector( fNGridPoints+2, ZElem, "CMassFrac", fp );
	
// // H
// 	elemMassFirst = GetElementMassFraction( Y[kPrev], "H", molarMassH );
// 	for ( k = -1; k <= fNGridPoints; ++k ) {
// 		ZElem[k+1] = ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirst ) / MAX( elemMassLastHInit - elemMassFirst, 1.0e-30 );
// 	}
// 	PrintFlameletVector( fNGridPoints+2, ZElem, "HMassFrac", fp );

// // O
// 	elemMassFirst = GetElementMassFraction( Y[kPrev], "O", molarMassO );
// 	for ( k = -1; k <= fNGridPoints; ++k ) {
// 		ZElem[k+1] = ( GetElementMassFraction( Y[k], "O", molarMassO ) - elemMassFirst ) / MAX( elemMassLastOInit - elemMassFirst, 1.0e-30 );
// 	}
// 	PrintFlameletVector( fNGridPoints+2, ZElem, "OMassFrac", fp );

// 	Free1DArray( ZElem );

// Mixture fraction following Barlows definition
	// Double	*ZBarlow = New1DArray( fNGridPoints+2 );
	// Double	elemMassFirstH = GetElementMassFraction( Y[kPrev], "H", molarMassH );
	// Double	elemMassFirstC = GetElementMassFraction( Y[kPrev], "C", molarMassC );
	// Double	denomBarlow = 0.5 * ( elemMassLastHInit - elemMassFirstH ) / molarMassH
	// 			+ 2.0 * ( elemMassLastCInit - elemMassFirstC ) / molarMassC;
	// for ( k = -1; k <= fNGridPoints; ++k ) {
	// 	ZBarlow[k+1] = ( 0.5 * ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirstH ) / molarMassH
	// 					+ 2.0 * ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirstC ) / molarMassC )
	// 					/ denomBarlow;
	// }

	// PrintFlameletVector( fNGridPoints+2, ZBarlow, "ZBarlow", fp );

	// Free1DArray( ZBarlow );
/*********************Commmented out by Rui 12/18/2018*********************************************/



// density
	Double	*density = New1DArray( fNGridPoints+2 );
	Double	*locchi = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
				, Y[k], molarMass, nSpeciesIn );
		density[k+1] = GetPressure( theTime ) * props->GetMixMolarMass() / ( RGAS * temp[k] );
		locchi[k+1] = GetDissRate( theTime, x[k] );
	}

	PrintFlameletVector( fNGridPoints+2, density, "density [kg/m^3]", fp );
	PrintFlameletVector( fNGridPoints+2, locchi, "chi [1/s]", fp );

	Free1DArray( locchi );

	Double	*dummy = New1DArray( fNGridPoints+2 );

	for ( k = -1; k <= fNGridPoints; ++k ) {
		dummy[k+1] = k+1;
	}

	Free1DArray( dummy );

	Double	*sootrad;
	Double	*gasrad;
	Double	*Kappa_rad;
	Double	*Emission_rad;


	// For 5 bands WSGG to consider TRI
	Double	*Kappa_rad1;
	Double	*Kappa_rad2;
	Double	*Kappa_rad3;
	Double	*Kappa_rad4;
	Double	*Emission_rad1;
	Double	*Emission_rad2;
	Double	*Emission_rad3;
	Double	*Emission_rad4;

	if ( fProperties->GetRadiation() ) {
		sootrad = New1DArray( fNGridPoints+2 );
		gasrad = New1DArray( fNGridPoints+2 );
		Kappa_rad = New1DArray( fNGridPoints+2 );
		Emission_rad = New1DArray( fNGridPoints+2 );
	}

	if ( GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar" ){
		Kappa_rad1 = New1DArray( fNGridPoints+2 );
		Kappa_rad2 = New1DArray( fNGridPoints+2 );
		Kappa_rad3 = New1DArray( fNGridPoints+2 );
		Kappa_rad4 = New1DArray( fNGridPoints+2 );
		Emission_rad1 = New1DArray( fNGridPoints+2 );
		Emission_rad2 = New1DArray( fNGridPoints+2 );
		Emission_rad3 = New1DArray( fNGridPoints+2 );
		Emission_rad4 = New1DArray( fNGridPoints+2 );
	}


	Double	*summh = New1DArray( fNGridPoints+2 );

	Double	*prodRate = fSpecies.GetProductionRate()->vec;
	Double	**prodRateGrid = fProdRate->mat;
	Double	*enth = fSpecies.GetEnthalpy()->vec;
	Double	press = GetPressure( theTime );
	Double	*cpMix = fHeatCpMix->vec;
	Double	*mu = fViscosity->vec;
	Double	*mixMolarMass = fMolarMassMix->vec;
	Double	*lambdaOverCp = fLambdaOverCpMix->vec;
	
	Double	**SpecSourceTerm = fSpecSource->mat;    // Rui

	int     currSpeciesIndex;

		
	for ( k = -1; k <= fNGridPoints; ++k ) {

#ifdef DELTAINEW
			UpdateThermoProps( k, Y[k], temp[k], press, fProperties->GetDensityRef()
										, kDensFromPress, (fSoot) ? theMom[k] : NULL );
#else 
			T0DFlame<Species>::UpdateThermoProps( Y[k], temp[k], press, fProperties->GetDensityRef()
										, kDensFromPress, (fSoot) ? theMom[k] : NULL );
#endif
			summh[k+1] = 0.0;
			for ( i = 0; i < nSpeciesIn; ++i ) {
				summh[k+1] += enth[i] * prodRate[i];
				prodRateGrid[k][i] = prodRate[i];
				SpecSourceTerm[k][i] = prodRate[i] * enth[i];   // Rui
			}

			cpMix[k] = fProperties->GetMixHeatCapacity();
			mu[k] = fProperties->GetMixViscosity();
			lambdaOverCp[k] = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();
			mixMolarMass[k] = fProperties->GetMixMolarMass();

			if ( fProperties->GetRadiation() ) {
				if (GasRadiationName == "Thin" ) {
					fProperties->GetRadiation()->SetRadiation( temp[k], Y[k], molarMass, density[k] );
					gasrad[k+1] = fProperties->GetRadiation()->GetRadiation();
				}

				else if (GasRadiationName == "WSGG" || GasRadiationName == "Grey" || GasRadiationName == "SNB" || 
					GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar") {

					if ( k == -1 || k == fNGridPoints ) {

						gasrad[k+1] = 0.0;
						Kappa_rad[k+1] = 0.0;
						Emission_rad[k+1] = 0.0;

						if (GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar") {

							Kappa_rad1[k+1] = 0.0;
							Kappa_rad2[k+1] = 0.0;
							Kappa_rad3[k+1] = 0.0;
							Kappa_rad4[k+1] = 0.0;
							Emission_rad1[k+1] = 0.0;
							Emission_rad2[k+1] = 0.0;
							Emission_rad3[k+1] = 0.0;
							Emission_rad4[k+1] = 0.0;
						}

					
					}

					else {

						gasrad[k+1] = - fProperties->GetRadiation()->dqrF[k] * fRadLossPercent; // Note dqrF starts from 0 and doesn't have B.C.
						Kappa_rad[k+1] = fProperties->GetRadiation()->kappa[k];
						Emission_rad[k+1] =  Kappa_rad[k+1]*5.670367E-8*pow(temp[k],4)/M_PI;


						if (GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar"){

							Kappa_rad1[k+1] = fProperties->GetRadiation()->kappa1[0][k];
							Kappa_rad2[k+1] = fProperties->GetRadiation()->kappa2[0][k];
							Kappa_rad3[k+1] = fProperties->GetRadiation()->kappa3[0][k];
							Kappa_rad4[k+1] = fProperties->GetRadiation()->kappa4[0][k];
							Emission_rad1[k+1] =  fProperties->GetRadiation()->kappa1[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
							Emission_rad2[k+1] =  fProperties->GetRadiation()->kappa2[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
							Emission_rad3[k+1] =  fProperties->GetRadiation()->kappa3[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
							Emission_rad4[k+1] =  fProperties->GetRadiation()->kappa4[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
						}
					}
				}

				else{
					cerr << "Error: gas radiation model: " << GasRadiationName << " is not found!!!" << NEWL;
				}

				if ( fSoot ) {
#ifdef ROSSELAND
				  if (k==-1||k==fNGridPoints) {
					sootrad[k+1] = 0.0;
				  }
				  else {
					Double rfM = density[k] * GetDissRate( theTime, x[k-1] ) * lambdaOverCp[k-1];
					Double rfC = density[k+1] * GetDissRate( theTime, x[k] ) * lambdaOverCp[k];
					Double rfP = density[k+2] * GetDissRate( theTime, x[k+1] ) * lambdaOverCp[k+1];

					sootrad[k+1] = -GetRosseRadiation( k, &temp[k], &theMom[k], rfM, rfC, rfP );
				  }
#else

					sootrad[k+1] = -fSoot->GetSootRadiation( temp[k], theMom[k] );
#endif
				}
			}

	}
	
		if ( fProperties->GetRadiation() ) {
			if ( fSoot ) {
				PrintFlameletVector( fNGridPoints+2, sootrad, "SootRadiation [J/m^3s]", fp );
			}
			PrintFlameletVector( fNGridPoints+2, gasrad, "GasRadiation [J/m^3s]", fp );
			PrintFlameletVector( fNGridPoints+2, Kappa_rad, "Kappa_rad [1/m]", fp );
			PrintFlameletVector( fNGridPoints+2, Emission_rad, "Emission_rad [J/s/m^3/sr]", fp );

			if (GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar"){

				PrintFlameletVector( fNGridPoints+2, Kappa_rad1, "Kappa_rad1 [1/m]", fp );
				PrintFlameletVector( fNGridPoints+2, Kappa_rad2, "Kappa_rad2 [1/m]", fp );
				PrintFlameletVector( fNGridPoints+2, Kappa_rad3, "Kappa_rad3 [1/m]", fp );
				PrintFlameletVector( fNGridPoints+2, Kappa_rad4, "Kappa_rad4 [1/m]", fp );
				PrintFlameletVector( fNGridPoints+2, Emission_rad1, "Emission_rad1 [J/s/m^3/sr]", fp );
				PrintFlameletVector( fNGridPoints+2, Emission_rad2, "Emission_rad2 [J/s/m^3/sr]", fp );
				PrintFlameletVector( fNGridPoints+2, Emission_rad3, "Emission_rad3 [J/s/m^3/sr]", fp );
				PrintFlameletVector( fNGridPoints+2, Emission_rad4, "Emission_rad4 [J/s/m^3/sr]", fp );

			}
		}

		PrintFlameletVector( fNGridPoints+2, summh, "SumMH [J/m^3s]", fp );
		PrintFlameletVector( fNGridPoints+2, &cpMix[kPrev], "Cp", fp );
		PrintFlameletVector( fNGridPoints+2, &mu[kPrev], "mu", fp );
		PrintFlameletVector( fNGridPoints+2, &lambdaOverCp[kPrev], "lambdaOverCp", fp );
		PrintFlameletVector( fNGridPoints+2, &mixMolarMass[kPrev], "W", fp );


 // Source term (Prodrate * enth) of all species         by Rui 01/12/2019
	// for ( i = 0; i < nOfSpecies; ++i ) {
	// 	fprintf( fp, "SpecSourceTerm-%s\n", names[i] );
	// 	fprintf( fp, "\t%-.6e", 0.0 );
	// 	for ( k = 0; k < fNGridPoints; ++k ) {
	// 		fprintf( fp, "\t%-.6e", SpecSourceTerm[k][i]);
	// 		if ( (k+2) % 5 == 0 ) {
	// 			fprintf( fp, "\n" );
	// 		}
	// 	}
	// 	fprintf( fp, "\t%-.6e\n", 0.0 );
	// }

 // Source term (Prodrate * enth) of C2H2, CO2, H2, H2O         by Rui 01/12/2019
	// int index = -1 ;

	// index = GetSpecies()->FindSpecies( "C2H2" );
	// fprintf( fp, "SpecSourceTerm-C2H2\n");
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", SpecSourceTerm[k][index] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );

	// index = GetSpecies()->FindSpecies( "CO2" );
	// fprintf( fp, "SpecSourceTerm-CO2\n");
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", SpecSourceTerm[k][index] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );


	// index = GetSpecies()->FindSpecies( "H2" );
	// fprintf( fp, "SpecSourceTerm-H2\n");
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", SpecSourceTerm[k][index] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );

	// index = GetSpecies()->FindSpecies( "H2O" );
	// fprintf( fp, "SpecSourceTerm-H2O\n");
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", SpecSourceTerm[k][index] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );

	// index = GetSpecies()->FindSpecies( "CO" );
	// fprintf( fp, "SpecSourceTerm-CO\n");
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", SpecSourceTerm[k][index] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );

 // write source term ProdVar
	// fprintf( fp, "ProdRateProgVar [kg/m^3s]\n" );
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", prodRateGrid[k][CO2] + prodRateGrid[k][CO] + prodRateGrid[k][H2O] + prodRateGrid[k][H2] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );
	
 // write source term CO2
	// fprintf( fp, "ProdRateCO2 [kg/m^3s]\n" );
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", prodRateGrid[k][CO2] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term H2O
	// fprintf( fp, "ProdRateH2O [kg/m^3s]\n" );
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", prodRateGrid[k][H2O] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term H2
	// fprintf( fp, "ProdRateH2 [kg/m^3s]\n" );
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", prodRateGrid[k][H2] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term CO
	// fprintf( fp, "ProdRateCO [kg/m^3s]\n" );
	// fprintf( fp, "\t%-.6e", 0.0 );
	// for ( k = 0; k < fNGridPoints; ++k ) {
	// 	fprintf( fp, "\t%-.6e", prodRateGrid[k][CO] );
	// 	if ( (k+2) % 5 == 0 ) {
	// 		fprintf( fp, "\n" );
	// 	}
	// }
	// fprintf( fp, "\t%-.6e\n", 0.0 );

		Free1DArray( summh );
		// Free1DArray( summh_C2H2 );
		if ( fProperties->GetRadiation() ) {
			Free1DArray( gasrad );
			Free1DArray( sootrad );
			Free1DArray( Kappa_rad );
			Free1DArray( Emission_rad );

			if (GasRadiationName == "WSGGJohansson" || GasRadiationName == "WSGGBordbar"){

				Free1DArray( Kappa_rad1 );
				Free1DArray( Kappa_rad2 );
				Free1DArray( Kappa_rad3 );
				Free1DArray( Kappa_rad4 );
				Free1DArray( Emission_rad1 );
				Free1DArray( Emission_rad2 );
				Free1DArray( Emission_rad3 );
				Free1DArray( Emission_rad4 );
			}
		}
		Free1DArray( density );
	
	if ( fSoot ) {
		for ( i = 0; i < nSootMoments; ++i ) {
			sprintf( tmpName, "M%d", i );
			PrintFlameletVector( fNGridPoints+2, &theMom[kPrev][i], tmpName
					, fp, fSootMomentsWork->phys_rows );
		}
		
		WriteSootInfo( theTime, temp, Y, GetPressure( theTime ), theMom, fp );		
	}

	Double	*totEnt = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k )
	{
		totEnt[k+1] = GetTotEnt( k, x, theTime );
	}
	PrintFlameletVector( fNGridPoints+2, totEnt, "TotalEnthalpy [J/kg]", fp );

	for ( k = -1; k <= fNGridPoints; ++k )
	{
		totEnt[k+1] = GetTotEnt( temp[k], Y[k] );
	}
	PrintFlameletVector( fNGridPoints+2, totEnt, "TotalEnthalpy2 [J/kg]", fp );

	Free1DArray( totEnt );

	PrintFlameletVector( fNGridPoints+2, &fMonFct[kPrev], "MonitorFunc", fp );

	
	// Production rate of all species         Rui 01/12/2019

	// int	index = -1;
	// for ( i = 0; i < nOfSpecies; ++i ) {
	// 	index = GetSpecies()->FindSpecies( names[i] );
	// 	if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// 	}
	// }

	// Production rate of 1-CH2 and 3-CH2         Rui 01/13/2019

	// int	index = -1;

	// index = GetSpecies()->FindSpecies( "1-CH2" );
	// if ( index > -1 ) {
	// PrintProdRate( theTime, index, fp );
	// }

	// index = GetSpecies()->FindSpecies( "3-CH2" );
	// if ( index > -1 ) {
	// PrintProdRate( theTime, index, fp );
	// }

	// Production Rate of C2H2
	// index = GetSpecies()->FindSpecies( "C2H2" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }

	// // Production Rate of H2
	// index = GetSpecies()->FindSpecies( "H2" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }

	// // Production Rate of H2O
	// index = GetSpecies()->FindSpecies( "H2O" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }

	// // Production Rate of C2H4
	// index = GetSpecies()->FindSpecies( "C2H4" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }

	// Production Rate of CO2
	// index = GetSpecies()->FindSpecies( "CO2" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }


	// // Production Rate of NO
	// index = GetSpecies()->FindSpecies( "NO" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }

	// Production Rate of NO
	// index = GetSpecies()->FindSpecies( "CH4" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }

	// Production Rate of NO
	// index = GetSpecies()->FindSpecies( "C2H2" );
	// if ( index > -1 ) {
	// 	PrintProdRate( theTime, index, fp );
	// }



/*********************Commmented out by Rui 12/18/2018*******************************************************************************************/

// #ifdef WRITEGLOBALPRODRATE
// 	PrintProdRateGlobalReac( theTime );	
// #endif
	
// 	fprintf( fp, "trailer\n" );
	
// #ifdef LEWISCHANGE
// 	Double					*LeOrig = fSpecies.GetLewisNumber()->vec;
// 	Double					*Le = new Double[nSpeciesIn];
// 	Double					LeFunc;
// 	Double					switchTime = LEWISSWITCHTIME;
// 	Double					switchPeriod = LEWISSWITCHPERIOD;
// 	const Double			Pi = 4.0 * atan( 1.0 );
// 	Double					omega = 2.0 * Pi / switchPeriod;
// 	Double					tStar = theTime - ( switchTime - 0.5 * switchPeriod );

// 	if ( tStar >= 0.0 ) {
// 		if ( tStar < switchPeriod ) {
// 			LeFunc = 0.5 * ( cos( omega * tStar ) + 1.0 );
// 			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
// 				Le[iLe] = LeFunc * ( LeOrig[iLe] - 1.0 ) + 1.0;
// 			}
// 		}
// 		else {
// 			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
// 				Le[iLe] = 1.0;
// 			}
// 		}
// 	}
// 	else {
// 		for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
// 			Le[iLe] = LeOrig[iLe];
// 		}
// 	}
// #else
// 	Double					*Le = fSpecies.GetLewisNumber()->vec;
// #endif
	// for ( i = 0; i < nSpeciesIn; ++i ) {
	// 	fprintf( fp, "%s\t%g\n", names[i], Le[i] );
	// }

/********************Added by Rui to check negative HRR*************************/
	// fprintf( fp, "Major Species for negative heat release:\n" );

	// for ( k = -1; k <= fNGridPoints; ++k ){
	// 	if ( MinHRRSpecies[k+1] != -1  ){
	// 		currSpeciesIndex = MinHRRSpecies[k+1];
	// 		fprintf( fp, "%s\n", names[currSpeciesIndex] );
	// 	}
	// }
/********************Added by Rui to check negative HRR*************************/

	if ( fpOpen ) {
		fclose( fp );
	}
#ifdef LEWISCHANGE
	delete Le;
#endif
	DisposeVector( ZBilgerVec );

/*********************Commmented out by Rui 12/18/2018*******************************************************************************************/

}

template<typename Species>
void TTransFlameSolver<Species>::PrintProdRate( Double t, int speciesIndex, FILE *fp )
{
	int			k;
	int			j;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	Double		*prodRate = fSpecies.GetProductionRate()->vec;
	Double		*reactionRate = fReaction.GetReactionRate()->vec;
	Double		pressure = GetPressure( t );
	Double		*prodOfZ = New1DArray( fNGridPoints+2 );
	Double		*prodPlusOfZ = New1DArray( fNGridPoints+2 );
	Double		*prodMinusOfZ = New1DArray( fNGridPoints+2 );
	Double		*prodMinusOverYOfZ = New1DArray( fNGridPoints+2 );
	char		buffer[128];
	sprintf( buffer, "ProdRate-%s", fSpecies.GetNames()[speciesIndex] );

	prodOfZ[0] = prodOfZ[fNGridPoints+1] = 0.0;
	prodPlusOfZ[0] = prodPlusOfZ[fNGridPoints+1] = 0.0;
	prodMinusOfZ[0] = prodMinusOfZ[fNGridPoints+1] = 0.0;
	prodMinusOverYOfZ[0] = prodMinusOverYOfZ[fNGridPoints+1] = 0.0;
	for ( k = 0; k < fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#else 
		T0DFlame<Species>::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#endif
		prodOfZ[k+1] = prodRate[speciesIndex];
		prodPlusOfZ[k+1] = fSpecies.GetPosNegProductionRate( speciesIndex, reactionRate, TRUE );
		prodMinusOfZ[k+1] = fSpecies.GetPosNegProductionRate( speciesIndex, reactionRate, FALSE );
		prodMinusOverYOfZ[k+1] = prodMinusOfZ[k+1] / MAX(1.0e-10, Y[k][speciesIndex] );
	}

	PrintFlameletVector( fNGridPoints+2, prodOfZ, buffer, fp );
	sprintf( buffer, "ProdRatePos-%s", fSpecies.GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodPlusOfZ, buffer, fp );
	sprintf( buffer, "ProdRateNeg-%s", fSpecies.GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodMinusOfZ, buffer, fp );
	sprintf( buffer, "ProdRateNegOverYNO-%s", fSpecies.GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodMinusOverYOfZ, buffer, fp );
	
	Free1DArray( prodMinusOverYOfZ );
	Free1DArray( prodPlusOfZ );
	Free1DArray( prodMinusOfZ );
	Free1DArray( prodOfZ );

// 	Double		*jReactionRate = New1DArray( fNGridPoints+2 );

// 	jReactionRate[0] = jReactionRate[fNGridPoints+1] = 0.0;

// 	int *usedReactions;
//   	int nOfUsedReactions;

//   	IntVectorPtr NofReactionsOfSpecies;   
//   	IntVectorPtr *jofSpeciesReactions;    
//   	NofReactionsOfSpecies = fSpecies.GetNOfUsedReactions();
//   	jofSpeciesReactions = fSpecies.GetUsedReactions();
//   	nOfUsedReactions = NofReactionsOfSpecies->vec[speciesIndex];
//   	usedReactions = jofSpeciesReactions[speciesIndex]->vec;

// 	for (j = 0; j < nOfUsedReactions; ++j) {

// 		sprintf( buffer, "jReactionRate-%s-%d", fSpecies.GetNames()[speciesIndex], usedReactions[j]);

// 		for ( k = 0; k < fNGridPoints; ++k ) {

// #ifdef DELTAINEW
// 		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
// 										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
// #else 
// 		T0DFlame<Species>::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
// 										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
// #endif	      
	
// 			jReactionRate[k+1] = reactionRate[usedReactions[j]];

// 	  	}

// 		PrintFlameletVector( fNGridPoints+2, jReactionRate, buffer, fp );

// 	}

}

template<typename Species>
Double TTransFlameSolver<Species>::ComputeMeanFromPDF( Double t, int speciesIndex, Double **PDF, Double *timePDF )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	Double		*Z = fSolGrid->vec;
	Double		*prodRate = fSpecies.GetProductionRate()->vec;
	Double		pressure = GetPressure( t );
	char		**names = fSpecies.GetNames();
	Double		sum, sumPDF;
	int 		indPDF = 0;

	Double	zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );

	if ( PDF ) {
		while ( timePDF[indPDF] < t ) ++indPDF;
	}
	else {
		return 0.0;
	}

	sumPDF = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k ) {
		sumPDF += ( Z[k+1] - Z[k-1] )*zR
		 	* Interpol( t, PDF[indPDF-1][k], timePDF[indPDF-1]
			, PDF[indPDF][k], timePDF[indPDF] );
	}
	sumPDF *= 0.5;

	sum = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k ) {
		sum += Y[k][speciesIndex] * ( Z[k+1] - Z[k-1] )*zR
		 	* Interpol( t, PDF[indPDF-1][k], timePDF[indPDF-1]
			, PDF[indPDF][k], timePDF[indPDF] );
	}
			
	return 0.5 * sum / sumPDF;
}

template<typename Species>
Double TTransFlameSolver<Species>::CheckPDF( Double t, Double **PDF, Double *timePDF )
{
	int			k;
	Double		*Z = fSolGrid->vec;
	Double		*prodRate = fSpecies.GetProductionRate()->vec;
	Double		pressure = GetPressure( t );
	char		**names = fSpecies.GetNames();
	Double		sum;
	int 		indPDF = 0;

	Double	zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );

	if ( PDF ) {
		while ( timePDF[indPDF] < t ) ++indPDF;
	}
	else {
		return 0.0;
	}

	sum = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k ) {
		sum += ( Z[k+1] - Z[k-1] )*zR
		 	* Interpol( t, PDF[indPDF-1][k], timePDF[indPDF-1]
			, PDF[indPDF][k], timePDF[indPDF] );
	}
			
	return 0.5 * sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::ComputeEmissionIndex( Double t, int speciesIndex, Double **PDF, Double *timePDF )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	Double		*Z = fSolGrid->vec;
	Double		*prodRate = fSpecies.GetProductionRate()->vec;
	Double		pressure = GetPressure( t );
	char		**names = fSpecies.GetNames();
	Double		sum;
	int 		indPDF = 0;

	Double	zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );

	if ( PDF ) {
		while ( timePDF[indPDF] < t ) ++indPDF;
	}

	sum = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#else 
		fprintf( fOutFilePtr, "UpdateThermoProps_ComputeEmissionIndex_4401\n" );
		T0DFlame<Species>::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#endif
		sum += prodRate[speciesIndex] * ( Z[k+1] - Z[k-1] )*zR
		 	* ( ( PDF ) ? Interpol( t, PDF[indPDF-1][k], timePDF[indPDF-1]
			, PDF[indPDF][k], timePDF[indPDF] ) : 1.0 );
	}
			
	return 0.5 * sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::ComputeEmissionIndexSoot( Double t, int which, Double **PDF, Double *timePDF )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	Double		*Z = fSolGrid->vec;
	Double		pressure = GetPressure( t );
	char		**names = fSpecies.GetNames();
	Double		sum;
	int 		indPDF = 0;

	Double	zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	
	if ( PDF ) {
		while ( timePDF[indPDF] < t ) ++indPDF;
	}

	sum = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#else 
		T0DFlame<Species>::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#endif // DELTAINEW
		sum += fSoot->NucleationNew( which, T[k], fSoot->GetPAHMoments()->vec ) 
				* ( Z[k+1] - Z[k-1] )*zR
				* ( ( PDF ) ? Interpol( t, PDF[indPDF-1][k], timePDF[indPDF-1], PDF[indPDF][k], timePDF[indPDF] ) : 1.0 );
	}
	
	return 0.5 * sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::TurbMeanZBarlow( Double t, Double ZMean, Double ZVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;
	
	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}


	Double		*pdf = New1DArray( nj );
	Double		sum;
	Double		**Y = fSolMassFracs->mat;

// Mixture fraction following Barlows definition
	static const Double	molarMassC = 12.01, 
						molarMassH = 1.008;
	static Double		elemMassLastCInit = 10000.0;
	static Double		elemMassLastHInit = 10000.0;
	if ( fZREnd >= 0.9999 ) {
		elemMassLastCInit = GetElementMassFraction( Y[fNGridPoints], "C", molarMassC );
		elemMassLastHInit = GetElementMassFraction( Y[fNGridPoints], "H", molarMassH );
	}
	Double	*ZBarlow = New1DArray( fNGridPoints+2 );
	ZBarlow = &ZBarlow[kNext];
	Double	elemMassFirstH = GetElementMassFraction( Y[kPrev], "H", molarMassH );
	Double	elemMassFirstC = GetElementMassFraction( Y[kPrev], "C", molarMassC );
	Double	denomBarlow = 0.5 * ( elemMassLastHInit - elemMassFirstH ) / molarMassH
				+ 2.0 * ( elemMassLastCInit - elemMassFirstC ) / molarMassC;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZBarlow[k] = ( 0.5 * ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirstH ) / molarMassH
						+ 2.0 * ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirstC ) / molarMassC )
						/ denomBarlow;
	}


	ZVar = MAX( 1e-10, ZVar );
	
	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newZBarlow = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	*ZBarlowBase0 = &ZBarlow[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newZBarlow[i*5+j] = ZBarlowBase0[i] + j * 0.2 * ( ZBarlowBase0[i+1] - ZBarlowBase0[i] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newZBarlow[i*5] = ZBarlowBase0[fNGridPoints + 1];


#ifdef CBETAPDF
	BetaPDF(newPoints, newGrid, ZMean, ZVar, newPdf, NULL);
#else
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newZBarlow[k];
	}
		
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free1DArray( newZBarlow );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	ZBarlow = &ZBarlow[kPrev];
	Free1DArray( ZBarlow );
	return sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::TurbMeanTemp( Double t, Double ZMean, Double ZVar, Double *tempVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;
	
	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
        Z[0]=0.0;
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	Double		*T = fSolTemp->vec;
	Double		*pdf = New1DArray( nj );
	Double		sum;


	ZVar = MAX( 1e-10, ZVar );
	
	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newTemp = New1DArray( newPoints ); 
	Double	*newPdf = New1DArray( newPoints );
	Double	*TBase0 = &T[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newTemp[i*5+j] = TBase0[i] + j * 0.2 * ( TBase0[i+1] - TBase0[i] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newTemp[i*5] = TBase0[fNGridPoints + 1];

#ifdef CBETAPDF
	BetaPDF(newPoints, newGrid, ZMean, ZVar, newPdf, NULL);
#else
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newTemp[k];
	}
	
	if ( tempVar ) {
		Double	sumVar = 0.0;
		for ( k = 0; k < newPoints; ++k ) {
			sumVar += newPdf[k] * ( newTemp[k] - sum ) * ( newTemp[k] - sum );
		}	
		*tempVar = sumVar;
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free1DArray( newTemp );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	return sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::TurbMeanTotEnergy( Double t, Double ZMean, Double ZVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	Double		*T = fSolTemp->vec;
	Double		**Y = fSolMassFracs->mat;
	Double		*pdf = New1DArray( nj );
	Double		sum;
	int 		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double		*molarMass = fSpecies.GetMolarMass()->vec;
	Double		dens;
	
	Double		*totEnt = New1DArray( fNGridPoints+2 );
	totEnt = &totEnt[kNext];

	for ( k = -1; k <= fNGridPoints; ++k )
	{
		fProperties->ComputeMixtureMolarMass( fProperties->GetMixMolarMassRef()
				, Y[k], molarMass, nSpeciesIn );
		dens = GetPressure( t ) * fProperties->GetMixMolarMass() / ( RGAS * T[k] );
		totEnt[k] = dens * ( GetTotEnt( k, fSolGrid->vec, t ) - GetTotEnt( T[k], Y[k] ) );
	}

	ZVar = MAX( 1e-10, ZVar );

	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newTotEnt = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	*HBase0 = &totEnt[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newTotEnt[i*5+j] = HBase0[i] + j * 0.2 * ( HBase0[i+1] - HBase0[i] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newTotEnt[i*5] = HBase0[fNGridPoints + 1];

#ifdef CBETAPDF
	BetaPDF(newPoints, newGrid, ZMean, ZVar, newPdf, NULL);
#else
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newTotEnt[k];
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	
	Free1DArray( newPdf );
	Free1DArray( newTotEnt );
	Free1DArray( newGrid );
	Free1DArray( &totEnt[kPrev] );
	Free1DArray( pdf );

	return sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::TurbMeanX( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = zeta[k] * zR;
		}
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	int 		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double		*molarMass = fSpecies.GetMolarMass()->vec;
	Double		**Y = fSolMassFracs->mat;
	Double		*pdf = New1DArray( nj );
	Double		sum;

	ZVar = MAX( 1e-10, ZVar );

	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	**newY = New2DArray( newPoints, nSpeciesIn );
	Double	*newPdf = New1DArray( newPoints );
	Double	**YBase0 = &Y[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			for ( k = 0; k < nSpeciesIn; ++k ) {
				newY[i*5+j][k] = YBase0[i][k] + j * 0.2 * ( YBase0[i+1][k] - YBase0[i][k] );
			}
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	for ( k = 0; k < nSpeciesIn; ++k ) {
		newY[i*5][k] = YBase0[fNGridPoints + 1][k];
	}

#ifdef CBETAPDF
	BetaPDF(newPoints, newGrid, ZMean, ZVar, newPdf, NULL);
#else
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		fProperties->ComputeMixtureMolarMass( fProperties->GetMixMolarMassRef()
				, newY[k], molarMass, nSpeciesIn );
		sum += newPdf[k] * newY[k][speciesIndex] * fProperties->GetMixMolarMass() / molarMass[speciesIndex];
	}

	if ( specVar ) {
		Double X;
		Double sumVar = 0.0;
		for ( k = 0; k < newPoints; ++k ) {
			X = newY[k][speciesIndex] * fProperties->GetMixMolarMass() / molarMass[speciesIndex];
			sumVar += newPdf[k] * ( X - sum ) * ( X - sum );
		}	
		*specVar = sumVar;
	}

	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free2DArray( newY );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	return sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::TurbMeanY( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
        Z[0]=0.0;
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	int 		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double		*molarMass = fSpecies.GetMolarMass()->vec;
	Double		**Y = fSolMassFracs->mat;
	Double		*pdf = New1DArray( nj );
	Double		sum;

	ZVar = MAX( 1e-10, ZVar );

	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newY = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	**YBase0 = &Y[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newY[i*5+j] = YBase0[i][speciesIndex] + j * 0.2 * ( YBase0[i+1][speciesIndex] - YBase0[i][speciesIndex] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newY[i*5] = YBase0[fNGridPoints + 1][speciesIndex];

#ifdef CBETAPDF
	BetaPDF(newPoints, newGrid, ZMean, ZVar, newPdf, NULL);
#else
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newY[k];
	}
	
	if ( specVar ) {
		Double	sumVar = 0.0;
		for ( k = 0; k < newPoints; ++k ) {
			sumVar += newPdf[k] * ( newY[k] - sum ) * ( newY[k] - sum );
		}	
		*specVar = sumVar;
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free1DArray( newY );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	return sum;
}

template<typename Species>
Double TTransFlameSolver<Species>::TurbMeanSoot( Double t, int sootIndex, Double ZMean, Double ZVar, Double *sootVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = zeta[k] * zR;
		}
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	int 		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double		*molarMass = fSpecies.GetMolarMass()->vec;
	Double		*pdf = New1DArray( nj );
	Double		sum;
	Double		**moments;

	if ( fSoot ) {
		moments = fSolSootMoments->mat;
	}
	else {
		return 0.0;
	}

	ZVar = MAX( 1e-10, ZVar );

#ifdef CBETAPDF
	BetaPDF(nj, Z, ZMean, ZVar, pdf, NULL);
#else
	F77BETAPDF( pdf, Z, &ZMean, &ZVar, &nj, &nj );
#endif

	sum = 0.0;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		sum += pdf[k+1] * moments[k][sootIndex];
	}

	if ( sootVar ) {
		Double	sumVar = 0.0;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			sumVar += pdf[k+1] * ( moments[k][sootIndex] - sum ) * ( moments[k][sootIndex] - sum );
		}	
		*sootVar = sumVar;
	}

	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( pdf );

	return sum;
}

template<typename Species>
FILE *TTransFlameSolver<Species>::GetOutputFile( Double theTime, const char *head, const char *tail, FileType type )
{
	int 		nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
	char		*name = new char[64];
	FILE 		*fp;
	
	sprintf( name, "%s%s%.8s_p%.2dt%9.3ems%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, fSpecies.GetNames()[GetFuelIndex()]
					, ( int ) floor( GetPressure( theTime ) * 1.0e-5 + 0.5 )	// in [bar]
					, theTime * 1000.0
					, ( tail ) ? tail : "" );
	
	fp = GetOutfile( name, type );
	delete[] name;
	
	return fp;
}

template<typename Species>
FILE *TTransFlameSolver<Species>::GetOutputFile( const char * /*head*/, const char * /*tail*/, const FileType /*type*/ )
{ 
	fprintf( fOutFilePtr, "wrong instance of function 'TTransFlameSolver<Species>::GetOutputFile'\n" );
	exit( 2 );

	return NULL;
}


template<typename Species>
void TTransFlameSolver<Species>::DoExit( void )
{
	FILE	*fp = GetOutfile( "interrupt", FileType::kText );
	WriteFlameletFile( fp, NULL, NULL );
	fprintf( fOutFilePtr, "\nprogram stopped by user\noutput written to %s\n"
				, GetOutfileName( "interrupt", FileType::kText ) );
	fclose( fp );
	exit( 2 );
}

template<typename Species>
Double TTransFlameSolver<Species>::GetCurrentTime( void )
{
	int	theTime = GetActualPoint( fTEnd );

	return ( theTime >= 0.0 ) ? fSolTime->vec[theTime] : fTEnd;
}

#ifdef NOCCLINK
extern "C" int _main();
#endif


template<typename Species>
Double TTransFlameSolver<Species>::GetDissRateReq( Double ZReq, Double dissRate, Double z )
{
	Double	zReqStar = ( ZReq - fZl ) / ( fZr - fZl );
	Double	zStar = ( z - fZl ) / ( fZr - fZl );
	Double	zReqNew = zReqStar - 0.5;
	Double	zReqNew2 = zReqNew * zReqNew;
	Double	zReqNew4 = zReqNew2 * zReqNew2;
	Double	zReqNew6 = zReqNew4 * zReqNew2;
	Double	zReqNew8 = zReqNew6 * zReqNew2;
	Double	zNew = zStar - 0.5;
	Double	zNew2 = zNew * zNew;
	Double	zNew4 = zNew2 * zNew2;
	Double	zNew6 = zNew4 * zNew2;
	Double	zNew8 = zNew6 * zNew2;
	Double	fitZReq = 1.00539 * ( 12.9041 - 
					82.123 * zReqNew2 +
                    115.29 * zReqNew4 -
                    201.898 * zReqNew6 +
                    912.136 * zReqNew8 );
	Double	fitZ = 1.00539 * ( 12.9041 - 
                    82.123  * zNew2 +
                    115.29  * zNew4 -
                    201.898 * zNew6 +
                    912.136 * zNew8 );
	
	
	return dissRate * fitZReq / fitZ;
}

#ifdef MOVEZRIGHT
template<typename Species>
Double TTransFlameSolver<Species>::GetDissFunc( Double z )
{
	Double	ZR = Interpol( GetCurrentTime(), fZRStart, fTStart, fZREnd, fTEnd );
	
	
	if ( z < 1.0e-10 || z > ZR ) {
		return 0.0;
	}

	return - ( z * z * log( z / ZR ) );
}

#else // !defined MOVEZRIGHT
template<typename Species>
Double TTransFlameSolver<Species>::GetDissFunc( Double z )
{
	Double	zStar = ( z - fZl ) / ( fZr - fZl );
	Double	zNew = zStar - 0.5;
	Double	zNew2 = zNew * zNew;
	Double	zNew4 = zNew2 * zNew2;
	Double	zNew6 = zNew4 * zNew2;
	Double	zNew8 = zNew6 * zNew2;

	return 1.00539 * ( 12.9041 - 
                    82.123  * zNew2 +
                    115.29  * zNew4 -
                    201.898 * zNew6 +
                    912.136 * zNew8 );
}
#endif // MOVEZRIGHT

template<typename Species>
void TTransFlameSolver<Species>::WriteSootInfo( Double theTime, Double *temp, Double **Y
				, Double pressure, Double **moments, FILE *fp )
{
	if ( fSoot ) {
		int		i, k;
		int		nSootMoments = fSoot->GetNSootMoments();
		const int	nPAHMolecules = 5;
		char	tmpName[127];
		
		// write soot volume fraction
		Double	fvFact = 24.0/1800.0; // fMolarMassSoot / fSootDensity;
		fprintf( fp, "fv\n" );
		for ( k = 0; k < fNGridPoints+2; ++k ) {
			fprintf( fp, "\t%-.6e", moments[k-1][1] * fvFact );
			if ( (k+1) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		if ( k % 5 ) {
			fprintf( fp, "\n" );
		}
				
		// write pah moments
		int			nPAHMoments = fSoot->GetNPAHMoments();
		int 		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
		VectorPtr	conv1Vec = NewVector( fNGridPoints+2 );
		Double		*conv1 = &conv1Vec->vec[kNext];
		VectorPtr	conv2Vec = NewVector( fNGridPoints+2 );
		Double		*conv2 = &conv2Vec->vec[kNext];
		VectorPtr	theCSootVec = NewVector( fNGridPoints+2 );
		Double		*theCSoot = &theCSootVec->vec[kNext];
		VectorPtr	theOxFactVec = NewVector( fNGridPoints+2 );
		Double		*theOxFact = &theOxFactVec->vec[kNext];
		MatrixPtr	thePAHMomMat = NewMatrix( nPAHMoments, fNGridPoints+2, kColumnPointers );
		Double		**thePAHMom = &thePAHMomMat->mat[kNext];
		TensorPtr	sourcesMat = NewTensor( nSootMoments, kNSootSources+2, fNGridPoints+2, kColumnPointers );
		Double		***sources = sourcesMat->tensor;
		MatrixPtr	pahConcsMat1 = NewMatrix( nPAHMolecules, fNGridPoints+2, kColumnPointers );
		Double		**pahConcs1 = &pahConcsMat1->mat[kNext];
		MatrixPtr	pahConcsMat2 = NewMatrix( nPAHMolecules, fNGridPoints+2, kColumnPointers );
		Double		**pahConcs2 = &pahConcsMat2->mat[kNext];
		Double		**si = NULL;
		Double		*sik = NULL;
		Double		*pahMoments = fSoot->GetPAHMoments()->vec;
		Double		**pahConcs = fSoot->GetPij()->mat;
		Double		*molarMass = fSpecies.GetMolarMass()->vec;
		Double		*Z = fSolGrid->vec;
		Double		density;

		for ( k = -1; k <= fNGridPoints; ++k )
		{
#ifdef DELTAINEW
			UpdateThermoProps( k, Y[k], temp[k], pressure, density, kDensFromPress, moments[k] ); 
#else
			T0DFlame<Species>::UpdateThermoProps( Y[k], temp[k], pressure, density, kDensFromPress, moments[k] );
#endif
			// save pah moments
			copy( nPAHMoments, pahMoments, 1, thePAHMom[k], 1 );
			// save pah concentrations
			copy( nPAHMolecules, &pahConcs[0][1], 1, pahConcs1[k], 1 );
			copy( nPAHMolecules, &pahConcs[1][1], 1, pahConcs2[k], 1 );
			// save oxidation factor
			theOxFact[k] = fSoot->GetSootOxCoeff( Y[k], density, molarMass );
			theCSoot[k] = fSoot->GetCCSootStar();
			// save souce terms
			for ( i = 0; i < nSootMoments; ++i ) {
				sik = sources[i][k+1];
				if ( fSoot->WithNucleation() ) {
					sik[kNuc] = fSoot->NucleationNew( i, temp[k], pahMoments ) / density;
				}
				if ( fSoot->WithCoagulation() ) {
					sik[kCoag] = fSoot->SourceCoagulationNew( i, temp[k], moments[k] ) / density;
				}
				if ( fSoot->WithCondensation() ) {
					sik[kCond] = fSoot->SourceCondensationNew( i, temp[k], fSoot->GetPAHMoments()->vec
									, moments[k], Y[k], density, molarMass ) / density;
				}
				if ( fSoot->WithSurfaceGrowth() ) {
					sik[kSG] = fSoot->SourceSurfGrowthNew( i, moments[k], Y[k], density, molarMass ) / density;
				}
				if ( fSoot->WithSurfaceOxidation() ) {
					sik[kOx] = fSoot->SourceSootOxidationNew( i, moments[k], Y[k], density, molarMass ) / density;
				}
				if ( k > -1 && k < fNGridPoints )
				{
					Double	lambdaOverCpCurr = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();
					Double	fracIndex = i - 2.0 / 3.0;
					Double	rhoPrev, rhoNext, mm;
					Double	lambdaOverCpPrev, lambdaOverCpNext;
					fProperties->ComputeMixtureMolarMass( mm, Y[k-1], molarMass, nSpeciesIn );
					rhoPrev = pressure * mm / ( RGAS * temp[k-1] );
					fSpecies.ComputeSpeciesProperties( temp[k-1] );
					ComputeDeltaI( k-1, Y[k-1], temp[k-1] );
					fProperties->CompMixtureProps( fSpecies.GetHeatCapacity()->vec
					, fSpecies.GetConductivity()->vec, fSpecies.GetViscosity()->vec, Y[k-1], temp[k-1]
					, pressure, fProperties->GetDensityRef()
					, kDensFromPress, nSpeciesIn, &fSpecies );
					lambdaOverCpPrev = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();

					fProperties->ComputeMixtureMolarMass( mm, Y[k+1], molarMass, nSpeciesIn );
					rhoNext = pressure * mm / ( RGAS * temp[k+1] );
					fSpecies.ComputeSpeciesProperties( temp[k+1] );
					ComputeDeltaI( k+1, Y[k+1], temp[k+1] );
					fProperties->CompMixtureProps( fSpecies.GetHeatCapacity()->vec
					, fSpecies.GetConductivity()->vec, fSpecies.GetViscosity()->vec, Y[k+1], temp[k+1]
					, pressure, fProperties->GetDensityRef()
					, kDensFromPress, nSpeciesIn, &fSpecies );
					lambdaOverCpNext = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();

					fProperties->ComputeMixtureMolarMass( mm, Y[k], molarMass, nSpeciesIn );
					rhoNext = pressure * mm / ( RGAS * temp[k+1] );
					fSpecies.ComputeSpeciesProperties( temp[k+1] );
					ComputeDeltaI( k+1, Y[k+1], temp[k+1] );
					fProperties->CompMixtureProps( fSpecies.GetHeatCapacity()->vec
					, fSpecies.GetConductivity()->vec, fSpecies.GetViscosity()->vec, Y[k+1], temp[k+1]
					, pressure, fProperties->GetDensityRef()
					, kDensFromPress, nSpeciesIn, &fSpecies );

#ifdef SOOTCONVECTION
					conv1[k] = 0.25 / density * ( 
										fFDWMinus[k] * rhoPrev * GetDissRate( theTime, Z[k-1] ) 
										+ fFDWCurr[k] * density * GetDissRate( theTime, Z[k] )
										+ fFDWPlus[k] * rhoNext * GetDissRate( theTime, Z[k+1] ) );
					conv2[k] = 0.25 / density * ( 
										density * GetDissRate( theTime, Z[k] ) / lambdaOverCpCurr
											* ( fFDWMinus[k] * lambdaOverCpPrev
												+ fFDWCurr[k] * lambdaOverCpCurr
												+ fFDWPlus[k] * lambdaOverCpNext ) );
					Double	convVelo = 0.25 / density * ( 
										fFDWMinus[k] * rhoPrev * GetDissRate( theTime, Z[k-1] ) 
										+ fFDWCurr[k] * density * GetDissRate( theTime, Z[k] )
										+ fFDWPlus[k] * rhoNext * GetDissRate( theTime, Z[k+1] )
										+ density * GetDissRate( theTime, Z[k] ) / lambdaOverCpCurr
											* ( fFDWMinus[k] * lambdaOverCpPrev
												+ fFDWCurr[k] * lambdaOverCpCurr
												+ fFDWPlus[k] * lambdaOverCpNext ) );
#	ifdef UPWINDSOOT
				if ( convVelo > 0.0 ) {
					sik[kNSootSources+1] = -convVelo * ( ( moments[k][i] / density - moments[k-1][i] / rhoPrev ) 
							/ ( Z[k] - Z[k-1] )
							- 1.0 / fSoot->GetLewis1() 
							* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, moments[k+1] ) / rhoNext 
								- fSoot->FracMom2( fracIndex, moments[k] ) / density ) 
								/ ( Z[k+1] - Z[k] ) );
				}
				else {
					sik[kNSootSources+1] = -convVelo * ( ( moments[k+1][i] / rhoNext - moments[k][i] / density ) 
							/ ( Z[k+1] - Z[k] )
							- 1.0 / fSoot->GetLewis1() 
							* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, moments[k] ) / density
							- fSoot->FracMom2( fracIndex, moments[k-1] ) / rhoPrev ) 
								/ ( Z[k] - Z[k-1] ) );
				}
#	else
				sik[kNSootSources+1] = -convVelo * ( 
					1.0 / fSoot->GetLewis1() 
					* FRACMOMFACT * ( fFDWMinus[k] * fSoot->FracMom2( fracIndex, moments[k-1] ) / rhoPrev
					+ fFDWCurr[k] * fSoot->FracMom2( fracIndex, moments[k] ) / density
					+ fFDWPlus[k] * fSoot->FracMom2( fracIndex, moments[k+1] ) / rhoNext )
										
					- ( fFDWMinus[k] * moments[k-1][i] / rhoPrev
					+ fFDWCurr[k] * moments[k][i] / density
					+ fFDWPlus[k] * moments[k+1][i] / rhoNext )
										);
#	endif
#endif
				}
				else {
					sik[kNSootSources+1] = 0.0;
				}

				if ( k > -1 && k < fNGridPoints ) {
#	ifdef PRONE
					Double	PrCurr = 1.0;
#	else
					Double	PrCurr = fProperties->GetMixViscosity() / ( fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity() );
#	endif
					sik[kTherPhor] = 0.275 * moments[k][i] / density / ( temp[k] ) * PrCurr * GetDissRate( theTime, Z[k] )
									 * ( fWMinus[k] * temp[k-1]
										+fWCurr[k] * temp[k]
										+fWPlus[k] * temp[k+1] );
				}
				if ( k > -1 && k < fNGridPoints ) {
					Double	rhoPrev, rhoNext, mm;
					fProperties->ComputeMixtureMolarMass( mm, Y[k-1], molarMass, nSpeciesIn );
					rhoPrev = pressure * mm / ( RGAS * temp[k-1] );
					fProperties->ComputeMixtureMolarMass( mm, Y[k+1], molarMass, nSpeciesIn );
					rhoNext = pressure * mm / ( RGAS * temp[k+1] );
#ifdef SIZEDEPDIFFUSION
					sik[kNSootSources] = GetDissRate( theTime, Z[k] ) / ( 2.0 * fSoot->GetLewis1() ) 
						* FRACMOMFACT * ( fWMinus[k] * fSoot->FracMom2( i - 2.0 / 3.0, moments[k-1] ) / rhoPrev
						+ fWCurr[k] * fSoot->FracMom2( i - 2.0 / 3.0, moments[k] ) / density
						+ fWPlus[k] * fSoot->FracMom2( i - 2.0 / 3.0, moments[k+1] ) / rhoNext );
#else
					sik[kNSootSources] = GetDissRate( theTime, Z[k] ) / ( 2.0 * fSoot->GetLewis1() ) 
						* ( fWMinus[k] * moments[k-1][i] / rhoPrev
						+ fWCurr[k] * moments[k][i] / density
						+ fWPlus[k] * moments[k+1][i] / rhoNext );
#endif
				} else {
					sik[kNSootSources] = 0.0;
				}
			}
		}

		PrintFlameletVector( fNGridPoints+2, theOxFact, "OxidationFactor"
				, fp, 1 );
		PrintFlameletVector( fNGridPoints+2, theCSoot, "CSootStar"
				, fp, 1 );

		PrintFlameletVector( fNGridPoints+2, conv1, "ConvVelo1", fp, 1 );
		PrintFlameletVector( fNGridPoints+2, conv2, "ConvVelo2", fp, 1 );
		
		for ( i = 0; i < nPAHMoments; ++i ) {
			sprintf( tmpName, "MPAH%d", i );
			PrintFlameletVector( fNGridPoints+2, &thePAHMom[kPrev][i], tmpName
					, fp, thePAHMomMat->phys_rows );
		}

		for ( i = 0; i < nSootMoments; ++i ) {
			si = sources[i];
			sprintf( tmpName, "So_Nuc_%d", i );
			PrintFlameletVector( fNGridPoints+2, &si[kCurr][kNuc], tmpName, fp, sourcesMat->phys_rows );
			if ( fSoot->WithCoagulation() ) {
				sprintf( tmpName, "So_Coag_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kCoag], tmpName, fp, sourcesMat->phys_rows );
			}
			if ( fSoot->WithCondensation() ) {
				sprintf( tmpName, "So_Cond_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kCond], tmpName, fp, sourcesMat->phys_rows );
			}
			if ( fSoot->WithSurfaceGrowth() ) {
				sprintf( tmpName, "So_SG_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kSG], tmpName, fp, sourcesMat->phys_rows );
			}
			if ( fSoot->WithSurfaceOxidation() ) {
				sprintf( tmpName, "So_Ox_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kOx], tmpName, fp, sourcesMat->phys_rows );
			}
//			if ( fSoot->WithThermoPhoresis() ) {
				sprintf( tmpName, "So_TherPhor_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kTherPhor], tmpName, fp, sourcesMat->phys_rows );
//			}
#ifdef SOOTCONVECTION
			sprintf( tmpName, "So_Conv_%d", i );
			PrintFlameletVector( fNGridPoints+2, &si[kCurr][kNSootSources+1], tmpName, fp, sourcesMat->phys_rows );
#endif // SOOTCONVECTION
			sprintf( tmpName, "So_Diff_%d", i );
			PrintFlameletVector( fNGridPoints+2, &si[kCurr][kNSootSources], tmpName, fp, sourcesMat->phys_rows );
		}

		for ( i = 0; i <= nPAHMolecules; ++i ) {
			sprintf( tmpName, "P0%d", i+1 );
			PrintFlameletVector( fNGridPoints+2, &pahConcs1[kPrev][i], tmpName
					, fp, pahConcsMat1->phys_rows );
		}

		for ( i = 0; i <= nPAHMolecules; ++i ) {
			sprintf( tmpName, "P1%d", i+1 );
			PrintFlameletVector( fNGridPoints+2, &pahConcs2[kPrev][i], tmpName
					, fp, pahConcsMat2->phys_rows );
		}

		DisposeMatrix( pahConcsMat2 );
		DisposeMatrix( pahConcsMat1 );
		DisposeTensor( sourcesMat );
		DisposeMatrix( thePAHMomMat );
		DisposeVector( theOxFactVec );
		DisposeVector( theCSootVec );
		DisposeVector( conv1Vec );
		DisposeVector( conv2Vec );
	}
}

#ifdef LOGNORMALCHI
template<typename Species>
void TTransFlameSolver<Species>::SetRandomNumber( void )
{
	static int	idum = 0;

	fRandomNumber = exp( gasdev( &idum ) );
}
#endif // LOGNORMALCHI

template<typename Species>
int TTransFlameSolver<Species>::GetSootSources( int whichMoment, ConstStringArray names, Double **sources, Double *grid, int gridPointsA, int vars )
{
	if ( !fSoot ) {
		return 1;
	}
	if ( whichMoment >= fSoot->GetNSootMoments() || whichMoment < 0 ) {
		fprintf( stderr, "###error in function 'GetSootSources': soot moment no. %d not allowed\n", whichMoment );
	}
	
	int			i, k, ind;
	int			nOfSpecies = fSpecies.GetNOfSpecies();
	int			nSootMoments = fSoot->GetNSootMoments();
	MMDataBag	bag( kNSootSources );
	Double		**Y = fMassFracsWork->mat;
	Double		*temp = fTempWork->vec;
	Double		**moments = fSootMomentsWork->mat;
	Double		theTime = GetCurrentTime();
	MatrixPtr	sourceWork = NewMatrix( kNSootSources, fNGridPoints+2, kColumnPointers );
	Double		**theSources = sourceWork->mat;
	VectorPtr	fOutSol = NewVector( gridPointsA );
	Flag		*outSet = new Flag[vars];
	if ( !outSet ) fprintf( stderr, "new failed\n" );
	char		theNames[kNSootSources][20];
	for ( i = 0; i < vars; ++i ) {
		outSet[i] = FALSE;
	}
	// init names
	strcpy( theNames[kNuc], "Nucleation" );
	strcpy( theNames[kCoag], "Coagulation" );
	strcpy( theNames[kCond], "Condensation" );
	strcpy( theNames[kSG], "SurfaceGrowth" );
	strcpy( theNames[kOx], "SurfaceOxidation" );
	strcpy( theNames[kTherPhor], "Thermophoresis" );
	
// copy solution to bag
	bag.Initialize();
	bag.SetOldInpedVar( &fSolGrid->vec[kPrev], fSolGrid->len, 1, "xIn" );

	SetOutSolution();
	
	SetSootSources( whichMoment, sourceWork
				, theTime, temp, Y
				, GetPressure( theTime ), moments );	
	
	for ( i = 0; i < kNSootSources; ++i ) {
		bag.Insert( &theSources[0][i], fNGridPoints+2, kNSootSources, theNames[i] );
	}
	
	bag.SetNewInpedVar( grid, gridPointsA, 1, "xNew" );

#ifdef DEBUGBAG
	std::cout << bag;
#endif // DEBUGBAG

	Double		**newSource = &sources[kNext];
	Double		*outWork = &fOutSol->vec[kNext];
	for ( i = 0; i < bag.NumElems(); ++i ) {
		ind = GetVariableIndex( bag[i].Name(), names, vars );
		if ( ind >= 0 ) {
			bag[i].Map( fOutSol );
			for ( k = -1; k < gridPointsA-1; ++k ) {
				newSource[k][ind] = outWork[k];
			}
			outSet[ind] = TRUE;
		}
		else {
		}
	}
	
	for ( i = 0; i < vars; ++i ) {
		if ( outSet[i] == FALSE ) {
			fprintf( fOutFilePtr, "#warning from function 'GetSolution': no match for source term %s\n", names[i] );
		}
	}
	
#ifdef DEBUGSOURCES
	FILE	*fp = GetOutfile( "SootSources", FileType::kData );
	PrintSolution( fp, gridPointsA-2, vars, &grid[kNext], newSource, names );
	fclose( fp );
#endif // DEBUGSOURCES

	// clean up
	delete[] outSet;
	DisposeVector( fOutSol );
	DisposeMatrix( sourceWork );
	
	return 0;
}

template<typename Species>
void TTransFlameSolver<Species>::SetSootSources( int whichMoment, MatrixPtr sourcesMat
				, Double theTime, Double *temp, Double **Y
				, Double pressure, Double **moments )
{
	int			k;
	int			nPAHMoments = fSoot->GetNPAHMoments();
	int 		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	MatrixPtr	thePAHMomMat = NewMatrix( nPAHMoments, fNGridPoints+2, kColumnPointers );
	Double		**sources = &sourcesMat->mat[kNext];
	Double		*sik = NULL;
	Double		*pahMoments = fSoot->GetPAHMoments()->vec;
	Double		*molarMass = fSpecies.GetMolarMass()->vec;
	Double		*Z = fSolGrid->vec;
	Double		density;

	for ( k = -1; k <= fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], temp[k], pressure
						, density, kDensFromPress, moments[k] );
#else 
		T0DFlame<Species>::UpdateThermoProps( Y[k], temp[k], pressure
						, density, kDensFromPress, moments[k] );
#endif // DELTAINEW

		// save souce terms
		sik = sources[k];
		sik[kNuc] = fSoot->NucleationNew( whichMoment, temp[k], pahMoments );
		
		if ( fSoot->WithCoagulation() ) {
			sik[kCoag] = fSoot->SourceCoagulationNew( whichMoment, temp[k], moments[k] );
		}
		else {
			sik[kCoag] = 0.0;
		}
		
		if ( fSoot->WithCondensation() ) {
			sik[kCond] = fSoot->SourceCondensationNew( whichMoment, temp[k], fSoot->GetPAHMoments()->vec
							, moments[k], Y[k], density, molarMass ); 
		}
		else {
			sik[kCond] = 0.0;
		}
		
		if ( fSoot->WithSurfaceGrowth() ) {
			sik[kSG] = fSoot->SourceSurfGrowthNew( whichMoment, moments[k], Y[k], density, molarMass );
		}
		else {
			sik[kSG] = 0.0;
		}
		
		if ( fSoot->WithSurfaceOxidation() ) {
			sik[kOx] = fSoot->SourceSootOxidationNew( whichMoment, moments[k], Y[k], density, molarMass );
		}
		else {
			sik[kOx] = 0.0;
		}
		
		if ( fSoot->WithThermoPhoresis() ) {
			if ( k > -1 && k < fNGridPoints )
			{
				Double	fracIndex = whichMoment - 2.0 / 3.0;
				Double	rhoPrev, rhoNext, mm;
				fProperties->ComputeMixtureMolarMass( mm, Y[k-1], molarMass, nSpeciesIn );
				rhoPrev = pressure * mm / ( RGAS * temp[k-1] );
				fProperties->ComputeMixtureMolarMass( mm, Y[k+1], molarMass, nSpeciesIn );
				rhoNext = pressure * mm / ( RGAS * temp[k+1] );
				Double	convVelo = 0.25 * ( fFDWMinus[k] * rhoPrev * GetDissRate( theTime, Z[k-1] ) 
									+ fFDWCurr[k] * density * GetDissRate( theTime, Z[k] )
									+ fFDWPlus[k] * rhoNext * GetDissRate( theTime, Z[k+1] ) );
				sik[kTherPhor] = convVelo * ( fFDWMinus[k] * moments[k-1][whichMoment] / rhoPrev
									+ fFDWCurr[k] * moments[k][whichMoment] / density
									+ fFDWPlus[k] * moments[k+1][whichMoment] / rhoNext
									- 1.0 / fSoot->GetLewis1() 
									* FRACMOMFACT * ( fFDWMinus[k] * fSoot->FracMom2( fracIndex, moments[k-1] ) / rhoPrev
									+ fFDWCurr[k] * fSoot->FracMom2( fracIndex, moments[k] ) / density
									+ fFDWPlus[k] * fSoot->FracMom2( fracIndex,moments[k+1] ) / rhoNext ) );
			}
			else {
				sik[kTherPhor] = 0.0;
			}
		}
	}
}

template<typename Species>
void TTransFlameSolver<Species>::ReadStartProfiles( TInputDataPtr inp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure
				, Double *chi, Double *theTime, Double *currTimeStep
				, int tempOff, int speciesOff, int sootOff )
{
	StartProfilePtr	sp = NULL;
	char 	*insp = inp->fStartProfileFile;
	FILE	*fpS = NULL;
	char	*fileName = GetFullPath( insp, kFileName );

	sp = new StartProfile;	
	if ( !sp ) FatalError( "new failed" );
	if ( !insp  || ( fpS = fopen( fileName, "r" ) ) == NULL ) {
		fprintf( stderr, "startprofiles file '%s' not found\n", fileName );
		exit( 2 );
	} 
	else {
		fprintf( stderr, "read initial solution from file '%s'\n", fileName );
		::ReadStartProfiles( sp, fpS );
		SetInitialValues( fInputData, sp, nGridPoints, nVars, x, y
						, pressure, chi, theTime, currTimeStep
						, tempOff, speciesOff, sootOff );
		CleanReadStartProfile();
		delete sp;
		fclose( fpS );
	}
	delete fileName;
}

template<typename Species>
void TTransFlameSolver<Species>::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure, Double *chi
				, Double *theTime, Double *currTimeStep
				, int tempOff, int speciesOff, int sootOff )
{
	x = &x[kNext];
	y = &y[kNext];

	int 				i, j, k;
	int					variable, speciesIndex;
	Flag				ZSet = FALSE;
	int					gridPointsIn = sp->gridPoints;	// including boundaries
	int					nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	Double				*yInFloat = sp->data;
	char				*string = sp->labels;
	struct _parameter	*param;
	
	if ( nGridPoints < gridPointsIn ) {
		fprintf( stderr, "nGridPointsIn = %d > nGridPoints = %d\n", gridPointsIn, nGridPoints );
		exit(2);
	}

// find independent coordinate
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "z", 1 ) == 0 && ZSet == FALSE ) {
			cerr << "choose inputGrid" << NEWL;
			for ( j = -1; j <= gridPointsIn-2; ++j ) {
				x[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
			}
			ZSet = TRUE;
		}
		string += strlen(string) + 1;
	}

// error checking
	if ( !ZSet ) {
		cerr << "error: can't find coordinate 'Z'" << NEWL;
		exit(2);
	}

	SetInitialValues( nGridPoints-2, nVars, x, y ); // default values
	
// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = tempOff;
		}
		else if ( strncmp( string, "massfraction-", 13 ) == 0 ){
			string += 13;
			UpperString( string );
			if ( ( speciesIndex = inp->FindSpecies( string ) ) >= 0 ) {
				if ( speciesIndex < nSpeciesInSystem ) {
					variable = speciesOff + speciesIndex;
				}
				else {
					string += strlen(string) + 1;
					continue;
				}
			}
			else {
				cerr << "warning: no match for species " << string << NEWL;
				string += strlen(string) + 1;
				continue;
			}
		}
		else if ( sootOff >= 0 && GetSoot() && strncmp( string, "m", 1 ) == 0 ) {
			string += 1;
			int	num = atoi( string );
			if ( isdigit( string[0] ) && num < GetSoot()->GetNSootMoments() ) {
				variable = num + sootOff;
			}
			else {
				string += strlen(string) + 1;
				continue;
			}
		}
		else {
			string += strlen(string) + 1;
			continue;
		}

		string += strlen(string) + 1;
		for ( k = -1; k <= gridPointsIn-2; ++k ) {
			y[k][variable] = yInFloat[i*gridPointsIn + k + 1];	// copy workspace to vector of solution
		}
	}
		
	param = GetParameter( "pressure" );
	if ( param ) {
		*pressure = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
			*pressure *= 1.0e5;
		}
	}

	param = GetParameter( "chi" );
	if ( param ) {
		*chi = (Double)param->what.quantity.value;
	}
	else {
		param = GetParameter( "chi_ref" );
		if ( param ) {
			*chi = (Double)param->what.quantity.value;
		}
	}

	param = GetParameter( "time" );
	if ( param ) {
		*theTime = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			*theTime *= 1.0e-3;
		}
	}

	param = GetParameter( "currtimestep" );
	if ( param ) {
		*currTimeStep = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			*currTimeStep *= 1.0e-3;
		}
	}
}

template<typename Species>
void TTransFlameSolver<Species>::SetInitialValues( int nGridPoints, int nVars, Double *x, Double **y )
{
	Double	*yLeft = y[kPrev];
	Double	*yRight = y[nGridPoints];
	Double	left = x[kPrev];
	Double	right = x[nGridPoints];

	for ( int k = 0; k < nGridPoints; ++k ) {
		for ( int j = 0; j < nVars; ++j ) {
			y[k][j] = Interpol( x[k], yLeft[j], left, yRight[j], right );
		}
	}
}

template<typename Species>
Double TTransFlameSolver<Species>::GetTotEnt( int k, Double *Z, Double t )
{
	Double	totentstart = Interpol( Z[k] , fTotEntStart->vec[kPrev], Z[kPrev], fTotEntStart->vec[fNGridPoints], Z[fNGridPoints] );
	Double	totentend = Interpol( Z[k] , fTotEntEnd->vec[kPrev], Z[kPrev] , fTotEntEnd->vec[fNGridPoints], Z[fNGridPoints] );

	return Interpol( t, totentstart, fTStart, totentend, fTEnd );
}

#ifdef INGOOPT
template<typename Species>
void TTransFlameSolver<Species>::UpdateThermoProps( int k, Double *Y, Double temp, Double &pressure
								, Double &density, EqOfState what, Double *sootMoments )
{
	int		i;
	int		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Flag	specChanged = FALSE;
		
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( Y[i] != fYGSave[k][i] ) {
			specChanged = TRUE;
			copy( nSpeciesIn, Y, 1, fYGSave[k], 1 );
			break;
		}
	}

	if ( specChanged || temp != fTempGSave[k] ) {
		ComputeDeltaI( k, Y, temp );
		T0DFlame<Species>::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments );
		return;
	}

	copy( nSpeciesIn, fDeltaI[k], 1, fSpecies.GetDeltaI()->vec, 1 );
	T0DFlame<Species>::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments );
}

#else // !defined INGOOPT

template<typename Species>
void TTransFlameSolver<Species>::UpdateThermoProps( int k, Double *Y, Double temp, Double &pressure
								, Double &density, EqOfState what, Double *sootMoments )
{
	ComputeDeltaI( k, Y, temp );
	T0DFlame<Species>::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments ); 
}
#endif // INGOOPT

template<typename Species>
void TTransFlameSolver<Species>::ComputeDeltaI( int k, Double *Y, Double temp )
{
	int 	i;
	int		nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	Double	*deltaI = fSpecies.GetDeltaI()->vec;
	
	CompDeltaIG( k, temp );

	for ( i = 0; i < nSpeciesIn; ++i ) {
		deltaI[i] = CompOneDeltaI( i, k, Y );
	}
	copy( nSpeciesIn, deltaI, 1, fDeltaI[k], 1 );
}

template<typename Species>
Double TTransFlameSolver<Species>::CompOneDeltaI( int i, int k, Double *Y )
{
	int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	Double	*M = fSpecies.GetMolarMass()->vec;
	Double	Delta_i = 0;

	for ( int j = 0; j < nSpeciesInSystem; ++j ){
		Delta_i += fG_ij[k][i][j] / M[j] * Y[j];
	}
	
	return Delta_i * M[i];
}

template<typename Species>
void TTransFlameSolver<Species>::CompDeltaIG( int k, Double temp )
{
	int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	Double	*mu = fSpecies.GetViscosity()->vec;
	
	if ( temp == fTempGSave[k] ) {
		return;
	}
	
	fTempGSave[k] = temp;

	fSpecies.ComputeSpeciesProperties( temp );

	Double	**sqrtMjOverMi = fSpecies.GetSqrtMjOverMi();
	Double	**sqrtMiOverMjPlusOne = fSpecies.GetSqrtMiOverMjPlusOne();
	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
		for ( int j = 0; j < nSpeciesInSystem; ++j ) {
			fG_ij[k][i][j] = sqrt( sqrtMjOverMi[j][i] * mu[i] / mu[j] ) + 1.0; // this is after bug fix
			fG_ij[k][i][j] *= fG_ij[k][i][j] * sqrtMiOverMjPlusOne[i][j];
		}
	}
}

template<typename Species>
Double TTransFlameSolver<Species>::GetU( Double t, Double Z )
{
	int	i = 1;
	int	j = 1;
	int	zlen = fZIn->len;
	int	tlen = fTimeIn->len;
	Double	*zIn = fZIn->vec;
	Double	*tIn = fTimeIn->vec;
	Double	**UIn = fUstOverU->mat;
	Double	UAtZt1, UAtZt2;

	while ( j < zlen-1 && Z > zIn[j] ) ++j;
	if ( j == zlen-1 && Z > zIn[j]+1.0e-8 ) {
		fprintf( stderr, "Z = %g out of range[%g,%g]\n", Z, zIn[0]
								, zIn[zlen-1] );
	}
	
	while ( i < tlen-1 && t > tIn[i] ) ++i;
	if ( i == tlen-1 && t > tIn[i] ) {
		fprintf( stderr, "t = %g out of range[%g,%g]: linear extrapolation\n", Z, tIn[0]
								, tIn[tlen-1] );
	}
	
	if ( zIn[j] == zIn[j-1] ) {
		fprintf( stderr, "something's wrong\n" );
	}
	else {
		UAtZt1 = Interpol( Z, UIn[i-1][j-1], zIn[j-1], UIn[i-1][j], zIn[j] );
		UAtZt2 = Interpol( Z, UIn[i][j-1], zIn[j-1], UIn[i][j], zIn[j] );
	}
	
	if ( tIn[i] == tIn[i-1] ) {
		fprintf( stderr, "something else's wrong\n" );
		return 0.0;
	}
	else {
		return Interpol( t, UAtZt1, tIn[i-1], UAtZt2, tIn[i] );
	}
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

static Double ran1( int *idum )
{
	static int		iff = 0;
	static long		ix1, ix2, ix3;
	static Double	r[98];
	Double			temp;
	int				j;

	if (*idum < 0 || iff == 0) {
		iff=1;
		ix1=(IC1-(*idum)) % M1;
		ix1=(IA1*ix1+IC1) % M1;
		ix2=ix1 % M2;
		ix1=(IA1*ix1+IC1) % M1;
		ix3=ix1 % M3;
		for (j=1;j<=97;j++) {
			ix1=(IA1*ix1+IC1) % M1;
			ix2=(IA2*ix2+IC2) % M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum=1;
	}
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	ix3=(IA3*ix3+IC3) % M3;
	j=1 + ((97*ix3)/M3);
	if (j > 97 || j < 1) fprintf(stderr,"RAN1: This cannot happen.");
	temp=r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3


static Double ran1( int *idum );
static Double expdev( int *idum )
{
	return -log(ran1(idum));
}


static Double gasdev( int *idum )
{
	static int iset=0;
	static Double gset;
	Double fac,r,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

template<typename Species>
Double TTransFlameSolver<Species>::GetRosseRadiation( int k, Double *nTemp, Double **moments, Double rfPrev
											, Double rfCurr, Double rfNext )
{
// rf is rho * chi * cp / lambda
  Double rossecutoff = 5.0e-7;
  if ( moments[kCurr][1]*24.0/1800.0 > rossecutoff ) { 
					// coeffs negative
	Double	coeff = fSoot->GetSootRadRossCoeff( nTemp[kCurr], moments[kCurr], rossecutoff );
	Double	coeffPlus = fSoot->GetSootRadRossCoeff( nTemp[kNext], moments[kNext], rossecutoff );
	Double	coeffMinus = fSoot->GetSootRadRossCoeff( nTemp[kPrev], moments[kPrev], rossecutoff );
	Double	delq_R = 0.5 * coeff * rfCurr
	  * ( fWMinus[k] * nTemp[kPrev]
		 + fWCurr[k] * nTemp[kCurr]
		 + fWPlus[k] * nTemp[kNext] )
		+ 0.25 
		  * ( fFDWMinus[k] * rfPrev * coeffMinus
			 + fFDWCurr[k] * rfCurr * coeff 
			 + fFDWPlus[k] * rfNext * coeffPlus )
			+ 0.25 * rfCurr
			  * ( fFDWMinus[k] * coeffMinus 
				 + fFDWCurr[k] * coeff 
				 + fFDWPlus[k] * coeffPlus );
	return delq_R;
  }
  else {
	return fSoot->GetSootRadiation( nTemp[kCurr], moments[kCurr] );
  }
}

#endif // __TTRANS_FLAME_SOLVER_HPP__
