#ifndef TTRANS_FLAMELET_HPP__
#define TTRANS_FLAMELET_HPP__

#include "Spline.h"

// define WRITETURBX to write mole fraction in Scalar output file, undefine for mass fractions
#undef WRITETURBX

// define to write NO on dry basis at 6% oxygen into Scalar output file
#undef DRYNOX6PERCENT

// write BARLOW definition of mixture fraction in Scalar output file
#undef BARLOWZ

// define CALCZR to compute ZR = ZMEAN + 2 SQRT(ZVAR), with ZMEAN and ZVAR from CA.in file
#undef CALCZR

template<typename Species>
void TTransFlamelet<Species>::InitTTransFlamelet( void )
{
	int	i;
	fDeltaStepsOut = -1;
	fNOutputs = fInputData->fNOutputs;
	fNGridPoints = fInputData->fInitialGridPoints-2;
	fEquidistant = fInputData->fEquidistant;
	fTStart = fInputData->fStart;
	fTEnd = fInputData->fEnd;
	fDeltaTStart = 0.0;
	fPressStart = 50e5;
	fTimeCut = 1.3e-3;

	const char *inputFile;
	FILE *fpIn = nullptr;
	if (!fInputData->fAddFileNo1 ) {
		inputFile = NULL;
	}
	else {
		inputFile = GetFullPath( fInputData->fAddFileNo1, kFileName );
		fpIn = fopen(inputFile, "r");
	}
	char 	dummy[128];
	int		conv;
	Double	*timeIn = New1DArray( 500000 );
	Double	*pressureIn = New1DArray( 500000 );
	Double	*tOxIn = New1DArray( 500000 );
	Double	*tFuelIn = New1DArray( 500000 );
	Double	*chiIn = New1DArray( 500000 );
	Double	*zRIn = New1DArray( 500000 );
	Double	*zmeanIn = New1DArray( 500000 );
	Double	*zvarIn = New1DArray( 500000 );
	Double	*xOverDIn = New1DArray( 500000 );
	
	if ( fpIn ) {
		fprintf( stderr, "use input file '%s' for boundary conditions\n", inputFile );
		if ( fscanf( fpIn, "RPM = %lg", &fRPM ) <= 0 ) {
			fprintf( stderr, "#error: missing variable RPM\n" );
		}
		if ( fscanf( fpIn, " VarsIn = %d", &fVarsIn ) <= 0 ) {
			fprintf( stderr, "#error: missing variable VarsIn\n" );
			exit(2);
		}
		else {
			fprintf( stderr, "read %d vars from file %s\n", fVarsIn, inputFile );
		}
		
		int	j;
		if ( fVarsIn == 4 ) {
			conv = fscanf( fpIn, "%s%s%s%s", dummy, dummy, dummy, dummy );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "##error: in file %s trying to read %d variables\n", inputFile, fVarsIn );
				exit(2);
			}
			for ( j = 0; j < 500000; ++j ) {
				conv = fscanf( fpIn, "%lg%lg%lg%lg", &timeIn[j], &pressureIn[j], &tOxIn[j], &tFuelIn[j] );
				fprintf( stderr, "time = %g\n", timeIn[j] );
				if ( !conv || conv == EOF ) {
					break;
				}
			}
		}
		else if ( fVarsIn == 5 ) {
			conv = fscanf( fpIn, "%s%s%s%s%s", dummy, dummy, dummy, dummy, dummy );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "##error: in file %s trying to read %d variables\n", inputFile, fVarsIn );
				exit(2);
			}
			for ( j = 0; j < 500000; ++j ) {
				conv = fscanf( fpIn, "%lg%lg%lg%lg%lg", &timeIn[j], &pressureIn[j], &tOxIn[j], &tFuelIn[j], &chiIn[j] );
				if ( !conv || conv == EOF ) {
					break;
				}
			}
		}
		else if ( fVarsIn == 6 ) {
			conv = fscanf( fpIn, "%s%s%s%s%s%s", dummy, dummy, dummy, dummy, dummy, dummy );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "##error: in file %s trying to read %d variables\n", inputFile, fVarsIn );
				exit(2);
			}
			for ( j = 0; j < 500000; ++j ) {
				conv = fscanf( fpIn, "%lg%lg%lg%lg%lg%lg", &timeIn[j], &pressureIn[j], &tOxIn[j], &tFuelIn[j], &chiIn[j], &zRIn[j] );
				if ( !conv || conv == EOF ) {
					break;
				}
			}
		}
		else if ( fVarsIn == 8 ) {
			conv = fscanf( fpIn, "%s%s%s%s%s%s%s%s", dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "##error: in file %s trying to read %d variables\n", inputFile, fVarsIn );
				exit(2);
			}
			for ( j = 0; j < 500000; ++j ) {
				conv = fscanf( fpIn, "%lg%lg%lg%lg%lg%lg%lg%lg"
					, &timeIn[j], &pressureIn[j], &tOxIn[j]
					, &tFuelIn[j], &chiIn[j], &zRIn[j]
					, &zmeanIn[j], &zvarIn[j] );
				if ( !conv || conv == EOF ) {
					break;
				}
			}
		}
		else if ( fVarsIn == 9 ) {
			conv = fscanf( fpIn, "%s%s%s%s%s%s%s%s%s", dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "##error: in file %s trying to read %d variables\n", inputFile, fVarsIn );
				exit(2);
			}
			for ( j = 0; j < 500000; ++j ) {
				conv = fscanf( fpIn, "%lg%lg%lg%lg%lg%lg%lg%lg%lg"
					, &timeIn[j], &xOverDIn[j], &pressureIn[j], &tOxIn[j]
					, &tFuelIn[j], &chiIn[j], &zRIn[j]
					, &zmeanIn[j], &zvarIn[j] );
#ifdef CALCZR
				zRIn[j] = zmeanIn[j] + 2.0 * sqrt( zvarIn[j] );
#endif
				if ( !conv || conv == EOF ) {
					break;
				}
			}
		}
		else {
			fprintf( stderr, "##error: reading %d variables from file %s is not implemented\n"
						, fVarsIn, inputFile );
			exit(2);
		}
		
		fclose( fpIn );
	
		fTimeIn = NewVector( j );
		fPressureIn = NewVector( j );
		fTOxIn = NewVector( j );
		fTFuelIn = NewVector( j );
		if ( fVarsIn >= 5 ) {
			fChiIn = NewVector( j );
			copy( fChiIn->len, chiIn, 1, fChiIn->vec, 1 );
		}
		else {
			fChiIn = NULL;
		}
		if ( fVarsIn >= 6 ) {
			fZRIn = NewVector( j );
			copy( fZRIn->len, zRIn, 1, fZRIn->vec, 1 );
		}
		else {
			fZRIn = NULL;
		}
		if ( fVarsIn >= 8 ) {
			fZMeanIn = NewVector( j );
			fZVarIn = NewVector( j );
			copy( fZMeanIn->len, zmeanIn, 1, fZMeanIn->vec, 1 );
			copy( fZVarIn->len, zvarIn, 1, fZVarIn->vec, 1 );
		}
		else {
			fZMeanIn = NULL;
			fZVarIn = NULL;
		}
		
		if ( fVarsIn >= 9 ) {
			fXOverD = NewVector( j );
			copy( fXOverD->len, xOverDIn, 1, fXOverD->vec, 1 );
		}
		else {
			fXOverD = NULL;
		}
	
		copy( fPressureIn->len, pressureIn, 1, fPressureIn->vec, 1 );
		copy( fTOxIn->len, tOxIn, 1, fTOxIn->vec, 1 );
		copy( fTFuelIn->len, tFuelIn, 1, fTFuelIn->vec, 1 );
		if ( fRPM > 0 ) {
			for ( j = 0; j < fTimeIn->len; ++j ) {
				fTimeIn->vec[j] = ( timeIn[j] - timeIn[0] ) / 360.0 * 60 / fRPM;
			}
		} else {
			for ( j = 0; j < fTimeIn->len; ++j ) {
				fTimeIn->vec[j] = ( timeIn[j]/* - timeIn[0]*/ );
			}
		}
	
		Free1DArray( xOverDIn );
		Free1DArray( zvarIn );
		Free1DArray( zmeanIn );
		Free1DArray( zRIn );
		Free1DArray( chiIn );
		Free1DArray( tFuelIn );
		Free1DArray( tOxIn );
		Free1DArray( pressureIn );
		Free1DArray( timeIn );
		fUseInput = TRUE;
	}
	else {
		fUseInput = FALSE;
		fTimeIn = NULL;
		fPressureIn = NULL;
		fTOxIn = NULL;
		fTFuelIn = NULL;
		fChiIn = NULL;
		fZMeanIn = NULL;
		fZVarIn = NULL;
		fXOverD = NULL;
	}
	
	if ( fChiIn ) {
		cerr << "use transient history of scalar dissipation rate from input file" << NEWL;
		// compress data	
		Double	timeOut = -1.0;
		FILE	*fp = fopen( "CA.inCompress", "w" );
		fprintf( fp, "RPM = -1\nVarsIn = %d\nTime(s)", fVarsIn );
		if ( fVarsIn == 9 ) {
			fprintf( fp, "\tx/D" );
		}
		fprintf( fp, "\tPressure(Pa)\tTOx(K)\tTFuel(K)\tSci(1/s)" );
		if ( fVarsIn >= 6 ) {
			fprintf( fp, "\tZR" );
		}
		else if ( fVarsIn >= 8 ) {
			fprintf( fp, "\tZMean\tZVar" );
		}
		fprintf( fp, "\n" );
		for ( i = 0; i < fTimeIn->len; i+=1 ) {
			if ( i == 0 || i == fTimeIn->len-1 || fTimeIn->vec[i] >= timeOut + 5.0e-9 ) {
				fprintf( fp, "%g", fTimeIn->vec[i] );
				if ( fVarsIn == 9 ) {
					fprintf( fp, "\t%g", fXOverD->vec[i] );
				}
				fprintf( fp, "\t%g\t%g\t%g\t%g"
						, fPressureIn->vec[i], fTOxIn->vec[i], fTFuelIn->vec[i]
						, fChiIn->vec[i] );
				if ( fVarsIn >= 6 ) {
					fprintf( fp, "\t%g", fZRIn->vec[i] );
				}
				if ( fVarsIn >= 8 ) {
					fprintf( fp, "\t%g\t%g", fZMeanIn->vec[i], fZVarIn->vec[i] );
				}
				fprintf( fp, "\n" );
				timeOut = fTimeIn->vec[i];
			}
		}
		fclose( fp );
	}
	else {
		if ( fInputData->fParameterComm >= 0.0 ) {
			fScalarDissRate = fInputData->fParameterComm;
			cerr << "initial scalar dissipation rate is " << fScalarDissRate << NEWL;
		}
		else {
			if ( fInputData->fDissRate ) {
				fScalarDissRate = fInputData->fDissRate->vec[0];
				if ( fInputData->fDissRate->len > 1 ) {
					cerr << "#warning: more than one scalar dissipation rate specified, use" 
							<< fScalarDissRate << NEWL;
				}
				if ( fScalarDissRate < 0.0 ) {
					cerr << "use transient history of scalar dissipation rate" << NEWL;
				}
				else {
					cerr << "initial scalar dissipation rate is " << fScalarDissRate << NEWL;
				}
			}
			else {
				cerr << "#error: no scalar dissipation rate specified" << NEWL;
				exit( 2 );
			}
		}
	}
	
	fGridSol = NewVector( fNGridPoints+2 );
	fStartSol = NewMatrix( fVariables, fNGridPoints+2, kColumnPointers );
	
	fGridSol->vec = &fGridSol->vec[kNext];
	fStartSol->mat = &fStartSol->mat[kNext];

	fBCFlagLeft = New1DIntArray( fNOfEquations );
	fBCFlagRight = New1DIntArray( fNOfEquations );

	
	int	nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
	fVariableNames = new String[fVariablesWithoutSpecies + nOfSpeciesIn];
	if ( !fVariableNames ) {
		cerr << "#error memory allocation of TTransFlamelet failed" << NEWL;
		exit( 2 );
	}
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies.GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies.GetNames()[i] );
	}
	if ( fSoot ) {
		/* attention: counter  GetOffsetSootMoments() includes grid which is not included here */
		int	offset = fSoot->GetOffsetSootMoments()-1;
		for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
			fVariableNames[offset + i] = new char[8];
			sprintf( fVariableNames[offset + i], "M%d", i );
		}
	}

	// Double	*CheckGrid = &fSolGrid->vec[kPrev];

	int kkk;

	for ( kkk = -1; kkk < fNGridPoints+1; ++kkk ) {
			cerr << "fGridSol->vec[i] : " << fGridSol->vec[kkk] << NEWL;
		}

	MakeGrid( &fGridSol->vec[kPrev], fNGridPoints+2, 0.0, 1.0, fEquidistant );

	for ( kkk = -1; kkk < fNGridPoints+1; ++kkk ) {
		cerr << "fGridSol->vec[i] : " << fGridSol->vec[kkk] << NEWL;
	}


	SetInitialBC();

	if ( fUseInput ) {
		fStartSol->mat[kPrev][fTemperature] = fTOxIn->vec[0];
		fStartSol->mat[fNGridPoints][fTemperature] = fTFuelIn->vec[0];
	}

	ReadStartProfiles( fInputData );

	fTempOxEnd = fStartSol->mat[kPrev][fTemperature];
	fTempFuelEnd = fStartSol->mat[fNGridPoints][fTemperature];

	fTCurr = fTStart;
		
	if ( !fUseInput ) {
		if ( fNOutputs < 1 ) {
			fNOutputs = 1;
		}
		fDeltaT = ( fTEnd - fTStart ) / ( ( Double ) fNOutputs );
		if ( fDeltaT <= 0.0 ) {
			cerr << "#error: t_start = " << fTStart << " greater than t_end = " 
						<< fTEnd << NEWL;
			exit( 2 );
		}
	}

	fPDF = NULL;
	ffpEI = NULL;
	if ( fUseInput ) {
		fTimePDFIn = NULL;
		
		ReadPDF();
	}

	char	buffer[32];
	sprintf( buffer, "Scalars" );
	ffpEI = GetOutfile( buffer, FileType::kData );

	fprintf( ffpEI, "*\ntime [s]" );
	if ( fXOverD ) fprintf( ffpEI, "\tx/D" );
#ifdef BARLOWZ
	fprintf( ffpEI, "\tZBarlowMean" );
#endif
	for ( i = 0; i < GetNFuels(); ++i ) {
		fprintf( ffpEI, "\tEI_%s [kg/sm^3]", fSpecies.GetNames()[GetFuelIndex( i )] );
	}
	if ( fSpecies.FindSpecies( "NO" ) > -1 && !fSpecies.IsSteadyState(fSpecies.FindSpecies( "NO" )) ) {
		fprintf( ffpEI, "\tEI_NO [kg/sm^3]" );
		fprintf( ffpEI, "\tYMean_NO" );
	}
	if ( fSpecies.FindSpecies( "NO2" ) > -1 && !fSpecies.IsSteadyState(fSpecies.FindSpecies( "NO2" )) ) {
		fprintf( ffpEI, "\tEI_NO2 [kg/sm^3]" );

		fprintf( ffpEI, "\tYMean_NO2" );
	}
	fprintf( ffpEI, "\tCheckPDF" );
	if ( fSpecies.FindSpecies( "N2O" ) > -1 && !fSpecies.IsSteadyState(fSpecies.FindSpecies( "N2O" )) ) {
		fprintf( ffpEI, "\tEI_N2O [kg/sm^3]" );
	}
	if ( fSoot ) {
		fprintf( ffpEI, "\tEI_Soot [kg/sm^3]" );
	}
	if ( fZMeanIn ) {
		for ( i = 0; i < nOfSpeciesIn; ++i ) {
#ifdef WRITETURBX
			fprintf( ffpEI, "\tX_%sMean\tX_%sVar", fSpecies.GetNames()[i], fSpecies.GetNames()[i] );
#else
			fprintf( ffpEI, "\tY_%sMean\tY_%sVar", fSpecies.GetNames()[i], fSpecies.GetNames()[i] );
#endif
		}
		fprintf( ffpEI, "\tTempMean [K]\tTempVar [K^2]" );
		if ( fSoot ) {
			fprintf( ffpEI, "\tNumDens [kmole/m^3]\tNumDensVar [kmole^2/m^6]\tfv\tfvVar" );
		}
	}
#ifdef DRYNOX6PERCENT
	if ( fSpecies.FindSpecies( "NO" ) > -1 ) {
		fprintf( ffpEI, "\tNO_dry_6perc" );
	}
#endif
	fprintf( ffpEI, "\n" );
	fflush( ffpEI );

	FILE	*fp = GetOutfile( "InitialValues", FileType::kData );
	PrintSolution( fp, fNGridPoints, fNOfEquations, fGridSol->vec
					, fStartSol->mat, fVariableNames );
	fclose( fp );

#ifdef HP
	InstallInterruptHP();
#endif
}

template<typename Species>
TTransFlamelet<Species>::~TTransFlamelet( void )
{
	if ( fPDF ) {
		for ( int i = 0; i < fPDF->rows; ++i ) {
			fPDF->mat[i] = &fPDF->mat[i][kPrev];
		}
		DisposeMatrix( fPDF );
		DisposeVector( fTimePDFIn );
	}

	Free1DIntArray( fBCFlagRight );
	Free1DIntArray( fBCFlagLeft );

	fStartSol->mat = &fStartSol->mat[kPrev];
	fGridSol->vec = &fGridSol->vec[kPrev];

	DisposeMatrix( fStartSol );
	DisposeVector( fGridSol );
	
	if ( fUseInput ) {
		DisposeVector( fTFuelIn );
		DisposeVector( fTOxIn );
		
		DisposeVector( fPressureIn );
		DisposeVector( fTimeIn );
	}
	
	if ( ffpEI ) {
		fclose( ffpEI );
	}
}

template<typename Species>
void TTransFlamelet<Species>::Solve( void )
{
	FILE	*fp;
	Flag 	leave = FALSE;
	Double	bound = fTEnd - 1.0e-20;
	Double	tempEnd;
	int 	start = 0;
	
	if ( fUseInput ) {
		if ( fTStart > 0.0 ) {
			while ( fTimeIn->vec[start] < fTStart + 1.0e-9 ) ++start;
		}
		if ( fTStart < fTimeIn->vec[start] ) {
			Initialize( fTStart, fVariableNames, &fStartSol->mat[kPrev], fNOfEquations
							, &fGridSol->vec[kPrev], fNGridPoints+2
							, fPressStart, GetScalarDiss( fTStart ), fDeltaTStart, -1, 1 );
			--start;
		}
		else {        // fTStart = fTimeIn->vec[start]
			Initialize( fTimeIn->vec[start], fVariableNames, &fStartSol->mat[kPrev]
							, fNOfEquations, &fGridSol->vec[kPrev], fNGridPoints+2
							, fPressureIn->vec[start], GetScalarDiss( fTimeIn->vec[start] )
							, fDeltaTStart, -1, 1 );
		}
	}
	else {
		Initialize( fTStart, fVariableNames, &fStartSol->mat[kPrev], fNOfEquations
						, &fGridSol->vec[kPrev], fNGridPoints+2
						, TFlame<Species>::GetPressure(), GetScalarDiss( fTStart ), fDeltaTStart, -1, 1 );
	}
 
	cerr << "dump output" << NEWL;
	// fp = GetOutputFile( fTCurr, NULL, NULL, FileType::kText );
	// WriteFlameletFileInitialize( fp, NULL, NULL );
	// fclose( fp );
	WriteFlameletFileInitialize( NULL, NULL, NULL );  // Replaced the above by Rui


	cerr << NEWL << "start computation" << NEWL << NEWL;

	if ( fUseInput ) {
		for ( int i = start+1; i < fTimeIn->len; ++i ) {
			fTCurr = fTimeIn->vec[i];
			Double	diss = GetScalarDiss( fTCurr );
			fprintf( stderr, "\ntend = %g", fTCurr );
			if ( fXOverD ) {
				fprintf( stderr, "\tx/D = %g", fXOverD->vec[i] );
			}
			fprintf( stderr, "\tPend = %g\nToxend = %g\tTfEnd = %g\tchi = %g\n"
					, fPressureIn->vec[i]
					, fTOxIn->vec[i]
					, fTFuelIn->vec[i]
					, diss );
			leave = TTransFlameSolver<Species>::Solve( fTCurr, fPressureIn->vec[i]
				, diss, fTOxIn->vec[i], fTFuelIn->vec[i]
				, ( fZRIn ) ? fZRIn->vec[i] : 1.0, fNOutputs );
			fprintf( stderr, "solved\n" );
		
			cerr << "dump output" << NEWL;
		
			// fp = GetOutputFile( fTCurr, NULL, NULL, FileType::kText );                        // Rui
			// if ( !fp ) {
			// 	fprintf( stderr, "#error: can't open outputfile\n" );
			// 	exit( 2 );
			// }
			// WriteFlameletFile( fp, NULL, NULL );
			// fclose( fp );
			WriteScalars( fTCurr, i );
		}
	}
	else {
		int	jCountOutputs = 0;
		do {
			fTCurr += fDeltaT;
			++jCountOutputs;
			tempEnd = fTempFuelEnd;
			leave = TTransFlameSolver<Species>::Solve( fTCurr, TFlame<Species>::GetPressure(), GetScalarDiss( fTCurr ), fTempOxEnd, tempEnd, fDeltaStepsOut );

			if ( jCountOutputs % 1 == 0 ) {
				cerr << "dump output" << NEWL;
				
				// fp = GetOutputFile( fTCurr, ( leave ) ? "Sol" : NULL, NULL, FileType::kText );    //Rui
				// WriteFlameletFile( fp, NULL, NULL );
				// fclose( fp );
			}
			WriteScalars( fTCurr, -1 );
	
		} while ( fTCurr < bound && !leave );
	}

	cerr << "done" << NEWL;

	return;
}

template<typename Species>
Double TTransFlamelet<Species>::GetScalarDiss( Double theTime )
{
	if ( fChiIn ) {
		int		i;
		Double	*timeIn = fTimeIn->vec;
		
		for ( i = 1; i < fTimeIn->len; ++i ) {
			if ( theTime <= timeIn[i] ) {
				break;
			}
		}
		Double	chiOld = fChiIn->vec[i-1];
		Double	tOld = timeIn[i-1];
		Double	chiNew = fChiIn->vec[i];
		Double	tNew = timeIn[i];
		
		return Interpol( theTime, chiOld, tOld, chiNew, tNew );
	}
	else {
		return fScalarDissRate;
	}
}

template<typename Species>
Double TTransFlamelet<Species>::GetTempEnd( Double theTime, Double timeStart, Double timeCut
									, Double tempStart, Double tempEnd )
{
	if ( theTime < timeCut ) {
		return Interpol( theTime, tempStart, timeStart, tempEnd, timeCut );
	}
	else {
		return tempEnd;
	}
}

template<typename Species>
void TTransFlamelet<Species>::SetInitialBC( void )
{
	int					i;
	Double				mixMolarMass;
	int					nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	//SpeciesPtr			species = fInputData->GetSpecies();				// ###
	BoundaryInputPtr	right = fInputData->leftBoundary;
	BoundaryInputPtr	left = fInputData->rightBoundary;
	int					inpTOffset = fInputData->fTemperatureOffset;
	int					*speciesIndexLeft = NULL;
	int					*speciesIndexRight = NULL;
	int					leftSpecifiedSpecies = left->fSpecifiedSpeciesBCs;
	int					rightSpecifiedSpecies = right->fSpecifiedSpeciesBCs;
	Double				*yleft = fStartSol->mat[kPrev];
	Double				*yright = fStartSol->mat[fNGridPoints];
	
	speciesIndexLeft = new int[left->fSpecifiedSpeciesBCs];
	if ( !speciesIndexLeft ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );
	speciesIndexRight = new int[right->fSpecifiedSpeciesBCs];
	if ( !speciesIndexRight ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );

	//	set speciesIndex
	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		speciesIndexLeft[i] = fInputData->FindSpecies( left->speciesName[i] );
	}
	for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
		speciesIndexRight[i] = fInputData->FindSpecies( right->speciesName[i] );
	}
	
	// set BCFlags
	fBCFlagLeft[fTemperature] = left->fBcFlag[inpTOffset];
	fBCFlagRight[fTemperature] = right->fBcFlag[inpTOffset];
	
	if ( left->fBcFlagSpecies == kNeumann || right->fBcFlagSpecies == kNeumann ) {
		cerr << "#error: invalid boundary condition type "
			<< "'Neumann' for species concentration" << NEWL;			
	}
	for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
		fBCFlagLeft[i] = left->fBcFlagSpecies;
		fBCFlagRight[i] = right->fBcFlagSpecies;
	}

	// set value
	yleft[fTemperature] = left->fValue[inpTOffset];
	yright[fTemperature] = right->fValue[inpTOffset];
	
	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		yleft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
	}
	if ( left->fMixtureSpecification == kMolarFraction ) {
		// first compute molar mass of mixture
		for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
			mixMolarMass += fSpecies.GetMolarMass()->vec[i] * yleft[i+fFirstSpecies];	// ###
		}
		// compute massfractions
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			yleft[i+fFirstSpecies] *= fSpecies.GetMolarMass()->vec[i] / mixMolarMass;	// ###
		}
		for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
			fBCFlagLeft[i] = kMassFraction;
		}
	}

	for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
		yright[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
	}
	if ( right->fMixtureSpecification == kMolarFraction ) {
		// first compute molar mass of mixture
		for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
			mixMolarMass += fSpecies.GetMolarMass()->vec[i] * yright[i+fFirstSpecies];	// ###
		}
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			yright[i+fFirstSpecies] *= fSpecies.GetMolarMass()->vec[i] / mixMolarMass;	// ###
		}
		for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
			fBCFlagRight[i] = kMassFraction;
		}
	}
	delete speciesIndexRight;
	delete speciesIndexLeft;
}

template<typename Species>
void TTransFlamelet<Species>::SetInitialValues( void )
{
	Double	**y = fStartSol->mat;
	Double	*x = fGridSol->vec;
	Double	*yLeft = fStartSol->mat[kPrev];
	Double	*yRight = fStartSol->mat[fNGridPoints];
	Double	left = fGridSol->vec[kPrev];
	Double	right = fGridSol->vec[fNGridPoints];

	for ( int k = 0; k < fNGridPoints; ++k ) {
		for ( int j = 0; j < fNOfEquations; ++j ) {
			y[k][j] = Interpol( x[k], yLeft[j], left, yRight[j], right );
		}
	}
}

template<typename Species>
void TTransFlamelet<Species>::ReadStartProfiles( TInputDataPtr inp )
{
	StartProfilePtr	sp = NULL;
	char 	*insp = inp->fStartProfileFile;
	FILE	*fpS = NULL;
	char	*fileName = GetFullPath( insp, kFileName );

	sp = new StartProfile;	

	if ( !insp || ( fpS = fopen( fileName, "r" ) ) == NULL ) {
		fprintf( stderr, "interpolate initial solution from boundary conditions\n" );
		SetInitialValues();
	} 
	else {
		if ( ( fpS = fopen( fileName, "r" ) ) == NULL ) {
			fprintf( stderr, "startprofiles file '%s' not found\n", fileName );
			exit( 2 );
		} 
		else {
			fprintf( stderr, "read initial solution from file '%s'\n", fileName );
			::ReadStartProfiles( sp, fpS );
			SetInitialValues( fInputData, sp );
			CleanReadStartProfile();
			delete sp;
			fclose( fpS );
		}
	}
	delete fileName;
	
}

template<typename Species>
void TTransFlamelet<Species>::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
{
	int 				i, j, k;
	int					variable, speciesIndex;
	Flag				ZSet = FALSE;
	Flag				chooseInputGrid = FALSE;
	int					gridPointsIn = sp->gridPoints;	// including boundaries
	int					nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	Double				*xIn = New1DArray( gridPointsIn );
	Double				*yIn =  New1DArray( gridPointsIn );
	Double				*yInFloat = sp->data;
	char				*string = sp->labels;
	int					oxidizerSide; // this program assumes that oxidizerSide = kLeft
	Double				**y = fStartSol->mat;
	Double				*yWork = New1DArray( fNGridPoints+2 );
	yWork = &yWork[kNext];
	SplinePtr			theSpline = NULL;
	Double				leftSlope;
	Double				rightSlope;
	Double				*x = fGridSol->vec;
	struct _parameter	*param = GetParameter( "pressure" );
	
	SetInitialValues(); // default values

//	choose grid from input or equidistant
	if ( gridPointsIn <= fNGridPoints+2 ) {
		fNGridPoints = gridPointsIn - 2;
		chooseInputGrid = TRUE;
	}
	else {
		chooseInputGrid = FALSE;
	}
	
// find independent coordinate
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "z", 1 ) == 0 && ZSet == FALSE ) {
			if ( yInFloat[i * gridPointsIn] < yInFloat[(i+1) * gridPointsIn - 1] ) {
				oxidizerSide = kLeft;
			}
			else {
				oxidizerSide = kRight;
			}
			if ( chooseInputGrid ) {
				cerr << "choose inputGrid" << NEWL;
				if ( oxidizerSide == kRight ) { // turn vector
					for ( j = -1; j <= gridPointsIn-2; ++j ) {
						x[gridPointsIn-j-3] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
					}
				}
				else { // copy vector
					for ( j = -1; j <= gridPointsIn-2; ++j ) {
						x[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
					}
				}
				ZSet = TRUE;
			}
			else { // choose own Grid, but read x for interpolation
				cerr << "choose own Grid" << NEWL;
				if ( oxidizerSide == kRight ) { // turn vector
					for ( j = 0; j < gridPointsIn; ++j ) {
						xIn[gridPointsIn - j - 1] = yInFloat[i*gridPointsIn + j];		// implicit cast from float to Double
					}
				}
				else { // copy vector
					for ( j = 0; j < gridPointsIn; ++j ) {
						xIn[j] = yInFloat[i*gridPointsIn + j];		// implicit cast from float to Double
					}
				}
				ZSet = TRUE;
			}
		}
		string += strlen(string) + 1;
	}

// error checking
	if ( !ZSet ) {
		cerr << "error: can't find coordinate 'Z'" << NEWL;
		exit(2);
	}

// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = fTemperature;
		}
		else if ( strncmp( string, "massfraction-", 13 ) == 0 ){
			string += 13;
			UpperString( string );
			if ( ( speciesIndex = inp->FindSpecies( string ) ) >= 0 ) {
				if ( speciesIndex < nSpeciesInSystem ) {
					variable = fFirstSpecies + speciesIndex;
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
		else if ( GetSoot() && strncmp( string, "m", 1 ) == 0 ) {
			string += 1;
			int	num = atoi( string );
			if ( isdigit( string[0] ) && num < GetSoot()->GetNSootMoments() ) {
				/* attention: counter  GetOffsetSootMoments() includes grid which is not included here */
				variable = num + GetSoot()->GetOffsetSootMoments() -1; 
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
		if ( chooseInputGrid ) {
			if ( oxidizerSide == kRight ) { // turn vector
				for ( k = -1; k <= gridPointsIn-2; ++k ) {
					y[k][variable] = yInFloat[(i+1)*gridPointsIn - k-2];	// copy workspace to vector of solution
				}
			}
			else {
				for ( k = -1; k <= gridPointsIn-2; ++k ) {
					y[k][variable] = yInFloat[i*gridPointsIn + k + 1];	// copy workspace to vector of solution
				}
			}
		}
		else {
			for ( k = 0; k < gridPointsIn; ++k ) {	// store vector in workspace
				yIn[k] = yInFloat[i * gridPointsIn + k];	// implicit cast from float to Double
			}
		
			leftSlope = ( yIn[1] - yIn[0] ) / ( xIn[1] - xIn[0] );
			rightSlope = ( yIn[gridPointsIn-1] - yIn[gridPointsIn-2] ) / ( xIn[gridPointsIn-1] - xIn[gridPointsIn-2] );
			theSpline = ComputeSimpleSpline( xIn, yIn, gridPointsIn, FALSE, leftSlope, FALSE, rightSlope, NULL, TRUE );
			SplineInterpolate( theSpline, &x[kPrev], &yWork[kPrev], fNGridPoints+2 );
			if ( oxidizerSide == kRight ) { // turn vector
				for ( k = -1; k <= fNGridPoints; ++k ) {
					y[k][variable] = yWork[fNGridPoints-k-1];	// copy workspace to vector of solution
				}
			}
			else {
				for ( k = -1; k <= fNGridPoints; ++k ) {
					y[k][variable] = yWork[k];	// copy workspace to vector of solution
				}
			}
		}
	}
	
	FreeSpline( theSpline );
	yWork = &yWork[kPrev];
	Free1DArray( yWork );
	Free1DArray( yIn );
	Free1DArray( xIn );
	
// set scalars
	if ( param ) {
		fPressStart = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
			fPressStart *= 1.0e5;
		}
	}
	else { // choose default
		cerr << "#warning: no input value for 'pressure', use " << fPressStart << NEWL;
	}
	

	if ( fInputData->fParameterComm < 0.0 ) {
	  param = GetParameter( "chi" );
	  if ( param ) {
	    fScalarDissRate = (Double)param->what.quantity.value;
	  }
	  else {
	    param = GetParameter( "chi_ref" );
	    if ( param ) {
	      fScalarDissRate = (Double)param->what.quantity.value;
	    }
	  }
	}

	param = GetParameter( "time" );
	if ( param ) {
		fTStart = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			fTStart *= 1.0e-3;
		}
	}
	else { // choose default
		cerr << "#warning: no input value for 'time', use " << fTStart << NEWL;
	}

	param = GetParameter( "currtimestep" );
	if ( param ) {
		fDeltaTStart = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			fDeltaTStart *= 1.0e-3;
		}
	}
	else { // choose default
		cerr << "#warning: no input value for 'CurrTimeStep', use default" << NEWL;
	}
	
	fTempOxStart = y[kPrev][fTemperature];
	fTempFuelStart = y[fNGridPoints][fTemperature];


	FILE	*fp = GetOutfile( "initialguess", FileType::kData );
	PrintSolution( fp, fNGridPoints, fNOfEquations, fGridSol->vec, y, fVariableNames );
	fclose(fp);
}

template<typename Species>
void TTransFlamelet<Species>::ReadPDF( void )
{
	const char	*inputFile = "PDF.in";
	Double 	ddummy;
	int		conv;
	int		i, j;
	Double	*timeIn = New1DArray( 5000 );
	Double	**Z;
	Double	**PDF;
	FILE 	*fpIn = fopen( inputFile, "r" );

	if ( fpIn ) {
		fprintf( stderr, "use input file '%s' for PDF\n", inputFile );
		if ( fscanf( fpIn, "gridpoints = %d", &fPDFGridPoints ) <= 0 ) {
			fprintf( stderr, "#error: missing variable PDFGridPoints\n" );
			exit(2);
		}
		Z = New2DArray( 5000, fPDFGridPoints );
		PDF = New2DArray( 5000, fPDFGridPoints );

		Double	*Zj;
		Double	*PDFj;
		for ( j = 0; j < 5000; ++j ) {
			conv = fscanf( fpIn, " TIME = %lg", &timeIn[j] );
			if ( !conv || conv == EOF ) {
				break;
			}
			conv = fscanf( fpIn, " CRANK %lg", &ddummy );
			conv = fscanf( fpIn, "%*s%*s" );
			Zj = Z[j];
			PDFj = PDF[j];
			for ( i = 0; i < fPDFGridPoints; ++i ) {
				conv = fscanf( fpIn, "%lg%lg%lg", &Zj[i], &ddummy, &PDFj[i] );
			}
			for ( i = fPDFGridPoints-1; i > 0; --i ) {
				Zj[i] = 0.5 * ( Zj[i] + Zj[i-1] );
				PDFj[i] *= fPDFGridPoints;
			}
			Zj[0] = 0.5 * ( Zj[i] + 0 );
			PDFj[0] *= fPDFGridPoints;
		}
		fclose( fpIn );
		fprintf( stderr, "%d time records read for pdf\n", j );
	
		fTimePDFIn = NewVector( j );
		fPDF = NewMatrix( j, fNGridPoints+2, kRowPointers );
	

		Double	**pdf = fPDF->mat;
		for ( i = 0; i < j; ++i ) {
			pdf[i] = &pdf[i][kNext];
			LinearInterpolate( Z[i], PDF[i], fPDFGridPoints,
						&fGridSol->vec[kPrev], &pdf[i][kPrev], fNGridPoints+2 );
		}

		copy( fTimePDFIn->len, timeIn, 1, fTimePDFIn->vec, 1 );
	
		FILE *fpPDF = GetOutfile( "PDF", FileType::kData );
		fprintf( fpPDF, "*\nZ" );
		for ( i = 0; i < fTimePDFIn->len; ++i ) {
			fprintf( fpPDF, "\t%g", fTimePDFIn->vec[i] );
		}
		for ( j = -1; j < fNGridPoints+1; ++j ) {
			fprintf( fpPDF, "\n%g", fGridSol->vec[j] );
			for ( i = 0; i < fTimePDFIn->len; ++i ) {
				fprintf( fpPDF, "\t%g", pdf[i][j] );
			}
		}
		fclose( fpPDF );

		Free2DArray( PDF );
		Free2DArray( Z );
		Free1DArray( timeIn );
	}
	else {
		fprintf( stderr, "#warning: no input file for PDF\n" );
	}
}

template<typename Species>
void TTransFlamelet<Species>::WriteScalars( Double time, int ind )
{
	int		i, index, indexNO, indexNO2;
	int		nOfSpeciesIn;
	Double	**pdf;
	Double	*tpdf;

	if ( fPDF ) {
		pdf = fPDF->mat;
		tpdf = fTimePDFIn->vec;
	}
	else {
		pdf = NULL;
		tpdf = NULL;
	}

	fprintf( ffpEI, "%g", time );

	if ( fXOverD ) fprintf( ffpEI, "\t%g", fXOverD->vec[ind] );

#ifdef BARLOWZ
	if ( fZMeanIn ) {
		fprintf( ffpEI, "\t%g", TurbMeanZBarlow( time, fZMeanIn->vec[ind], fZVarIn->vec[ind] ) );
	}
#endif

	Double	EIFuel = 0.0;
	for ( i = 0; i < GetNFuels(); ++i ) {
		EIFuel = ComputeEmissionIndex( time, GetFuelIndex( i ), pdf, tpdf );
		fprintf( ffpEI, "\t%g", EIFuel );
	}
	
	indexNO = fSpecies.FindSpecies( "NO" );
	if ( indexNO > -1 && !fSpecies.IsSteadyState(indexNO) ) {
		fprintf( ffpEI, "\t%g", ComputeEmissionIndex( time, indexNO, pdf, tpdf ) );
		fprintf( ffpEI, "\t%g", ComputeMeanFromPDF( time, indexNO, pdf, tpdf ) );
	}
	
	indexNO2 = fSpecies.FindSpecies( "NO2" );
	if ( indexNO2 > -1 && !fSpecies.IsSteadyState(indexNO2) ) {
		fprintf( ffpEI, "\t%g", ComputeEmissionIndex( time, indexNO2, pdf, tpdf ) );
		fprintf( ffpEI, "\t%g", ComputeMeanFromPDF( time, indexNO2, pdf, tpdf ) );
	}
	fprintf( ffpEI, "\t%g", CheckPDF( time, pdf, tpdf ) );
	index = fSpecies.FindSpecies( "N2O" );
	if ( index > -1 && !fSpecies.IsSteadyState(index) ) {
		fprintf( ffpEI, "\t%g", ComputeEmissionIndex( time, index, pdf, tpdf ) );
	}

	if ( fSoot ) {
		fprintf( ffpEI, "\t%g", ComputeEmissionIndexSoot( time, 1, pdf, tpdf ) * 24 );
	}

	if ( fZMeanIn ) {
		Double	mean, variance;
		nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
		for ( i = 0; i < nOfSpeciesIn; ++i ) {
#ifdef WRITETURBX
			mean = TurbMeanX( time, i, fZMeanIn->vec[ind], fZVarIn->vec[ind], &variance );
#else
			mean = TurbMeanY( time, i, fZMeanIn->vec[ind], fZVarIn->vec[ind], &variance );
#endif
			fprintf( ffpEI, "\t%g\t%g", mean, variance );
		}
		mean = TurbMeanTemp( time, fZMeanIn->vec[ind], fZVarIn->vec[ind], &variance );
		fprintf( ffpEI, "\t%g\t%g", mean, variance );
		if ( fSoot ) {
			mean = TurbMeanSoot( time, 0, fZMeanIn->vec[ind], fZVarIn->vec[ind], &variance );
			fprintf( ffpEI, "\t%g\t%g", mean, variance );
			mean = TurbMeanSoot( time, 1, fZMeanIn->vec[ind], fZVarIn->vec[ind], &variance );
			fprintf( ffpEI, "\t%g\t%g", mean / 75.0, variance / ( 75.0 * 75.0 ) );
		}
#ifdef DRYNOX6PERCENT
		indexNO = fSpecies.FindSpecies( "NO" );
		Double H2O = TurbMeanX( time, fSpecies.FindSpecies( "H2O" ) , fZMeanIn->vec[ind], fZVarIn->vec[ind], NULL );
		Double XNO_dry = TurbMeanX( time, indexNO, fZMeanIn->vec[ind], fZVarIn->vec[ind], NULL ) / (1.0 - H2O);
		Double XO2_dry = TurbMeanX( time, fSpecies.FindSpecies( "O2" ) , fZMeanIn->vec[ind], fZVarIn->vec[ind], NULL ) / (1.0 - H2O);
		Double NODry6 = 1.0e6*XNO_dry*(4.76*(3+(1-XO2_dry)*2)/(1-4.76*XO2_dry)-2)/(4.76*(3+(1-0.06)*2)/(1-4.76*0.06)-2);
		if ( indexNO > -1 ) {
			fprintf( ffpEI, "\t%g", NODry6 );
			fprintf( stderr, "H2O = %g\tNO = %g\tO2 = %g\tNO6 = %g\n", H2O, XNO_dry, XO2_dry, NODry6 );
		}
#endif
	}

	fprintf( ffpEI, "\n" );
	fflush( ffpEI );

	static int	counter = 0;
	char 	dummy[128];
	Double	xOverDs[100];
	int		j, conv;
	Double	radius;
	Double	ZM;
	Double	ZVar;

// *** CUSTOMIZE HERE
/*  Liangyu */
	xOverDs[0] = .98;
	xOverDs[1] = 1.16;
	xOverDs[2] = 1.42;
	xOverDs[3] = 100000;

	if ( !fXOverD ) {
		return;
	}

fprintf( stderr, "start radial\n" );

	if ( fXOverD->vec[ind] > xOverDs[counter] ) {
		fprintf(stderr, "# print radial profiles at %g\n", xOverDs[counter] );
		char	buffer[32];
// *** CUSTOMIZE HERE
		sprintf( buffer, "rad_profile_%02d.dout", (counter+1) );

		fprintf( stderr, "open %s\n", buffer );
		FILE *fpIn = fopen( buffer, "r" );

		if ( !fpIn ) {
			fprintf( stderr, "cannot open file '%s'\n", buffer );
			return;
		}
		else {
			if (fscanf( fpIn, "%s", dummy ) == 0) 
		      printf("Could not write to file\n"); 
			if (fscanf( fpIn, "%s%s%s", dummy, dummy, dummy ))
			  printf("Could not write to file\n"); 
		}

		sprintf( buffer, "RadMean_%3.1f", xOverDs[counter] );
		FILE *fpOutMean = GetOutfile( buffer, FileType::kData );
		sprintf( buffer, "RadVar_%3.1f", xOverDs[counter] );
		FILE *fpOutVar = GetOutfile( buffer, FileType::kData );

// *** CUSTOMIZE HERE
		fprintf( fpOutMean, "*\nr/D" );
#ifdef BARLOWZ
		fprintf( fpOutMean, "\tZBarlowMean" );
		fprintf( fpOutVar, "\tZBarlowMean" );
#endif
		for ( i = 0; i < nOfSpeciesIn; ++i ) {
#ifdef WRITETURBX
			fprintf( fpOutMean, "\tX_%s", fSpecies.GetNames()[i] );
#else
			fprintf( fpOutMean, "\tX_%s", fSpecies.GetNames()[i] );
#endif
		}
		fprintf( fpOutMean, "\tTemperature [K]" );
		if ( fSoot ) {
			fprintf( fpOutMean, "\tNumDens [kmole/m^3]\tfv" );
		}
		fprintf( fpOutMean, "\tTotalEnergy [J/kg]" );
#ifdef DRYNOX6PERCENT
		fprintf( fpOutMean, "\tDryNO_6 [ppm]" );
#endif
		fprintf( fpOutMean, "\n" );

// *** CUSTOMIZE HERE
		fprintf( fpOutVar, "*\nr/D" );
		for ( i = 0; i < nOfSpeciesIn; ++i ) {
#ifdef WRITETURBX
			fprintf( fpOutVar, "\tY_%s", fSpecies.GetNames()[i] );
#else
			fprintf( fpOutVar, "\tX_%s", fSpecies.GetNames()[i] );
#endif
		}
		fprintf( fpOutVar, "\tTemperature [K]" );
		if ( fSoot ) {
			fprintf( fpOutVar, "\tNumDens [kmole^2/m^6]\tfv" );
		}
		fprintf( fpOutVar, "\n" );

		for ( j = 0; j < 5000; ++j ) {
			conv = fscanf( fpIn, "%lg%lg%lg", &radius, &ZM, &ZVar );
			if ( !conv || conv == EOF ) {
				break;
			}
			else {
				fprintf( fpOutMean, "%g", radius );
				fprintf( fpOutVar, "%g", radius );
				Double	variance;
#ifdef BARLOWZ
					Double	meanZBarlow;
					meanZBarlow = TurbMeanZBarlow( time, ZM, ZVar );
					fprintf( fpOutMean, "\t%g", meanZBarlow );
					fprintf( fpOutVar, "\t%g", meanZBarlow );
#endif
				for ( i = 0; i < nOfSpeciesIn; ++i ) {
#ifdef WRITETURBX
					fprintf( fpOutMean, "\t%g", TurbMeanX( time, i, ZM, ZVar, &variance ) );
#else
					fprintf( fpOutMean, "\t%g", TurbMeanY( time, i, ZM, ZVar, &variance ) );
#endif
					fprintf( fpOutVar, "\t%g", variance );
				}
				fprintf( fpOutMean, "\t%g", TurbMeanTemp( time, ZM, ZVar, &variance ) );
				fprintf( fpOutVar, "\t%g", variance );
				if ( fSoot ) {
					fprintf( fpOutMean, "\t%g", TurbMeanSoot( time, 0, ZM, ZVar, &variance ) );
					fprintf( fpOutVar, "\t%g", variance );
					fprintf( fpOutMean, "\t%g", TurbMeanSoot( time, 1, ZM, ZVar, &variance ) / 75.0 );
					fprintf( fpOutVar, "\t%g", variance / ( 75.0 * 75.0 ) );
				}
				fprintf( fpOutMean, "\t%g", TurbMeanTotEnergy( time, ZM, ZVar ) );
#ifdef DRYNOX6PERCENT
				Double H2O = TurbMeanX( time, fSpecies.FindSpecies( "H2O" ) , ZM, ZVar, NULL );
				Double XNO_dry = TurbMeanX( time, fSpecies.FindSpecies( "NO" ) , ZM, ZVar, NULL ) / (1.0 - H2O);
				Double XO2_dry = TurbMeanX( time, fSpecies.FindSpecies( "O2" ) , ZM, ZVar, NULL ) / (1.0 - H2O);
				Double NODry6 = 1.0e6*XNO_dry*(4.76*(3+(1-XO2_dry)*2)/(1-4.76*XO2_dry)-2)/(4.76*(3+(1-0.06)*2)/(1-4.76*0.06)-2);
				fprintf( fpOutMean, "\t%g", NODry6 );
#endif
				fprintf( fpOutMean, "\n" );
				fprintf( fpOutVar, "\n" );
			}
		}
		fclose( fpOutVar );
		fclose( fpOutMean );
		fclose( fpIn );
		++counter;
	}
}

#endif // TTRANS_FLAMELET_HPP__
