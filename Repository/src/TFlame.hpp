#ifndef TFLAME_HPP__
#define TFLAME_HPP__

#include <vector>

#include "MapMan.h"
#include "FlameUtilities.h"
#include "RadiationModel.h"


// if source code for production rates is generated with ScanMan, it can be tested
// here and compared with the solution of the mechanism itself
// For using production rate files produced by ScanMan, define PRODRATEFILE
// If Fortran files should be used, define PRODRATEFILEF77, for C files undefine PRODRATEFILEF77
// Needs to be set also in TSpecies.C
#undef PRODRATEFILE
#undef CHECKPRODRATE

#define MOLARDIFFUSION

// optimize
#define DELTAINEW
#define OPTIMIZEDELTAI

#ifdef PRODRATEFILE
void ComputeTheProductionRates( Double *productionRate, Double *reactionRate
				, Double temp, Double pressure, Double *concs );
#endif

template<class Species>
void TFlame<Species>::InitTFlame( const FirstInput* firstInp )
{
   	if ( !( fInputData = new TInputData( firstInp ) ) ) FatalError( "memory allocation of TInputData failed" );

	fClipNegativeConcs = fInputData->fClipNegativeConcs;
	
	fFuelIndex = NewIntVector( fInputData->fFuelIndex->len );
	for ( int i = 0; i < fFuelIndex->len; ++i ) {
		fFuelIndex->vec[i] = fInputData->fFuelIndex->vec[i];
	}
	fFromSpeciesIndex = fInputData->fFromSpeciesIndex;
	fToSpeciesIndex = fInputData->fToSpeciesIndex;
	fDateCreated = fInputData->fDateCreated;
	fAuthor = fInputData->fAuthor;

	fContinType = fInputData->fContinType;
	fContBound = fInputData->fContBound;
	fContInc = fInputData->fContInc;
	
	fOutputPath = GetFullPath( fInputData->fOutputPath, kPath );
	fprintf( stderr, "%s%s%s\n", "outpath is '", fOutputPath, "'"  );
	fOutFile = new char[strlen( fOutputPath ) + 128];


    // radiation

    fRadiationName = fInputData->fRadiationName;
    fRadiativeFrac = fInputData->fRadiativeFrac;

    // WSGG

    cout << " Radiation Modeling / Model " << fRadiationName << " Radiative Frac " << fRadiativeFrac << endl;

	if ( fInputData->fPressure ) {
		fPressure = NewVector( fInputData->fPressure->len );
		copy_vec( fPressure, fInputData->fPressure );
		fPressure->len = 0;
		if ( GetPressure() > 0.0 ) {
			fprintf( stderr, "%s%g%s\n", "initial pressure is ", GetPressure()/1.0e5, " bar"  );
		}
	}
	else {
		fPressure = NewVector( 1 );
		fPressure->vec[0] = -1.0;
		fPressure->len = 0;
	}

	fUseNumericalJac = fInputData->fUseNumericalJac;
	fUseNumericalDM = fInputData->fUseNumericalDM;

	// copy objects for sensitivity analysis
	fSensObjAll = fInputData->fSensObjAll;
	fSensAnal = fInputData->fSensAnal;
	fSensAnalSpec = fInputData->fSensAnalSpec;
	fSensAnalReac = fInputData->fSensAnalReac;
	//cai: 24/08/2012
	fSensAnalClas = fInputData->fSensAnalClas;
	fSensAnalFac = fInputData->fSensAnalFac;
	fSensMax = fInputData->fSensMax;
	fSensFinal = fInputData->fSensFinal;

	if (!fSensObjAll)
	  {
	    fNSensObj = fInputData->fNSensObj;
	    fSensObj = new char *[fNSensObj];
	    CopyStringArray (fSensObj, fInputData->fSensObj, fNSensObj);
	  }

	fNSensObj = fInputData->fNSensObj;
	fSensObj = new char* [fNSensObj];
	CopyStringArray( fSensObj, fInputData->fSensObj, fNSensObj );
	fReactionFluxes = fInputData->fReactionFluxes;
}

template<class Species>
TFlame<Species>::~TFlame( void )
{
  if (!fSensObjAll)
    {
      // sensitivity stuff
      for ( int i = 0; i < fNSensObj; ++i ) 
	{
	  delete fSensObj[i];
	}
      delete fSensObj;
    }


  DisposeVector( fPressure );
  delete fOutFile;
  delete fInputData;

}

template<class Species>
FILE *TFlame<Species>::GetOutfile( const char *name, FileType type )
{
	char	*fullName = GetOutfileName( name, type );
	FILE	*fp;
	
	//krithika
	if  (type==FileType::kAppend)	{		//open in append mode
		if ( !( fp = fopen( fullName, "a") ) ) { 
			fprintf( stderr, "%s%s\n", "#warning: unable to open file ", fullName  );
			exit(2);
		}
	}
	else {
		if ( !( fp = fopen( fullName, "w") ) ) { 
			fprintf( stderr, "%s%s\n", "#warning: unable to open file ", fullName  );
			exit(2);
		}
	}
	
#ifdef MACOSX
	if ( type == FileType::kData ) {
		char	buffSystem[128];
		sprintf( buffSystem, "/Developer/Tools/SetFile -t \"TEXT\" -c \"QKPT\" %s", fullName );
		system( buffSystem );
	}
#endif

	return fp;
}

template<class Species>
char *TFlame<Species>::GetOutfileName( const char *name, FileType type )
{
	if ( type == FileType::kNone ) {
	  sprintf( fOutFile, "%s%s", fOutputPath, name );
	}
	else {
	  //krithika
	//sprintf( fOutFile, "%s%s.%cout", fOutputPath, name, ( type == FileType::kData ) ? 'd' : 't' );
	sprintf( fOutFile, "%s%s.%cout", fOutputPath, name, ((type == FileType::kData)||(type==FileType::kAppend) ) ? 'd' : 't' );
	}
	
	return fOutFile;
}

template<class Species>
int TFlame<Species>::CheckSolution( Double &temp, Double *Y, int YLen )
{
	int			i;
	
#ifdef HP
		if ( isnan( temp ) ) {
			fprintf( stderr, "%s\n", "#error: temperature = NAN"  );
			return 1;
		}
#endif
	if ( !( temp > 10.0 ) ) {
		temp = 10.0;
	}
	else if ( !( temp < 20000) ) {
		fprintf( stderr, "%s\n", "#error: temperature > 10000"  );
		temp = 10000;
	}
	for ( i = 0; i < YLen; ++i ) {
#ifdef HP
		if ( isnan( Y[i] ) ) {
			fprintf( stderr, "%s%d%s\n", "#error: Y[", i, "] = NAN"  );
			return 1;
		}
#endif
		if ( fClipNegativeConcs ) {
			Y[i] = ( Y[i] <= 0.0 ) ? 1.0e-60 : Y[i];
		}
		Y[i] = MIN( Y[i], 1.0 );
	}
	
	return 0;
}

template<class Species>
Double TFlame<Species>::GetNu( char *globalReaction, const char *name )
{
//gb
	char 	*start;
	char 	*startf;
	char 	ptr[256] = "";
	char	ptr1[256] = "";
	char	*species;
	char	*coeff;
	char	*species_name;
	const char ch = '+';
	int 	len_startf;
// find name in global reaction
// a1F1 + a2F2 + ... ->(==) b1P1 + b2P2 + ...
// find first '+' in globalReac.
	if ( !( startf = strchr( globalReaction, ch ) ) ) {
		return -1.0 ;
	}
// copy left side of globalReaction in start
	startf = strstr( globalReaction, "==" );
	len_startf = startf - globalReaction;
	strncpy( ptr, globalReaction, len_startf );

	if ( !startf ) {
		startf = strstr( globalReaction, ">" );
		--startf;
		len_startf = startf - globalReaction;
		strncpy ( ptr, globalReaction, len_startf );
	}
	strcpy( ptr1, ptr );
	start = ptr1;
// Find species on the left hand side.
	while(*startf != '\0') {
		ptr[0] = '\0';
		startf = strchr(start, ch);
		len_startf = startf - start;

			if ( !startf ) {	// last string in left side of globalReac.
				species = start;
				len_startf = strlen(species);
				coeff = species;
				while ( isdigit(*(coeff)) || *(coeff) == '.' ) ++coeff;
				species_name = coeff;
			if ( strcmp( species_name, name ) == 0 ) {
				len_startf = species_name - species;
				if ( len_startf == 0 ) { 
					return 1.0;
				}
						strncpy( ptr, species, len_startf);
						ptr[len_startf] = '\0';
						return (Double) atof(ptr);
			
					}
			else	
				return -1.0;
			}

		strncpy( ptr, start, len_startf);
		ptr[len_startf] = '\0';
		species = ptr;
		coeff = species;
		while ( isdigit(*(coeff)) || *(coeff) == '.' ) ++coeff;
		species_name = coeff;
		if ( strcmp( species_name, name ) == 0 ) {
			len_startf = species_name - species;
			if ( len_startf == 0 ) { 
				return 1.0;
			}
			ptr1[0] = '\0';
			strcpy ( ptr1, species );
			ptr1[len_startf] = '\0';
			return (Double) atof(ptr1);	
		}
		
		++startf;
		start = startf;
	}

	return -1.0;
}

template<class Species>
Double TFlame<Species>::GetNuProduct( char *globalReaction, const char *name )
{
//gb
	char 	*start;
	char 	*startf;
	char 	ptr[256] = "";
	char	ptr1[256] = "";
	char	*species;
	char	*coeff;
	char	*species_name;
	const char ch = '+';
	int 	len_startf;
// find name in global reaction
// a1F1 + a2F2 + ... ->(==) b1P1 + b2P2 + ...
// find first '+' in globalReac.
	if ( !( startf = strchr( globalReaction, ch ) ) ) {
		return -1.0 ;
	}
// copy right side of globalReaction in start
	startf = strstr( globalReaction, "==" );
	if ( startf ) {
	  startf = startf + 2;
	} else {
	  startf = strstr( globalReaction, ">" );
	  startf = startf + 1;
	}
	len_startf = strlen(globalReaction) - ( startf - globalReaction );
	strncpy( ptr, startf, len_startf );

	strcpy( ptr1, ptr );
	start = ptr1;
// Find species on the left hand side.
	while(*startf != '\0') {
		ptr[0] = '\0';
		startf = strchr(start, ch);
		len_startf = startf - start;

			if ( !startf ) {	// last string in left side of globalReac.
				species = start;
				len_startf = strlen(species);
				coeff = species;
				while ( isdigit(*(coeff)) || *(coeff) == '.' ) ++coeff;
				species_name = coeff;
			if ( strcmp( species_name, name ) == 0 ) {
				len_startf = species_name - species;
				if ( len_startf == 0 ) { 
					return 1.0;
				}
						strncpy( ptr, species, len_startf);
						ptr[len_startf] = '\0';
						return (Double) atof(ptr);
			
					}
			else	
				return -1.0;
			}

		strncpy( ptr, start, len_startf);
		ptr[len_startf] = '\0';
		species = ptr;
		coeff = species;
		while ( isdigit(*(coeff)) || *(coeff) == '.' ) ++coeff;
		species_name = coeff;
		if ( strcmp( species_name, name ) == 0 ) {
			len_startf = species_name - species;
			if ( len_startf == 0 ) { 
				return 1.0;
			}
			ptr1[0] = '\0';
			strcpy ( ptr1, species );
			ptr1[len_startf] = '\0';
			return (Double) atof(ptr1);	
		}
		
		++startf;
		start = startf;
	}

	return -1.0;
}

template<class Species>
void TFlame<Species>::PrintFlameletVector( int len, Double *mat, const char *name, FILE *fp, int inc )
{
	int k;
	Double val;
	fprintf( fp, "%s\n", name );
	for ( k = 0; k < len; ++k ) {
		if ( mat[k*inc] > 0.0 && mat[k*inc] < 1.0e-99 ) {
			val = 1.0e-99;
		}
		else {
			val = mat[k*inc];
		}
		fprintf( fp, "\t%-.6e", val );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}
}

template<class Species>
void T0DFlame<Species>::InitT0DFlame( void )
{ 
	fRadiationName = fInputData->fRadiationName;

   	//if ( !( fReaction = new T0DReaction( fInputData ) ) ) 
	//	FatalError( "memory allocation of T0DReaction failed" );
	//if ( !( fSpecies = new T0DSpecies( fInputData, fReaction ) ) ) 
	//	FatalError( "memory allocation of T0DSpecies failed" );
   	if ( !( fProperties = new T0DProperties( fInputData, &fSpecies ) ) ) 
		FatalError( "memory allocation of T0DProperties failed" );

	if ( fInputData->fWithSoot ) {
		if ( !( fSoot = new T0DSoot( fInputData, &fSpecies ) ) ) {
			FatalError( "memory allocation of TSoot failed" );
		}
	}
	else {
		fSoot = NULL;
	}
}

template<class Species>
T0DFlame<Species>::~T0DFlame( void )
{
	delete fProperties;
	//delete fSpecies;
	//delete fReaction;
}

/**
template<class Species>
void T1DFlame<Species>::ChemEquilSolver( Double temp  )
{	
	
	int 		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	fprintf( stderr, "Temp = %f\n",temp);
	
	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
			fprintf( stderr, "Species %d = %s\n",i, fSpecies.GetNames()[i]);
	}
}*/


template<class Species>
void T0DFlame<Species>::CompLewisNumbers( const char *lewisFile )
{	
	if ( !lewisFile ) {	// Lewis = 1
		int 		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
		Double		*lewis = fSpecies.GetLewisNumber()->vec;
	
		for ( int i = 0; i < nSpeciesInSystem; ++i ) {
			lewis[i] = 1.0;
		}
	}
	else {
		fSpecies.ReadLewisNumbers( lewisFile, fSpecies.GetLewisNumber() );
	}
}

template<class Species>
void T0DFlame<Species>::ComputeProperties( Double temp, Double *Y, Double &pressure
				  , Double &density, EqOfState what, Double * /*sootMoments*/ )
{
//	computes molar mass of mixture, 
//	enthalpy and heat capacity of species
//	heat capacity of mixture and pressure or density, depending on 'EqOfState what'

	int	nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	
	
#if defined (applec) || defined (powerc)
	SpinCursor( 32 );
#endif

	//  First compute molar mass of mixture
	fProperties->ComputeMixtureMolarMass( fProperties->GetMixMolarMassRef(), Y, fSpecies.GetMolarMass()->vec
				, nSpeciesIn );
	
	//  compute properties of Species
	fSpecies.ComputeSpeciesProperties( temp );

	//	compute Delta_i, which is used by CompMixtureProps and Compute_DTherm
#ifndef DELTAINEW
	fSpecies.ComputeDeltaI( fSpecies.GetDeltaI()->vec, fSpecies.GetViscosity()->vec );
#endif

	//  compute properties of the mixture
	fProperties->CompMixtureProps( fSpecies.GetHeatCapacity()->vec
						, fSpecies.GetConductivity()->vec, fSpecies.GetViscosity()->vec, Y, temp
						, pressure, density, what, nSpeciesIn, &fSpecies );
	
	//  compute properties of soot
}

template<class Species>
void T0DFlame<Species>::UpdateThermoProps( Double *Y, Double temp, Double &pressure
				  , Double &density, EqOfState what, Double *sootMoments )
{
	Double	*molarMass = fSpecies.GetMolarMass()->vec;
	Double	*tbConc = ( fReaction.GetTBConc() ) ? fReaction.GetTBConc()->vec : NULL;
	Double	*rateCoeffs = fReaction.GetRateCoefficients()->vec;
	Double	*reactionRate = fReaction.GetReactionRate()->vec;
	Flag	&kNewerThanW = fReaction.GetKNewerThanW();

	ComputeProperties( temp, Y, pressure, density, what, sootMoments );//cai: above subroutine

//cai:undef
#ifdef CHECKPRODRATE
	int ii;
		for ( ii = 0; ii < fSpecies.GetNSpeciesInSystem(); ++ii ) {
		  		Y[ii] = 1.0/fSpecies.GetNSpeciesInSystem();
		}
#endif

//cai:undef
#ifdef PRODRATEFILE
	int		i, j;
	Double	*concs = fReaction.GetMolarConcs()->vec;
	Double	*prodRate = fSpecies.GetProductionRate()->vec;
	fReaction.ComputeConcs( concs, Y, molarMass, density );
	fSpecies.ComputeTheProductionRates( prodRate, reactionRate
						, temp, pressure, concs, rateCoeffs, tbConc );
	for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
		prodRate[i] *= molarMass[i];
	}
	for ( i = fSpecies.GetNSpeciesInSystem(); i < fSpecies.GetNOfSpecies(); ++i ) {
		Y[i] = molarMass[i] / density * concs[i];
	}
#	ifdef CHECKPRODRATE
	for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
		for ( j = 0; j < fReaction.GetNOfReactions(); ++j ) {
			fprintf( stderr, "k[%s] = %.3g\n", fReaction.GetLabels()[j], rateCoeffs[j] );
		}
		for ( j = 0; j < fReaction.GetNOfReactions(); ++j ) {
			fprintf( stderr, "w[%s] = %.3g\n", fReaction.GetLabels()[j], reactionRate[j] );
		}
		for ( i = 0; i < fSpecies.GetNOfSpecies(); ++i ) {
			fprintf( stderr, "C[%s] = %.3g\n", fSpecies.GetNames()[i], concs[i] );
		}
		for ( i = 0; i < fSpecies.GetNOfSpecies(); ++i ) {
			fprintf( stderr, "Y[%s] = %.3g\n", fSpecies.GetNames()[i], Y[i] );
		}
		for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
			fprintf( stderr, "W[%s] = %.3g\n", fSpecies.GetNames()[i], molarMass[i] );
		}
		for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
			fprintf( stderr, "mdot[%s] = %.3g\n", fSpecies.GetNames()[i], fSpecies.GetProductionRate()->vec[i] );
		}
		exit( 2 );
	}
#	endif
#else
	fReaction.CompThirdBodyConcs( tbConc, Y, molarMass, density );
	fReaction.ComputeRateCoefficients( rateCoeffs, fReaction.GetCurrRateCoefficients()
					, fReaction.GetKNewerThanW(), temp, fReaction.GetCurrTemp()
					, pressure, fReaction.GetCurrPressure(), tbConc, &fSpecies );
	fReaction.ComputeReactionRates( reactionRate, kNewerThanW, fReaction.GetCurrReactionRate()
					, rateCoeffs, tbConc, density, Y
					, fReaction.GetCurrConc(), molarMass, &fSpecies );

	fSpecies.ComputeProductionRates( fSpecies.GetProductionRate()->vec, reactionRate );
#	ifdef CHECKPRODRATE
	int i,j;
		for ( j = 0; j < fReaction.GetNOfReactions(); ++j ) {
			fprintf( stderr, "k[%s] = %.3g\n", fReaction.GetLabels()[j], rateCoeffs[j] );
		}
		for ( j = 0; j < fReaction.GetNOfReactions(); ++j ) {
			fprintf( stderr, "w[%s] = %.3g\n", fReaction.GetLabels()[j], reactionRate[j] );
		}
		for ( i = 0; i < fSpecies.GetNOfSpecies(); ++i ) {
			fprintf( stderr, "C[%s] = %.3g\n", fSpecies.GetNames()[i], density*Y[i]/molarMass[i] );
		}
		for ( i = 0; i < fSpecies.GetNOfSpecies(); ++i ) {
			fprintf( stderr, "Y[%s] = %.3g\n", fSpecies.GetNames()[i], Y[i] );
		}
		for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
			fprintf( stderr, "W[%s] = %.3g\n", fSpecies.GetNames()[i], molarMass[i] );
		}
		for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
			fprintf( stderr, "mdot[%s] = %.3g\n", fSpecies.GetNames()[i], fSpecies.GetProductionRate()->vec[i] );
		}
		exit(2);
#	endif
#endif

/**************commented out by Rui to add more radiation models. 10/11/2018***************************/
	// if ( fProperties->GetRadiation() ) {
	// 	fProperties->GetRadiation()->SetRadiation( temp, Y, molarMass, density); 
	// }

	if ( fSoot && sootMoments != NULL ) {
		fSoot->UpdateSoot( &fReaction, &fSpecies, sootMoments, temp, Y, density
							, fProperties->GetMixMolarMass() );
		fSoot->ComputePolymereConcs( Y, temp, density, fSpecies.GetMolarMass()->vec
				, fSoot->GetPij()->mat, fSoot->GetSumPi()->vec, fSoot->GetPAHMoments()->vec
				, sootMoments, &fReaction );
		fSoot->UpdateProductionRates( &fSpecies, &fReaction, fSpecies.GetProductionRate()->vec
					, density, Y, temp, fSoot->GetSumPi()->vec, sootMoments
					, fSoot->GetReactionRate()->vec, fProperties->GetMixMolarMass() );
	}

#undef DEBUGTHERMOPROPS
#ifdef DEBUGTHERMOPROPS
	fprintf( stderr, "temp = %g\tdensity = %g\tM = %g\n", temp, density, fProperties->GetMixMolarMass() );
	for ( int i = 0; i < fReaction.GetNOfReactions(); ++i ) {
		fprintf( stderr, "\treactionRate_%s = %g\n", fReaction.GetLabels()[i]
				, reactionRate[i] );
	}
	for ( i = 0; i < fReaction.GetNOfReactions(); ++i ) {
		if ( reactionRate[i] != 0.0 ) {
			fprintf( stderr, "\treactionRate_%s = %g\n", fReaction.GetLabels()[i]
					, reactionRate[i] );
		}
	}
	
#endif
}

template<class Species>
Double T0DFlame<Species>::ComputeZBilger( Double *Y, Double *YFuelSide, Double *YOxSide )
{
	static const Double	molarMassC = 12.01, 
				molarMassO = 16.0,
				molarMassH = 1.008;
	Double			z, zC, zO, zH, zOF, zCF, zHF, zOO, zCO, zHO;
	TInputDataPtr		inp = GetInputData();
	int			CNum = inp->FindAtomIndex( "C" );
	int			HNum = inp->FindAtomIndex( "H" );
	int			ONum = inp->FindAtomIndex( "O" );
	//SpeciesPtr		species = inp->GetSpecies();		// ### composition can be determined by calling fSpecies.GetComposition
	Double			nuC, nuH, nuO;
	int			oxIndex = inp->FindSpecies( "O2" );
	
	if ( oxIndex < 0 ) {
		return -1.0;
	}

	nuC = ( CNum >= 0 ) ? fSpecies.GetComposition()->mat[ GetFuelIndex() ][ CNum ] : 0;	// ###
	nuH = ( HNum >= 0 ) ? fSpecies.GetComposition()->mat[ GetFuelIndex() ][ HNum ] : 0;	// ###
	nuO = ( ONum >= 0 ) ? fSpecies.GetComposition()->mat[ oxIndex ][ ONum ] : 0;		// ###
	zC = GetElementMassFraction( Y, "C", molarMassC );
	zO = GetElementMassFraction( Y, "O", molarMassO );
	zH = GetElementMassFraction( Y, "H", molarMassH );
	zOF = GetElementMassFraction( YFuelSide, "O", molarMassO );
	zCF = GetElementMassFraction( YFuelSide, "C", molarMassC );
	zHF = GetElementMassFraction( YFuelSide, "H", molarMassH );
	zOO = GetElementMassFraction( YOxSide, "O", molarMassO );
	zCO = GetElementMassFraction( YOxSide, "C", molarMassC );
	zHO = GetElementMassFraction( YOxSide, "H", molarMassH );
	
	if ( nuC == 0 ) {
		z = ( ( zH - zHO ) / ( nuH * molarMassH ) 
				+ 2.0 * ( zOO - zO ) / ( nuO * molarMassO ) )
				/ ( ( zHF - zHO ) / ( nuH * molarMassH ) 
				+ 2.0 * ( zOO - zOF ) / ( nuO * molarMassO ) );
	}
	else if ( nuH == 0 ) {
		z = ( ( zC - zCO ) / ( nuC * molarMassC )
				+ 2.0 * ( zOO - zO ) / ( nuO * molarMassO ) )
				/ ( ( zCF - zCO ) / ( nuC * molarMassC )
				+ 2.0 * ( zOO - zOF ) / ( nuO * molarMassO ) );
	}
	else {
		z = ( ( zC - zCO ) / ( nuC * molarMassC ) + ( zH - zHO ) / ( nuH * molarMassH ) 
				+ 2.0 * ( zOO - zO ) / ( nuO * molarMassO ) )
				/ ( ( zCF - zCO ) / ( nuC * molarMassC ) + ( zHF - zHO ) / ( nuH * molarMassH ) 
				+ 2.0 * ( zOO - zOF ) / ( nuO * molarMassO ) );
	}
	
	return z;
}

template<class Species>
Double T0DFlame<Species>::GetElementMassFraction( Double *Y, const char *const atomName, Double atomMolarMass )
{
	TInputDataPtr		inp = GetInputData();
	int			i, atomNumber = inp->FindAtomIndex( atomName );
	int 			nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	//SpeciesPtr		species = inp->GetSpecies();	// ### composition can be determined by calling fSpecies.GetComposition
	Double			z = 0.0;

	if ( atomNumber > -1 ) {
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			z += fSpecies.GetComposition()->mat[i][atomNumber] * atomMolarMass / fSpecies.GetMolarMass()->vec[i] * Y[i];
		}
	}
	
	return z;
}

template<class Species>
void T1DFlame<Species>::FilldMdYOnePointAnal( TFlameNodePtr flameNode )
{
	fReaction.FilldMdYOnePointAnal( flameNode, fSpecies.GetMolarMass()->vec );
}

template<class Species>
void T1DFlame<Species>::FilldMdTOnePointAnal( TFlameNodePtr flameNode )
{
	fReaction.FilldMdTOnePointAnal( flameNode, GetPressure(), fSpecies.GetMolarMass()->vec );
}

template<class Species>
void T1DFlame<Species>::CompLewisNumbers( const char *lewisFile, const std::string localLewisFileName )
{	
	if ( !lewisFile ) {	// compute Lewis number
		//	Le = lamba / ( rho * cp * D );
		int 		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
		int			maxTempPoint = LocationOfMax( fSolTemp->len, &fSolTemp->vec[kPrev] ) - 1;
		Double		temp = fSolTemp->vec[maxTempPoint];
		Double		*Y = NULL;
		Double		*lewisNumber = fSpecies.GetLewisNumber()->vec;
		Double		*diffusivity = NULL;
		Double		num;
	
#ifdef WRITE_LOCAL_LEWIS
		auto grid = fSolver->bt->GetGrid()->GetCurrentGrid();
		auto nGridPoints = grid->GetNGridPoints();
		auto x = grid->GetX()->vec;
		auto speciesNames = fSpecies.GetNames();
		auto locLewisMat = NewMatrix( nGridPoints, nSpeciesInSystem, kColumnPointers);
		auto locLewis = locLewisMat->mat;
		auto filename = GetOutFileBuff();
		for ( decltype(nGridPoints) k = 0; k < nGridPoints; ++k ) {
			SetFlameNode( k );
			diffusivity = fFlameNode->diffusivity;
			Y = fSolMassFracs->mat[k];
			temp = fSolTemp->vec[k];
			fProperties->ComputeMixtureMolarMass( *fFlameNode->mixMolarMass, Y, 
				fSpecies.GetMolarMass()->vec, nSpeciesInSystem );
			//fSpecies.ComputeSpeciesProperties( fFlameNode, temp );
			fSpecies.ComputeSpeciesProperties( fFlameNode, temp, GetPressure() );

#ifdef OPTIMIZEDELTAI
			fSpecies.ComputeDeltaIOpt( fFlameNode, Y, fFlameNode->GijOverWj, TRUE );
#else
			fSpecies.ComputeDeltaI( fFlameNode->deltaI, Y, fFlameNode->viscosity );
#endif
			//fSpecies.ComputeDeltaI( Y, fFlameNode->viscosity );
			fProperties->CompMixtureProps( fFlameNode, Y, temp, GetPressure(), &fSpecies );
			fSpecies.Compute_D( fFlameNode, temp, Y, GetPressure(), TRUE );
			num = *fFlameNode->mixConductivity / ( *fFlameNode->mixDensity * *fFlameNode->mixHeatCapacity );
			for ( decltype(nSpeciesInSystem) i = 0; i < nSpeciesInSystem; ++i )
				locLewis[i][k] = num / diffusivity[i];
		}
		if(localLewisFileName.empty())
			sprintf( filename, "%s%s", GetOutputPath(), "localLewisNumbers" );
		else
			sprintf( filename, "%s%s", GetOutputPath(), localLewisFileName.c_str() );
		SaveArray( locLewis, nGridPoints, nSpeciesInSystem, kColumnPointers, x, speciesNames, filename );
		DisposeMatrix( locLewisMat );
#endif // WRITE_LOCAL_LEWIS

		SetFlameNode( maxTempPoint );
		diffusivity = fFlameNode->diffusivity;
		Y = fSolMassFracs->mat[maxTempPoint];
		
		fProperties->ComputeMixtureMolarMass( *fFlameNode->mixMolarMass, Y, fSpecies.GetMolarMass()->vec, nSpeciesInSystem );
		fSpecies.ComputeSpeciesProperties( fFlameNode, temp, GetPressure() );
#ifdef OPTIMIZEDELTAI
	fSpecies.ComputeDeltaIOpt( fFlameNode, Y, fFlameNode->GijOverWj, TRUE );
#else
	fSpecies.ComputeDeltaI( fFlameNode->deltaI, Y, fFlameNode->viscosity );
#endif
		fProperties->CompMixtureProps( fFlameNode, Y, temp, GetPressure(), &fSpecies );
		fSpecies.Compute_D( fFlameNode, temp, Y, GetPressure(), TRUE );
		num = *fFlameNode->mixConductivity / ( *fFlameNode->mixDensity * *fFlameNode->mixHeatCapacity );
		for ( int i = 0; i < nSpeciesInSystem; ++i ) {
			lewisNumber[i] = num / diffusivity[i];
		}
	}
	else {
		fSpecies.ReadLewisNumbers( lewisFile, fSpecies.GetLewisNumber() );
	}
}

template<class Species>
void T1DFlame<Species>::ComputeProperties( TFlameNodePtr flameNode, Double temp
									, Double *Y, Double pressure )
{
	int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	Double	*molarMass = fSpecies.GetMolarMass()->vec;
	Flag	newTemp;
	
	//  First compute molar mass of mixture
	fProperties->ComputeMixtureMolarMass( *flameNode->mixMolarMass, Y, molarMass, nSpeciesInSystem );
	
	//  compute properties of Species
	newTemp = fSpecies.ComputeSpeciesProperties( flameNode, temp, pressure );

	//	compute Delta_i, which is used by CompMixtureProps and Compute_DTherm
#ifdef OPTIMIZEDELTAI
	fSpecies.ComputeDeltaIOpt( flameNode, Y, flameNode->GijOverWj, newTemp );
#else
	fSpecies.ComputeDeltaI( flameNode->deltaI, Y, flameNode->viscosity );
#endif

	//  compute properties of the mixture
	fProperties->CompMixtureProps( flameNode, Y, temp, pressure, &fSpecies );

	if ( fSpecies.IsConstantLewisNumber() ) {
		fSpecies.Compute_D( flameNode );
	}
	else {
		fSpecies.Compute_D( flameNode, temp, Y, pressure, newTemp );
		if ( fThermoDiffusion ) {
			fSpecies.Compute_DTherm( flameNode, newTemp );
		}
	}

	//  compute properties of soot
	if ( fSoot ) {
		fSoot->ComputePolymereConcs( Y, temp, flameNode->mixDensity[kCurr]
				, molarMass, flameNode->Pij, flameNode->sumPi, flameNode->pahMoments
				, flameNode->moments, &fReaction );
		fSoot->ComputeDiffusivity( this );
	}
}

template<class Species>
void T1DFlame<Species>::FilldMdYOnePointNum( TFlameNodePtr flameNode )
{
	int				nOfSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	Double			savedVar;
	Double			DeltaY;
	Double			*Y = flameNode->Y[kCurr];
	Double			T = flameNode->temp[kCurr];
	Double			*tBodyConc = flameNode->tBodyConc;
	Double			&tempReaction = (*flameNode->tempReaction);
	Double			&pressureReaction = (*flameNode->pressureReaction);
	Double			*rateCoeff = flameNode->rateCoeff;
	Double			*currRateCoeff = flameNode->currRateCoeff;
	Flag			&kNewerThanW = (*flameNode->kNewerThanW);
	Double			*reactionRate = flameNode->reactionRate;
	Double			*YReaction = flameNode->YReaction;
	Double			*currReacRate = flameNode->currReacRate;
	Double			*productionRate = flameNode->productionRate;
	Double			**dMdY = flameNode->dMdY;
	Double			density;
	Double			*Mmod = fMmod->vec;
	Double			pressure = GetPressure();
	Double			*molarMass = fSpecies.GetMolarMass()->vec;

	// species
	for ( int i = 0; i < nOfSpeciesInSystem; ++i ) {
		savedVar = Y[i];
		DeltaY = ( fabs( Y[i] ) > 1.0e-10) ? .00001 * Y[i] : 1.e-15;
	    Y[i] += DeltaY;
		fProperties->ComputeMixtureMolarMass( *flameNode->mixMolarMass, Y, molarMass, nOfSpeciesInSystem );
		density = pressure * (*flameNode->mixMolarMass) / ( RGAS * T );
		flameNode->mixDensity[kCurr] = density;
#ifdef PRODRATEFILE
		Double	*concs = fReaction.GetMolarConcs()->vec;
		fReaction.ComputeConcs( concs, Y, molarMass, density );
		fSpecies.ComputeTheProductionRates( Mmod, reactionRate
							, T, pressure, concs, rateCoeff, tBodyConc );
		for ( int k = 0; k < fSpecies.GetNSpeciesInSystem(); ++k ) {
			Mmod[k] *= molarMass[k];
		}
#else
		fReaction.CompThirdBodyConcs( tBodyConc, Y, molarMass, density );
		fReaction.ComputeRateCoefficients( rateCoeff, currRateCoeff, kNewerThanW, T
				, tempReaction, pressure, pressureReaction, tBodyConc, &fSpecies );
		fReaction.ComputeReactionRates( reactionRate, kNewerThanW, currReacRate, rateCoeff, tBodyConc, density, Y, YReaction, molarMass, &fSpecies );
		fSpecies.ComputeProductionRates( Mmod, reactionRate );
#endif
		if ( fSoot ) {
			int		nMoments = fSoot->GetNSootMoments();
			Double	**dMdx = flameNode->dMdx;
			Double	*source = flameNode->sootSource;
			Double	*sourceMod = fSoot->GetSourceMod()->vec;

			fSoot->ComputePolymereConcs( Y, T, density
					, molarMass, flameNode->Pij, flameNode->sumPi, flameNode->pahMoments
					, flameNode->moments, &fReaction );
			fSoot->UpdateProductionRates( &fSpecies, &fReaction, Mmod, density, Y, T
						, flameNode->sumPi, flameNode->moments, flameNode->pahReactionRate
						, flameNode->mixMolarMass[kCurr] );
			fSoot->FillSource( sourceMod, this );
	
			for ( int j = 0; j < nMoments; ++j ) {
				dMdx[nMoments+1+i][j] = ( sourceMod[j] - source[j] ) / ( Y[i] - savedVar );
			}
		}
		for ( int j = 0; j < nOfSpeciesInSystem; ++j ) {
			dMdY[i][j] = ( Mmod[j] - productionRate[j] ) / ( Y[i] - savedVar );
		}
		Y[i] = savedVar;
	}

	// clean up
	fProperties->ComputeMixtureMolarMass( *flameNode->mixMolarMass, Y, molarMass, nOfSpeciesInSystem );
	density = pressure * (*flameNode->mixMolarMass) / ( RGAS * T );
	flameNode->mixDensity[kCurr] = density;
	if ( fSoot ) {
		fSoot->ComputePolymereConcs( Y, T, density
				, molarMass, flameNode->Pij, flameNode->sumPi, flameNode->pahMoments
				, flameNode->moments, &fReaction );
	}
}	

template<class Species>
void T1DFlame<Species>::FilldMdTOnePointNum( TFlameNodePtr flameNode )
{
	int				nOfSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	Double			savedVar;
	Double			DeltaT;
	Double			*Y = flameNode->Y[kCurr];
	Double			T = flameNode->temp[kCurr];
	Double			*tBodyConc = flameNode->tBodyConc;
	Double			&tempReaction = (*flameNode->tempReaction);
	Double			&pressureReaction = (*flameNode->pressureReaction);
	Double			*rateCoeff = flameNode->rateCoeff;
	Double			*currRateCoeff = flameNode->currRateCoeff;
	Flag			&kNewerThanW = (*flameNode->kNewerThanW);
	Double			*reactionRate = flameNode->reactionRate;
	Double			*YReaction = flameNode->YReaction;
	Double			*currReacRate = flameNode->currReacRate;
	Double			*productionRate = flameNode->productionRate;
	Double			**dMdY = flameNode->dMdY;
	Double			density;
	Double			*Mmod = fMmod->vec;
	Double			pressure = GetPressure();
	Double			*molarMass = fSpecies.GetMolarMass()->vec;

	// temperature
	savedVar = T;
	DeltaT = ( fabs( T ) > 0.1) ? .00001 * T : 1.e-6;
	T += DeltaT;
	density = pressure * (*flameNode->mixMolarMass) / ( RGAS * T );
	flameNode->mixDensity[kCurr] = density;
#ifdef PRODRATEFILE
	Double	*concs = fReaction.GetMolarConcs()->vec;
	fReaction.ComputeConcs( concs, Y, molarMass, density );
	fSpecies.ComputeTheProductionRates( Mmod, reactionRate
						, T, pressure, concs, rateCoeff, tBodyConc );
	for ( int k = 0; k < fSpecies.GetNSpeciesInSystem(); ++k ) {
		Mmod[k] *= molarMass[k];
	}
#else
	fReaction.CompThirdBodyConcs( tBodyConc, Y, molarMass, density );
	fReaction.ComputeRateCoefficients( rateCoeff, currRateCoeff, kNewerThanW, T
			, tempReaction, pressure, pressureReaction, tBodyConc, &fSpecies );
	fReaction.ComputeReactionRates( reactionRate, kNewerThanW, currReacRate, rateCoeff, tBodyConc, density, Y, YReaction, molarMass, &fSpecies );
	fSpecies.ComputeProductionRates( Mmod, reactionRate );
#endif
	if ( fSoot ) {
		int		nMoments = fSoot->GetNSootMoments();
		Double	**dMdx = flameNode->dMdx;
		Double	*source = flameNode->sootSource;
		Double	*sourceMod = fSoot->GetSourceMod()->vec;
		fSoot->ComputePolymereConcs( Y, T, density
				, molarMass, flameNode->Pij, flameNode->sumPi, flameNode->pahMoments
				, flameNode->moments, &fReaction );
		fSoot->UpdateProductionRates( &fSpecies, &fReaction, Mmod, density, Y, T
				, flameNode->sumPi, flameNode->moments, flameNode->pahReactionRate
				, flameNode->mixMolarMass[kCurr]);
		fSoot->FillSource( sourceMod, this );

		for ( int j = 0; j < nMoments; ++j ) {
			dMdx[nMoments][j] = ( sourceMod[j] - source[j] ) / ( T - savedVar );
		}
	}
	for ( int j = 0; j < nOfSpeciesInSystem; ++j ) {
		dMdY[nOfSpeciesInSystem][j] = ( Mmod[j] - productionRate[j] ) / ( T - savedVar );
	}
	T = savedVar;
		
	// clean up
	density = pressure * (*flameNode->mixMolarMass) / ( RGAS * T );
	flameNode->mixDensity[kCurr] = density;
	if ( fSoot ) {
		fSoot->ComputePolymereConcs( Y, T, density
				, molarMass, flameNode->Pij, flameNode->sumPi, flameNode->pahMoments
				, flameNode->moments, &fReaction );
	}
}	

template<class Species>
void T1DFlame<Species>::FilldMomdMomOnePointNum( TFlameNodePtr flameNode )
{
	int				nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int				nMoments = fSoot->GetNSootMoments();
	int				sootOff = fSoot->GetOffsetSootMoments();
	Double			savedVar;
	Double			DeltaM;
	Double			*Y = flameNode->Y[kCurr];
	NodeInfoPtr		nodeInfo = fSolver->bt->GetNodeInfo();
	Double			*MOverRho = &nodeInfo->y[sootOff];
	Double			T = flameNode->temp[kCurr];
	Double			**dMdx = flameNode->dMdx;
	Double			density = flameNode->mixDensity[kCurr];
	Double			*source = flameNode->sootSource;
	Double			*sourceMod = fSoot->GetSourceMod()->vec;
	Double			pressure = GetPressure();
	Double			*molarMass = fSpecies.GetMolarMass()->vec;

	// species
	for ( int i = 0; i < nMoments; ++i ) {
		savedVar = MOverRho[i];
		DeltaM = ( fabs( MOverRho[i] ) > 1.0e-10 ) ? .00001 * MOverRho[i] : 1.e-15;
	    MOverRho[i] += DeltaM;

		fSoot->MOverRhoToM( MOverRho, flameNode->moments, nMoments
				, Y, T, pressure, molarMass, nSpeciesIn, fProperties );

		fSoot->ComputePolymereConcs( Y, T, density
				, molarMass, flameNode->Pij, flameNode->sumPi, flameNode->pahMoments
				, flameNode->moments, &fReaction );
				

		fSoot->FillSource( sourceMod, this );
		
		for ( int j = 0; j < nMoments; ++j ) {
			dMdx[i][j] = ( sourceMod[j] - source[j] ) / ( MOverRho[i] - savedVar );
		}
		MOverRho[i] = savedVar;

		fSoot->MOverRhoToM( MOverRho, flameNode->moments, nMoments
				, Y, T, pressure, molarMass, nSpeciesIn, fProperties );
	}
}	

template<class Species>
void T1DFlame<Species>::InitT1DFlame( void )
{ 
   	//if ( !( fReaction = new T1DReaction( fInputData ) ) ) 
	//	FatalError( "memory allocation of T1DReaction failed" );
	//if ( !( fSpecies = new T1DSpecies( fInputData, fReaction ) ) ) 
	//	FatalError( "memory allocation of T1DSpecies failed" );
   	if ( !( fProperties = new T1DProperties( fInputData, &fSpecies ) ) ) 
		FatalError( "memory allocation of T1DProperties failed" );

	if ( fInputData->fWithSoot ) {
		if ( !( fSoot = new T1DSoot( fInputData, &fSpecies ) ) ) {
			FatalError( "memory allocation of TSoot failed" );
		}
		fSoot->PrintPAHReactions( this, &fSpecies );
		fSoot->PrintSootReactions( this, &fSpecies );
	}
	else {
		fSoot = NULL;
	}

	if ( fInputData->fStrainRate ) {
		fStrainRate = NewVector( fInputData->fStrainRate->len );
		copy_vec( fStrainRate, fInputData->fStrainRate );
		fStrainRate->len = 0;
		fprintf( stderr, "%s%g\n", "initial strainrate is ", GetStrainRate()  );
	}
	else {
		fStrainRate = NewVector( 1 );
		fStrainRate->len = 0;
		fStrainRate->vec[0] = 0.0;
	}
	
	fGeometry = ( fInputData->fIsAxiSymmetric ) ? 1.0 : 0.0;
	fThermoDiffusion = fInputData->fThermoDiffusion;

	fPrintRHSSpecies = fInputData->fPrintRHSSpecies;
	fPrintRHSTemp = fInputData->fPrintRHSTemp;

	fContinSide = fInputData->fContinSide;

	TBVPInputPtr input = new TBVPInput( fInputData->fNVariables );
	if ( !input ) FatalError( "memory allocation of TBVPInput failed" );
	SetBVPInput( input );
	
	fSolver = new TBVPSolver( input );
	if ( !fSolver ) FatalError( "memory allocation of TBVPSolver failed" );

	delete input;
	
	TNewtonPtr	bt = GetSolver()->bt;
	
	int maxGridPoints = bt->GetMaxGridPoints();

	fSolTemp = NewVector( maxGridPoints + 2 );
	fSolMassFracs = NewMatrix( fSpecies.GetNOfSpecies(), maxGridPoints + 2, kColumnPointers );
	fSolTemp->vec = &fSolTemp->vec[kNext];
	fSolTemp->len -= 2;
	fSolMassFracs->mat = &fSolMassFracs->mat[kNext];
	fSolMassFracs->cols -= 2;

	fSavedTemp = NewVector( maxGridPoints + 2 );
	fSavedMassFracs = NewMatrix( fSpecies.GetNOfSpecies(), maxGridPoints + 2, kColumnPointers );
	fSavedTemp->vec = &fSavedTemp->vec[kNext];
	fSavedTemp->len -= 2;
	fSavedMassFracs->mat = &fSavedMassFracs->mat[kNext];
	fSavedMassFracs->cols -= 2;

	fSavedGrid = NewVector( maxGridPoints + 2 );
	fSavedGrid->len -= 2;
	fSavedGrid->vec = &fSavedGrid->vec[kNext];

	fMmod = NewVector( fSpecies.GetNSpeciesInSystem() + 1 );

	fNoDiffCorr = fInputData->fNoDiffCorr;
	fDiffusivityCorrection = NewVector( maxGridPoints+2 );
	fDiffusivityCorrection->vec = &fDiffusivityCorrection->vec[kNext];

	fFlameNode = new TFlameNode();
	if ( fUseNumericalJac ) {
		fprintf( stderr, "%s\n", "numerical evaluation of Jacobian"  );
		fFlameNodeSaved = NewTFlameNodeSaved();
	}
	
	if ( fUseNumericalDM || fReaction.IsReducedMech() ) {
		fprintf( stderr, "%s\n", "numerical evaluation of jacobian matrix entries of source term"  );
		FilldMdYOnePointPtr = &T1DFlame<Species>::FilldMdYOnePointNum;
		FilldMdTOnePointPtr = &T1DFlame<Species>::FilldMdTOnePointNum;
		if ( fSoot ) {
			FilldMomdMomOnePointPtr = &T1DFlame<Species>::FilldMomdMomOnePointNum;
		}
		else {
			FilldMomdMomOnePointPtr = NULL;
		}
	}
	else {
		FilldMdYOnePointPtr = &T1DFlame<Species>::FilldMdYOnePointAnal;
		FilldMdTOnePointPtr = &T1DFlame<Species>::FilldMdTOnePointAnal;
		FilldMomdMomOnePointPtr = NULL;
	}

    // radiation
    fRadiationName = fInputData->fRadiationName;
    fRadiativeFrac = fInputData->fRadiativeFrac;


}

template<class Species>
void T1DFlame<Species>::FilldMdYOnePoint( TFlameNodePtr flameNode )
{
	(this->*T1DFlame<Species>::FilldMdYOnePointPtr)( flameNode );
}

template<class Species>
void T1DFlame<Species>::FilldMdTOnePoint( TFlameNodePtr flameNode )
{
	(this->*T1DFlame<Species>::FilldMdTOnePointPtr)( flameNode );
}

template<class Species>
void T1DFlame<Species>::FilldMomdMomOnePoint( TFlameNodePtr flameNode )
{
	(this->*T1DFlame<Species>::FilldMomdMomOnePointPtr)( flameNode );
}

template<class Species>
void T1DFlame<Species>::UpdateDimensions( int len )
{
	fSolTemp->len = len;
	fSolMassFracs->cols = len;

	if ( fSoot ) {
		fSoot->UpdateDimensions( len );
	}
}

template<class Species>
void T1DFlame<Species>::UpdateSolution( Double *y, int gridPoint )
{
	int		i;
	int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	int		firstSpeciesOff = GetOffsetFirstSpecies();
	Double	*massFracs = fSolMassFracs->mat[gridPoint];

	//std::cout << "y[GetOffsetTemperature()]: " << y[GetOffsetTemperature()] << '\n';

	fSolTemp->vec[gridPoint] = y[GetOffsetTemperature()];
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		massFracs[i] = y[i+firstSpeciesOff];
	}

	if ( fSoot ) {
		fSoot->UpdateSolution( this, y, gridPoint );
	}
}

template<class Species>
void T1DFlame<Species>::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		i, k, nGridPoints = yMat->cols;
	int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	int		tempOff = GetOffsetTemperature();
	int		firstSpeciesOff = GetOffsetFirstSpecies();
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;
	Double	*temp = fSolTemp->vec;
	Double	**massFracs = fSolMassFracs->mat;

//	set boundary values
	//std::cout << "yLeft[tempOff]: " << yLeft[tempOff] << '\n';
	//std::cout << "yRight[tempOff]: " << yRight[tempOff] << '\n';
	temp[kPrev] = yLeft[tempOff];
	temp[nGridPoints] = yRight[tempOff];
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		massFracs[kPrev][i] = yLeft[i+firstSpeciesOff];
		massFracs[nGridPoints][i] = yRight[i+firstSpeciesOff];
		//std::cout << "massFracs[" << kPrev << "][" << i << "]: " << massFracs[kPrev][i] << '\n';
		//std::cout << "massFracs[" << nGridPoints << "][" << i << "]: " << massFracs[nGridPoints][i] << '\n';
	}

//	set inner values
	for ( k = 0; k < nGridPoints; ++k ) {
		temp[k] = y[k][tempOff];
		//std::cout << "temp[" << k << "]: " << temp[k] << '\n';
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			massFracs[k][i] = y[k][i+firstSpeciesOff];
			//if (k%100 == 0 && i%2==0)
			//std::cout << "massFracs[" << k << "][" << i << "]: " << massFracs[k][i] << '\n';
		}
	}

	if ( fSoot ) {
		fSoot->UpdateSolution( this, yMat, yLeftVec, yRightVec );
	}
}

template<class Species>
void T1DFlame<Species>::SolutionToSolver( void )
{
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int		i, k, nGridPoints = fSolTemp->len;
	int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
	int		tempOff = GetOffsetTemperature();
	int		firstSpeciesOff = GetOffsetFirstSpecies();
	Double	**y = grid->GetY()->mat;
	Double	*temp = fSolTemp->vec;
	Double	**massFracs = fSolMassFracs->mat;

//	set inner values
	for ( k = 0; k < nGridPoints; ++k ) {
		y[k][tempOff] = temp[k];
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			y[k][i+firstSpeciesOff] = massFracs[k][i];
		}
	}
	
	if ( fSoot ) {
		fSoot->SolutionToSolver( y );
	}
}

template<class Species>
void T1DFlame<Species>::SaveSolution( void )
{
	int		i, k;
	int		len = fSolMassFracs->cols;
	int 	nSpecies = fSolMassFracs->rows;
	Double	*temp = fSolTemp->vec;
	Double	**massFracs = fSolMassFracs->mat;
	Double	*savedTemp = fSavedTemp->vec;
	Double	**savedMassFracs = fSavedMassFracs->mat;

	SaveGrid();

	fSavedTemp->len = fSolTemp->len;
	fSavedMassFracs->cols = fSolMassFracs->cols;
	fSavedMassFracs->rows = fSolMassFracs->rows;

	for ( k = -1; k <= len; ++k ) {
		savedTemp[k] = temp[k];
		for ( i = 0; i < nSpecies; ++i ) {
			savedMassFracs[k][i] = massFracs[k][i];
		}
	}
	
	if ( fSoot ) {
		fSoot->SaveSolution();
	}
}

template<class Species>
void T1DFlame<Species>::RestoreSolution( void )
{
	int		i, k;
	int		len = fSavedMassFracs->cols;
	int 	nSpecies = fSavedMassFracs->rows;
	Double	*temp = fSolTemp->vec;
	Double	**massFracs = fSolMassFracs->mat;
	Double	*savedTemp = fSavedTemp->vec;
	Double	**savedMassFracs = fSavedMassFracs->mat;

	RestoreGrid();

	for ( k = -1; k <= len; ++k ) {
		temp[k] = savedTemp[k];
		for ( i = 0; i < nSpecies; ++i ) {
			massFracs[k][i] = savedMassFracs[k][i];
		}
	}
	
	if ( fSoot ) {
		fSoot->RestoreSolution();
	}
}

template<class Species>
void T1DFlame<Species>::SaveGrid( void )
{
	int			k;
	Double		*savedGrid = fSavedGrid->vec;
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int			len = grid->GetNGridPoints();
	Double		*x = grid->GetX()->vec;

	fSavedGrid->len = len;
	savedGrid[kPrev] = bt->GetLeft();
	
	for ( k = 0; k < len; ++k ) {
		savedGrid[k] = x[k];
	}
	
	savedGrid[len] = bt->GetRight();
}

template<class Species>
void T1DFlame<Species>::RestoreGrid( void )
{
	int			k;
	int			len = fSavedGrid->len;
	Double		*savedGrid = fSavedGrid->vec;
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	Double		*x = grid->GetX()->vec;

	bt->SetLeft( savedGrid[kPrev] );
	
	for ( k = 0; k < len; ++k ) {
		x[k] = savedGrid[k];
	}
	
	bt->SetRight( savedGrid[len] );

	grid->AdjustNGridPoints( len );
	fSolver->UpdateAllDimensions( len );	
}

template<class Species>
T1DFlame<Species>::~T1DFlame( void )
{
	fDiffusivityCorrection->vec = &fDiffusivityCorrection->vec[kPrev];
	DisposeVector( fDiffusivityCorrection );
	if ( fUseNumericalJac ) {
		DisposeTFlameNodeSaved( fFlameNodeSaved );
	}
	delete fFlameNode;

	DisposeVector( fMmod );

	fSavedGrid->vec = &fSavedGrid->vec[kPrev];
	DisposeVector( fSavedGrid );
	
	fSavedMassFracs->mat = &fSavedMassFracs->mat[kPrev];
	fSavedTemp->vec = &fSavedTemp->vec[kPrev];
	DisposeMatrix( fSavedMassFracs );
	DisposeVector( fSavedTemp );

	fSolMassFracs->mat = &fSolMassFracs->mat[kPrev];
	fSolTemp->vec = &fSolTemp->vec[kPrev];
	DisposeMatrix( fSolMassFracs );
	DisposeVector( fSolTemp );

	delete fSolver;
	if ( fStrainRate ) {
		DisposeVector( fStrainRate );
	}
	delete fProperties;
	// delete fSpecies;
	// delete fReaction;
}
template<class Species>
void T1DFlame<Species>::UpdateThermoProps( vector<double> physicalGrid, double Delta, double chiRef, vector <double> HRR, string radiationName)
{
	int					k;
	int					nSpeciesIn = fSpecies.GetNSpeciesInSystem();
	int					tempOffset = GetOffsetTemperature();
	int					speciesOffset = GetOffsetFirstSpecies();
	TNewtonPtr			bt = fSolver->bt;
	int					currentGridPoints = bt->GetCurrentGridPoints();
	Double				density;
	Double				pressure = GetPressure();
	Double				*molarMass = fSpecies.GetMolarMass()->vec;
	Double				**rateCoeffs = fReaction.GetRateCoefficients()->mat;
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
	T1DRadiationPtr		radiation = fProperties->GetRadiation();
	Double 				**reactionRate = fReaction.GetReactionRate()->mat;
	Double				**tBConc = fReaction.GetTBConcentrations()->mat;
	Double				*temp = fSolTemp->vec;
	Double				**Y = fSolMassFracs->mat;
    Double			**massFracs = fSolMassFracs->mat;
    Double			*mixMolarMass = fProperties->GetMolarMass()->vec;
   

    // index of H2O and CO2 in data vectors
	int fH2OIndex = fInputData->fH2OIndex;
	int fCO2Index = fInputData->fCO2Index;
    // compute molar fraction
      
    vector <double> locMoleMassH2O (physicalGrid.size(), 0.);
    vector <double> locMoleMassCO2 (physicalGrid.size(), 0.);

    // left bc`
    
    locMoleMassH2O[0] = massFracs[-1][fH2OIndex] * mixMolarMass[-1] / molarMass[fH2OIndex];
    locMoleMassCO2[0] = massFracs[-1][fCO2Index] * mixMolarMass[-1] / molarMass[fCO2Index];

    // right bc
    
    locMoleMassH2O[currentGridPoints+1] = massFracs[currentGridPoints][fH2OIndex] * mixMolarMass[currentGridPoints] / molarMass[fH2OIndex];
    locMoleMassCO2[currentGridPoints+1] = massFracs[currentGridPoints][fCO2Index] * mixMolarMass[currentGridPoints] / molarMass[fCO2Index];

    // grid
    
    for ( int i = 0; i < currentGridPoints; ++i){ 
        locMoleMassH2O[i+1] = massFracs[i][fH2OIndex] * mixMolarMass[i] / molarMass[fH2OIndex];
        locMoleMassCO2[i+1] = massFracs[i][fCO2Index] * mixMolarMass[i] / molarMass[fCO2Index];
    }

    fProperties->GetRadiation()->ComputeRadiationRadiation(temp, locMoleMassH2O, locMoleMassCO2, nSpeciesIn, physicalGrid, Delta, chiRef, HRR, radiationName);

}


    template<class Species>
void T1DFlame<Species>::UpdateThermoPropsRadiation()
{
    int					k;
    int					nSpeciesIn = fSpecies.GetNSpeciesInSystem();
    int					tempOffset = GetOffsetTemperature();
    int					speciesOffset = GetOffsetFirstSpecies();
    TNewtonPtr			bt = fSolver->bt;
    int					currentGridPoints = bt->GetCurrentGridPoints();
    Double				density;
    Double				pressure = GetPressure();
    Double				*molarMass = fSpecies.GetMolarMass()->vec;
    Double				**rateCoeffs = fReaction.GetRateCoefficients()->mat;
    NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    T1DRadiationPtr		radiation = fProperties->GetRadiation();
    Double 				**reactionRate = fReaction.GetReactionRate()->mat;
    Double				**tBConc = fReaction.GetTBConcentrations()->mat;
    Double				*temp = fSolTemp->vec;
    Double				**Y = fSolMassFracs->mat;
#define BOGYCHECK
#ifdef BOGYCHECK
    Double				*ant = New1DArray( currentGridPoints + 2 );
    ant = &ant[kNext];
#endif

    for ( k = -1; k <= currentGridPoints; ++k ) {
        CheckSolution( temp[k], Y[k], nSpeciesIn );
    }

    SetFlameNode( kPrev );
    ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], pressure );


#ifdef BOGYCHECK
    if ( fSoot ) ant[-1] = fSoot->GetGamma();
#endif

    for ( k = 0; k < currentGridPoints; ++k ) {
    	
#if defined (applec) || defined (powerc)
        RotateCursor( 32 );
#endif
        bt->SetNodeInfo( this, k );
        SetFlameNode( k );
        ComputeProperties( fFlameNode, temp[k], Y[k], pressure );
#ifdef BOGYCHECK
        if ( fSoot ) ant[k] = fSoot->GetGamma();
#endif
        density = GetProperties()->GetDensity()->vec[k];
#ifdef PRODRATEFILE
        int		i;
        Double	*concs = fReaction.GetMolarConcs()->vec;
        Double	*prodRate = fFlameNode->productionRate;
        fReaction.ComputeConcs( concs, Y[k], molarMass, density );
        fSpecies.ComputeTheProductionRates( prodRate, reactionRate[k]
                , temp[k], pressure, concs, rateCoeffs[k], tBConc[k] );
        for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
            prodRate[i] *= molarMass[i];
        }
        for ( i = fSpecies.GetNSpeciesInSystem(); i < fSpecies.GetNOfSpecies(); ++i ) {
            Y[k][i] = molarMass[i] / density * concs[i];
        }

#else
        fReaction.CompThirdBodyConcs( tBConc[k], Y[k], molarMass, density );
        fReaction.ComputeRateCoefficients( rateCoeffs[k]
                , fReaction.GetCurrRateCoeff()->mat[k], fReaction.GetKNewerThanW()[k]
                , temp[k], fReaction.GetTempReaction()->vec[k], pressure
                , fReaction.GetPressureReaction()->vec[k], tBConc[k], &fSpecies );
        fReaction.ComputeReactionRates( reactionRate[k], fReaction.GetKNewerThanW()[k], fReaction.GetCurrReacRate()->mat[k]
                , rateCoeffs[k], tBConc[k], density
                , Y[k], fReaction.GetYReaction()->mat[k], molarMass, &fSpecies );
        fSpecies.ComputeProductionRates( fFlameNode->productionRate, reactionRate[k] );
#endif
#undef CHECKPRODRATE
#ifdef CHECKPRODRATE
        int j;
        for ( j = 0; j < fReaction.GetNOfReactions(); ++j ) {
            fprintf( stdout, "k[%d] = %.6g\n", j, rateCoeffs[k][j] );
        }
        for (  j = 0; j < fReaction.GetNOfReactions(); ++j ) {
            fprintf( stdout, "w[%d] = %.6g\n", j, reactionRate[k][j] );
        }
        for ( j = 0; j < fSpecies.GetNSpeciesInSystem(); ++j ) {
            fprintf( stdout, "mdot[%s] = %.6g\n", fSpecies.GetNames()[j], fFlameNode->productionRate[j] );
        }
        exit( 2 );
#endif
#undef CHECKPRODRATE
        ComputeDiffusivityCorrection( &Y[k], nodeInfo );
        Double HeatReleaseRate = 0.0;
        for (int j = 0; j < fSpecies.GetNSpeciesInSystem(); ++j ) {
            //Double *ent = fSpecies.GetEnthalpy()->vec;
                                                                  
            HeatReleaseRate += fFlameNode->productionRate[j]*fFlameNode->enthalpy[j];
        }

        HeatReleaseRate = -HeatReleaseRate;

        if ( radiation ) {

            if ( fRadiationName == string("Adiabatic")){

            }
            else if (fRadiationName == string("Thin")){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, temp[k]
                        , Y[k], molarMass, density);
            }
            else if (fRadiationName == "RadiativeFrac"){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, temp[k]
                        , Y[k], molarMass, density, HeatReleaseRate, fRadiativeFrac, fFlameNode->kappa);
            }
            else if (fRadiationName == "WSGG"){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[k]);
            }
            else if (fRadiationName == "WSGGJohansson"){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[k]);
            }
            else if (fRadiationName == "WSGGBordbar"){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[k]);
            }
            else if (fRadiationName == "Grey"){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[k]);
            }
            else if (fRadiationName == "SNB"){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[k]);
            }
            else{
                cout << "#1 Error Radiation Model " <<  fRadiationName << "not found " << endl;
            }
        }
        if ( fSoot ) {
            fSoot->UpdateProductionRates( &fSpecies, &fReaction, fFlameNode->productionRate, density
                    , Y[k], temp[k], fFlameNode->sumPi, fFlameNode->moments, fFlameNode->pahReactionRate
                    , fFlameNode->mixMolarMass[kCurr]);
            fSoot->FillSource( fFlameNode->sootSource, this );
        }

    }

#ifdef CHECKPRODRATE
    exit(2);
#endif
    SetFlameNode( currentGridPoints );
    ComputeProperties( fFlameNode, temp[currentGridPoints], Y[currentGridPoints], pressure );
#ifdef BOGYCHECK
    if ( fSoot ) ant[currentGridPoints] = fSoot->GetGamma();
#endif


    if ( fSoot ) {
#ifdef BOGYCHECK
        Double		dummy;
        int counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
        FILE	*fp = NULL;
        int		k;
        int			nGridPoints = bt->GetCurrentGridPoints();
        Double		left = bt->GetLeft();
        Double		right = bt->GetRight();
        TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
        Double		*x = grid->GetX()->vec;
        sprintf( bt->GetOutFileBuff(), "%sGamma_%d.dout", bt->GetOutputPath(), counter );
        if ( !( fp = fopen( bt->GetOutFileBuff(), "w") ) ) {
            cerr << "#warning: unable to open file" << bt->GetOutFileBuff() << NEWL;
            return;
        }
        fprintf( fp, "*\n" );
        fprintf( fp, "%-12s\t%-12s", "y", "Gamma" );

        fprintf( fp, "\n%-9E", left );
        fprintf( fp, "\t%-9E", ant[-1] );
        for ( k = 0; k < nGridPoints; ++k ) {
            fprintf( fp, "\n%-9E", x[k] );
            fprintf( fp, "\t%-9E", ant[k] );
        }
        fprintf( fp, "\n%-9E", right );
        fprintf( fp, "\t%-9E", ant[nGridPoints] );

        fclose( fp );
#endif
        fSoot->PrintPAHOneStage( 0, bt );
        fSoot->PrintPAHMoments( bt );
        if ( fPrintRHSSpecies ) {
            fSoot->PrintRHSSoot( bt, this );
        }
    }
#ifdef BOGYCHECK
    ant = &ant[kPrev];
    Free1DArray( ant );
#endif
}

    template<class Species>
void T1DFlame<Species>::UpdateThermoProps()
{
    int					k;
    int					nSpeciesIn = fSpecies.GetNSpeciesInSystem();
    int					tempOffset = GetOffsetTemperature();
    int					speciesOffset = GetOffsetFirstSpecies();
    TNewtonPtr			bt = fSolver->bt;
    int					currentGridPoints = bt->GetCurrentGridPoints();
    Double				density;
    Double				pressure = GetPressure();
    Double				*molarMass = fSpecies.GetMolarMass()->vec;
    Double				**rateCoeffs = fReaction.GetRateCoefficients()->mat;
    NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    T1DRadiationPtr		radiation = fProperties->GetRadiation();
    Double 				**reactionRate = fReaction.GetReactionRate()->mat;
    Double				**tBConc = fReaction.GetTBConcentrations()->mat;
    Double				*temp = fSolTemp->vec;
    Double				**Y = fSolMassFracs->mat;
#define BOGYCHECK
#ifdef BOGYCHECK
    Double				*ant = New1DArray( currentGridPoints + 2 );
    ant = &ant[kNext];
#endif

    for ( k = -1; k <= currentGridPoints; ++k ) {
        CheckSolution( temp[k], Y[k], nSpeciesIn );
    }

    SetFlameNode( kPrev );
    ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], pressure );


#ifdef BOGYCHECK
    if ( fSoot ) ant[-1] = fSoot->GetGamma();
#endif
    for ( k = 0; k < currentGridPoints; ++k ) {
#if defined (applec) || defined (powerc)
        RotateCursor( 32 );
#endif
        bt->SetNodeInfo( this, k );
        SetFlameNode( k );
        ComputeProperties( fFlameNode, temp[k], Y[k], pressure );
#ifdef BOGYCHECK
        if ( fSoot ) ant[k] = fSoot->GetGamma();
#endif
        density = GetProperties()->GetDensity()->vec[k];
#ifdef PRODRATEFILE
        int		i;
        Double	*concs = fReaction.GetMolarConcs()->vec;
        Double	*prodRate = fFlameNode->productionRate;
        fReaction.ComputeConcs( concs, Y[k], molarMass, density );
        fSpecies.ComputeTheProductionRates( prodRate, reactionRate[k]
                , temp[k], pressure, concs, rateCoeffs[k], tBConc[k] );
        for ( i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
            prodRate[i] *= molarMass[i];
        }
        for ( i = fSpecies.GetNSpeciesInSystem(); i < fSpecies.GetNOfSpecies(); ++i ) {
            Y[k][i] = molarMass[i] / density * concs[i];
        }

#else
        fReaction.CompThirdBodyConcs( tBConc[k], Y[k], molarMass, density );
        fReaction.ComputeRateCoefficients( rateCoeffs[k]
                , fReaction.GetCurrRateCoeff()->mat[k], fReaction.GetKNewerThanW()[k]
                , temp[k], fReaction.GetTempReaction()->vec[k], pressure
                , fReaction.GetPressureReaction()->vec[k], tBConc[k], &fSpecies );
        fReaction.ComputeReactionRates( reactionRate[k], fReaction.GetKNewerThanW()[k], fReaction.GetCurrReacRate()->mat[k]
                , rateCoeffs[k], tBConc[k], density
                , Y[k], fReaction.GetYReaction()->mat[k], molarMass, &fSpecies );
        fSpecies.ComputeProductionRates( fFlameNode->productionRate, reactionRate[k] );
#endif
#undef CHECKPRODRATE
#ifdef CHECKPRODRATE
        int j;
        for ( j = 0; j < fReaction.GetNOfReactions(); ++j ) {
            fprintf( stdout, "k[%d] = %.6g\n", j, rateCoeffs[k][j] );
        }
        for (  j = 0; j < fReaction.GetNOfReactions(); ++j ) {
            fprintf( stdout, "w[%d] = %.6g\n", j, reactionRate[k][j] );
        }
        for ( j = 0; j < fSpecies.GetNSpeciesInSystem(); ++j ) {
            fprintf( stdout, "mdot[%s] = %.6g\n", fSpecies.GetNames()[j], fFlameNode->productionRate[j] );
        }
        exit( 2 );
#endif
#undef CHECKPRODRATE
        ComputeDiffusivityCorrection( &Y[k], nodeInfo );
        Double HeatReleaseRate = 0.0;
        for (int j = 0; j < fSpecies.GetNSpeciesInSystem(); ++j ) {
            //Double *ent = fSpecies.GetEnthalpy()->vec;

            HeatReleaseRate += fFlameNode->productionRate[j]*fFlameNode->enthalpy[j];
        }

        HeatReleaseRate = -HeatReleaseRate;

        if ( radiation ) {

            if ( fRadiationName == string("Adiabatic")){
            }
            else if (fRadiationName == string("Thin")){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, temp[k]
                        , Y[k], molarMass, density);
            }
            else if (fRadiationName == "RadiativeFrac"){
                fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, temp[k]
                        , Y[k], molarMass, density, HeatReleaseRate, fRadiativeFrac, fFlameNode->kappa);
            }
            else if (fRadiationName == "WSGG"){
            }
            else if (fRadiationName == "WSGGJohansson"){
            }
            else if (fRadiationName == "WSGGBordbar"){
            }
            else if (fRadiationName == "Grey"){
            }
            else if (fRadiationName == "SNB"){
            }
            else{
                cout << "#1 Error Radiation Model " <<  fRadiationName << "not found " << endl;
            }
        }
        if ( fSoot ) {
            fSoot->UpdateProductionRates( &fSpecies, &fReaction, fFlameNode->productionRate, density
                    , Y[k], temp[k], fFlameNode->sumPi, fFlameNode->moments, fFlameNode->pahReactionRate
                    , fFlameNode->mixMolarMass[kCurr]);
            fSoot->FillSource( fFlameNode->sootSource, this );
        }

    }

#ifdef CHECKPRODRATE
    exit(2);
#endif
    SetFlameNode( currentGridPoints );
    ComputeProperties( fFlameNode, temp[currentGridPoints], Y[currentGridPoints], pressure );
#ifdef BOGYCHECK
    if ( fSoot ) ant[currentGridPoints] = fSoot->GetGamma();
#endif


    if ( fSoot ) {
#ifdef BOGYCHECK
        Double		dummy;
        int counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
        FILE	*fp = NULL;
        int		k;
        int			nGridPoints = bt->GetCurrentGridPoints();
        Double		left = bt->GetLeft();
        Double		right = bt->GetRight();
        TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
        Double		*x = grid->GetX()->vec;
        sprintf( bt->GetOutFileBuff(), "%sGamma_%d.dout", bt->GetOutputPath(), counter );
        if ( !( fp = fopen( bt->GetOutFileBuff(), "w") ) ) {
            cerr << "#warning: unable to open file" << bt->GetOutFileBuff() << NEWL;
            return;
        }
        fprintf( fp, "*\n" );
        fprintf( fp, "%-12s\t%-12s", "y", "Gamma" );

        fprintf( fp, "\n%-9E", left );
        fprintf( fp, "\t%-9E", ant[-1] );
        for ( k = 0; k < nGridPoints; ++k ) {
            fprintf( fp, "\n%-9E", x[k] );
            fprintf( fp, "\t%-9E", ant[k] );
        }
        fprintf( fp, "\n%-9E", right );
        fprintf( fp, "\t%-9E", ant[nGridPoints] );

        fclose( fp );
#endif
        fSoot->PrintPAHOneStage( 0, bt );
        fSoot->PrintPAHMoments( bt );
        if ( fPrintRHSSpecies ) {
            fSoot->PrintRHSSoot( bt, this );
        }
    }
#ifdef BOGYCHECK
    ant = &ant[kPrev];
    Free1DArray( ant );
#endif
}

    template<class Species>
void T1DFlame<Species>::ComputeDiffusivityCorrection( Double **Y, NodeInfoPtr nodeInfo )
{
    int		i;
    int 	nOfSpecies = fSpecies.GetNOfSpecies();
    int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    Double	h = nodeInfo->h;
    Double	hm = nodeInfo->hm;
    Double	*diffusivity = fFlameNode->diffusivity;
    Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
    Double	*diffusivityNext = fFlameNode->diffusivityNext;
    Double	*y = Y[kCurr];
    Double	*yPrev = Y[kPrev];
    Double	*yNext = Y[kNext];
#ifdef MOLARDIFFUSION
    Double	M = fFlameNode->mixMolarMass[kCurr];
    Double	MPrev = fFlameNode->mixMolarMass[kPrev];
    Double	MNext = fFlameNode->mixMolarMass[kNext];
    Double	dMdyOverM = 0.0;
#endif
    Double	*diffCorr = fFlameNode->diffCorr;
    Double	locDiffCorr = 0.0;
    if ( nodeInfo->firstPoint ) {
#ifdef MOLARDIFFUSION
        dMdyOverM = FirstDerivUpwind( M, MPrev, hm ) / MPrev;
#endif
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            if (hm == 0.0) {
                //			fprintf(stderr, "## hm = 0:x0: %g\txL %g\n"
                //					, fSolver->bt->GetGrid()->GetCurrentGrid()->GetX()->vec[0]
                //					, fSolver->bt->GetLeft());
            }
            locDiffCorr += diffusivityPrev[i] * FirstDerivUpwind( y[i], yPrev[i], hm );
#ifdef MOLARDIFFUSION
            locDiffCorr += diffusivityPrev[i] * yPrev[i] * dMdyOverM;
#endif
        }
        diffCorr[kPrev] = locDiffCorr;
    }
    else if ( nodeInfo->lastPoint ) {
#ifdef MOLARDIFFUSION
        dMdyOverM = FirstDerivUpwind( MNext, M, h ) / MNext;
#endif
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            locDiffCorr += diffusivityNext[i] * FirstDerivUpwind( yNext[i], y[i], h );
#ifdef MOLARDIFFUSION
            locDiffCorr += diffusivityNext[i] * yNext[i] * dMdyOverM;
#endif
        }
        diffCorr[kNext] = locDiffCorr;
    }

    locDiffCorr = 0.0;
#ifdef MOLARDIFFUSION
    dMdyOverM = FirstDerivUpwind( M, MPrev, hm ) / M;//FirstDeriv( MPrev, M, MNext, hm, h ) / M;
#endif
    for ( i = 0; i < nSpeciesInSystem; ++i ) {
        locDiffCorr += diffusivity[i] * FirstDeriv( yPrev[i], y[i], yNext[i], hm, h );
#ifdef MOLARDIFFUSION
        locDiffCorr += diffusivity[i] * y[i] * dMdyOverM;
#endif
    }
    diffCorr[kCurr] = locDiffCorr;
}

    template<class Species>
void T1DFlame<Species>::WriteRoggFiles( TNewtonPtr /*bt*/ )
{
    fSpecies.WriteRoggsSymbolsData();
    fReaction.WriteRoggsMechanismData( &fSpecies, fInputData );
}

    template<class Species>
void T1DFlame<Species>::PostConvergence( void *object )
{
    TNewtonPtr 			bt = GetSolver()->bt;
    int					isConverged = bt->GetConvergeNewton();
    TAdaptiveGridPtr	grid = bt->GetGrid();
    NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int					firstSpecies = GetOffsetFirstSpecies();
    int					tempOffset = GetOffsetTemperature();
    char				tail[32];
    int					toSpeciesIndex = GetToSpeciesIndex() + firstSpecies;
    VectorPtr			aVec = GetStrainRateVector();
    int 				contVar;
    Double				contSmall;
    static Double		aConv = -1.0;	//means not set
    static Double		aNotConv = -1.0;	//means not set
    static Double		pConv = -1.0;	//means not set
    static Double		pNotConv = -1.0;	//means not set
    static Double		contConv = -1.0;	//	means not set
    static Double		contNotConv = -1.0;	//	means not set

    if ( isConverged ) {
        if ( toSpeciesIndex-firstSpecies >= 0 ) {
            int 	fromSpeciesIndex = GetFromSpeciesIndex() + firstSpecies;
            char	**names = GetSpecies()->GetNames();
            Double	*yBound;
            Double	*bcBound;

            if ( fromSpeciesIndex-firstSpecies < 0 ) {
                fromSpeciesIndex = GetFuelIndex() + firstSpecies;
            }

            if ( fContInc == 0.0 || fContinSide == kNoSide ) {
                fprintf( stderr, "%s%s%s\n", "#warning: increment of continuation is zero", NEWL, 
                        "          or no continuation side specified"  );
                bt->SetLeaveContin();
                return;
            }

            if ( fContinSide == kLeft ) {
                yBound = grid->GetFine()->GetYLeft()->vec;
                bcBound = grid->GetFine()->GetBcLeft()->vec;
            }
            else {
                yBound = grid->GetFine()->GetYRight()->vec;
                bcBound = grid->GetFine()->GetBcRight()->vec;
            }

            // write Output
            if ( bcBound[toSpeciesIndex] < 0.01 ) {
                sprintf( tail, "%s_%g", names[toSpeciesIndex-firstSpecies], bcBound[toSpeciesIndex] );
            }
            else {
                sprintf( tail, "%s_%d", names[toSpeciesIndex-firstSpecies], (int) ( bcBound[toSpeciesIndex] * 100 ) );
            }
            bt->WriteOutput( object, NULL, tail );

            if ( yBound[toSpeciesIndex] * fContInc < fContBound * fContInc ) {
                // muss geaendert werden
                // 1.		0 <= Y <= 1
                // 2.		deltaYTo = deltaYFrom
                fSolver->ReInit();
                Double	incTo = MIN( 1.0, MAX( 1.0e-60, yBound[toSpeciesIndex] + fContInc ) ) 
                    - yBound[toSpeciesIndex];
                incTo = MIN( fContBound - yBound[toSpeciesIndex], incTo );
                Double	incFrom = -MIN( 1.0, MAX( 1.0e-60, yBound[fromSpeciesIndex] - fContInc ) ) 
                    + yBound[fromSpeciesIndex];
                Double	inc = MIN( fabs( incFrom ), fabs( incTo ) ) * ( ( fContInc >= 0 ) ? 1.0 : -1.0 );
                if ( inc < fContInc - 1.0e-10 ) {
                    fContBound = yBound[toSpeciesIndex];
                }
                fprintf( stderr, "incTo = %g\tincFrom = %g\tinc = %g\n", incTo, incFrom, inc );
                fprintf( stderr, "yBound[toSpeciesIndex] = %g\n", yBound[toSpeciesIndex] );
                fprintf( stderr, "yBound[fromSpeciesIndex] = %g\n", yBound[fromSpeciesIndex] );
                bcBound[toSpeciesIndex] = yBound[toSpeciesIndex] 
                    = yBound[toSpeciesIndex] + inc;
                bcBound[fromSpeciesIndex] = yBound[fromSpeciesIndex] 
                    = yBound[fromSpeciesIndex] - inc;
                fprintf( stderr, "%s%s%s%g%s%s%s%g\n", "Y_", names[toSpeciesIndex-firstSpecies]
                        , " is now ", bcBound[toSpeciesIndex], " and Y_"
                        , names[fromSpeciesIndex-firstSpecies], " is now ", bcBound[fromSpeciesIndex] );
            }
            else {
                bt->SetLeaveContin();
            }
        }
        else if ( bt->GetNEquations() < bt->GetNVariables() ) {
            ConstStringArray	varNames = GetVariableNames();
            int	inc = ( 1 > (int)fContInc ) ? 1 : (int)fContInc;
            sprintf( tail, "eq%d", bt->GetNEquations() );
            bt->WriteOutput( object, NULL, tail );
            fSolver->ReInit();

            bt->SetNOfEquations( bt->GetNEquations() + inc );//( fVariablesWithoutSpecies + nOfSpecies - 1 );
            fprintf( stderr, "%s%s%s%d%s\n", "solving the system up to "
                    , varNames[bt->GetNEquations()-1]
                    , " yields ", bt->GetNEquations(), " equations"  );
        }
        else {
            bt->WriteOutput( object, NULL, "" );
            if ( aVec && ( aVec->len < aVec->phys_len - 1 || aNotConv > 0.0 ) ) {
                fSolver->ReInit();
                aConv = GetStrainRate();
                if ( aNotConv < 0.0 ) {
                    ++aVec->len;
                }
                else {
                    aVec->vec[aVec->len] = aNotConv;
                    aNotConv = -1.0;
                }
                fprintf( stderr, "%s%g\n", "strainRate is now ", GetStrainRate()  );
                return;
            }
            if ( fPressure && ( fPressure->len < fPressure->phys_len - 1 || pNotConv > 0.0 ) ) {
                fSolver->ReInit();
                pConv = GetPressure();
                if ( pNotConv < 0.0 ) {
                    ++fPressure->len;
                }
                else {
                    fPressure->vec[fPressure->len] = pNotConv;
                    pNotConv = -1.0;
                }
                fprintf( stderr, "%s%g bar\n", "pressure is now ", GetPressure()/1.0e5  );
                return;
            }
            if ( fContinType == kTemperature || fContinType == kVelocity 
                    || fContinType == kPolytrope || fContinType == kPressure 
                    || fContinType == kMomentum || fContinType == kSeshaTemp ) {
                if ( fContInc == 0.0 || ( fContinSide == kNoSide && fContinType != kPressure ) ) {
                    fprintf( stderr, "%s%s%s\n", "#warning: increment of continuation is zero", NEWL, 
                            "          or no continuation side specified"  );
                    bt->SetLeaveContin();
                    return;
                }

                switch( fContinType ) {
                    case kTemperature:
                        contVar = GetOffsetTemperature();
                        contSmall = 1.0;
                        break;
                    case kVelocity:
                        contVar = GetOffsetVVelocity();
                        contSmall = 0.001;
                        break;
                    case kMomentum:
                        contVar = GetOffsetVVelocity();
                        contSmall = 0.01;
                        fContinSide = kBothSides;
                        break;
                    case kSeshaTemp:
                        contVar = GetOffsetTemperature();
                        contSmall = 1.0;
                        break;
                    case kPolytrope:
                        contVar = GetOffsetTemperature();
                        contSmall = 1.0;
                        break;
                    case kPressure:
                        if ( fContInc < 1.0 ) {
                            fprintf( stderr, "not yet implemented feature kPressure && fContInc < 1.0\n" );
                        }
                        contVar = -1;
                        contSmall = 1.0;
                        fContinSide = kNoSide;
                        break;
                    default:
                        fprintf( stderr, "%s\n", "#something wrong in T1DFlame<Species>::PostIter"  );
                        exit(2);
                }

                if ( fContinSide == kBothSides && fContinType != kMomentum ) {
                    Double *bcLeft = grid->GetFine()->GetBcLeft()->vec;
                    Double *bcRight = grid->GetFine()->GetBcRight()->vec;

                    if ( bcLeft[contVar] != bcRight[contVar] ) {
                        fprintf( stderr, "%s%s%s\n", "#error: if ContinuationSide is BothSides, boundary " 
                                , NEWL, "        conditions for left and right side have to be the same" 
                               );
                        exit(2);
                    }
                }

                Double	*yBound;
                Double	*bcBound;
                if ( fContinSide == kLeft || fContinSide == kBothSides ) {
                    yBound = grid->GetFine()->GetYLeft()->vec;
                    bcBound = grid->GetFine()->GetBcLeft()->vec;
                }
                else if ( fContinSide == kRight ) {
                    yBound = grid->GetFine()->GetYRight()->vec;
                    bcBound = grid->GetFine()->GetBcRight()->vec;
                }

                if ( fContinType == kPolytrope ) {
                    // save temperature in pConv;
                    pConv = yBound[contVar];
                }

                if ( fContinType != kPressure ) {
                    contConv = bcBound[contVar];
                }
                else {
                    contConv = GetPressure();
                }

                if ( contNotConv < 0.0 ) {
                    if ( fContinType == kPressure ) {
                        if ( GetPressure() * fContInc < fContBound ) {
                            fSolver->ReInit();
                            SetPressure( MIN( GetPressure() * fContInc
                                        , fContBound ) );
                            fprintf( stderr, "pressure is now %g bar\n", GetPressure()*1.0e-5 );
                        }
                        else {
                            bt->SetLeaveContin();
                        }
                    }
                    else if ( fContinType == kMomentum ) {
                        if ( bcBound[contVar] * fContInc < fContBound * fContInc ) {
                            fSolver->ReInit();
                            Double	*yBoundLeft = grid->GetFine()->GetYLeft()->vec;
                            Double	*bcBoundLeft = grid->GetFine()->GetBcLeft()->vec;
                            Double	*yBoundRight = grid->GetFine()->GetYRight()->vec;
                            Double	*bcBoundRight = grid->GetFine()->GetBcRight()->vec;
                            bcBoundLeft[contVar] *= 1.0 + fContInc;
                            bcBoundRight[contVar] *= 1.0 + fContInc;
                            yBoundLeft[contVar] = bcBoundLeft[contVar];
                            yBoundRight[contVar] = bcBoundRight[contVar];
                            fprintf( stderr, "%s%g\n"
                                    , "V at left boundary is now ", yBoundLeft[contVar] );
                            fprintf( stderr, "%s%g\n"
                                    , "V at right boundary is now ", yBoundRight[contVar] );
                        }
                    }
                    else if ( fContinType == kSeshaTemp ) {
                        if ( bcBound[contVar] * fContInc < fContBound * fContInc ) {
                            fSolver->ReInit();
                            Double	*yBoundLeft = grid->GetFine()->GetYLeft()->vec;
                            Double	*bcBoundLeft = grid->GetFine()->GetBcLeft()->vec;
                            Double	*yBoundRight = grid->GetFine()->GetYRight()->vec;
                            Double	*bcBoundRight = grid->GetFine()->GetBcRight()->vec;
                            Double oldTemp = yBound[contVar];
                            int vVar = GetOffsetVVelocity();

                            yBound[contVar] = bcBound[contVar] = MIN( ( bcBound[contVar] + fContInc ) * fContInc
                                    , fContBound * fContInc ) / fContInc;
                            Double  rhoLeft = fProperties->GetDensity()->vec[kPrev];
                            Double  rhoRight = fProperties->GetDensity()->vec[bt->GetCurrentGridPoints()];

                            bcBoundRight[vVar] *= oldTemp / yBound[contVar];
                            yBoundRight[vVar] = bcBoundRight[vVar];

                            yBoundLeft[vVar] = sqrt( yBoundRight[vVar] * yBoundRight[vVar] / rhoRight * rhoLeft );
                            bcBoundLeft[vVar] = yBoundLeft[vVar];

                            fprintf( stderr, "%s%s%s%s%g\n"
                                    , GetVariableNames()[contVar], " at " 
                                    , ( ( fContinSide == kRight ) ? "right" : "left" )
                                    , " boundary is now ", bcBound[contVar]  );
                            fprintf( stderr, "MomentumLeft is now %g, momentum right is now %g\n", yBoundLeft[vVar] * yBoundLeft[vVar] / rhoLeft
                                    , bcBoundRight[vVar] *bcBoundRight[vVar] / rhoRight );
                        }
                        else {
                            bt->SetLeaveContin();
                        }

                    }
                    else {
                        if ( bcBound[contVar] * fContInc < fContBound * fContInc ) {
                            fSolver->ReInit();
                            yBound[contVar] = bcBound[contVar] = MIN( ( bcBound[contVar] + fContInc ) * fContInc
                                    , fContBound * fContInc ) / fContInc;
                            fprintf( stderr, "%s%s%s%s%g\n"
                                    , GetVariableNames()[contVar], " at " 
                                    , ( ( fContinSide == kRight ) ? "right" : "left" )
                                    , " boundary is now ", bcBound[contVar]  );
                        }
                        else {
                            bt->SetLeaveContin();
                        }
                    }
                }
                else {
                    if ( fContinType == kPressure ) {
                        fSolver->ReInit();
                        SetPressure( contNotConv );
                        fprintf( stderr, "pressure is now %g bar\n", GetPressure()*1.0e-5 );
                        contNotConv = -1.0;
                    }
                    else if ( fContinType == kMomentum ) {
                        if ( bcBound[contVar] * fContInc < fContBound * fContInc ) {
                            fSolver->ReInit();
                            Double	*yBoundLeft = grid->GetFine()->GetYLeft()->vec;
                            Double	*bcBoundLeft = grid->GetFine()->GetBcLeft()->vec;
                            Double	*yBoundRight = grid->GetFine()->GetYRight()->vec;
                            Double	*bcBoundRight = grid->GetFine()->GetBcRight()->vec;
                            fContInc = 0.5 * fContInc;
                            bcBoundLeft[contVar] *= 1.0 + fContInc;
                            bcBoundRight[contVar] *= 1.0 + fContInc;
                            yBoundLeft[contVar] = bcBoundLeft[contVar];
                            yBoundRight[contVar] = bcBoundRight[contVar];
                            fprintf( stderr, "%s%g\n"
                                    , "V at left boundary is now ", yBoundLeft[contVar] );
                            fprintf( stderr, "%s%g\n"
                                    , "V at right boundary is now ", yBoundRight[contVar] );
                        }
                    }
                    else if ( fContinType == kSeshaTemp ) {
                        fSolver->ReInit();

                        Double	*yBoundLeft = grid->GetFine()->GetYLeft()->vec;
                        Double	*bcBoundLeft = grid->GetFine()->GetBcLeft()->vec;
                        Double	*yBoundRight = grid->GetFine()->GetYRight()->vec;
                        Double	*bcBoundRight = grid->GetFine()->GetBcRight()->vec;
                        Double oldTemp = yBound[contVar];
                        int vVar = GetOffsetVVelocity();

                        bcBound[contVar] = yBound[contVar] = 0.5 * ( contNotConv + yBound[contVar] );
                        fContInc = 0.5 * fContInc;

                        Double  rhoLeft = fProperties->GetDensity()->vec[kPrev];
                        Double  rhoRight = fProperties->GetDensity()->vec[bt->GetCurrentGridPoints()];

                        bcBoundRight[vVar] *= oldTemp / yBound[contVar];
                        yBoundRight[vVar] = bcBoundRight[vVar];

                        yBoundLeft[vVar] = sqrt( yBoundRight[vVar] * yBoundRight[vVar] / rhoRight * rhoLeft );
                        bcBoundLeft[vVar] = yBoundLeft[vVar];

                        fprintf( stderr, "%s%s%s%s%g\n"
                                , GetVariableNames()[contVar], " at " 
                                , ( ( fContinSide == kRight ) ? "right" : "left" )
                                , " boundary is now ", bcBound[contVar]  );
                        fprintf( stderr, "MomentumLeft is now %g, momentum right is now %g\n", yBoundLeft[vVar] * yBoundLeft[vVar] / rhoLeft
                                , bcBoundRight[vVar] *bcBoundRight[vVar] / rhoRight );

                        contNotConv = -1.0;
                    }
                    else {
                        fSolver->ReInit();
                        bcBound[contVar] = yBound[contVar] = contNotConv;
                        fprintf( stderr, "%s%s%s%s%g\n", GetVariableNames()[contVar], " at " 
                                , ( ( fContinSide == kRight ) ? "right" : "left" )
                                , " boundary is now ", bcBound[contVar]  );
                        contNotConv = -1.0;
                    }
                }

                if ( fContinSide == kBothSides && !bt->GetLeaveContin()
                        && fContinType != kPolytrope && fContinType != kMomentum ) {
                    Double	*yRight = grid->GetFine()->GetYRight()->vec;
                    Double	*bcRight = grid->GetFine()->GetBcRight()->vec;

                    yRight[contVar] = bcRight[contVar] = yBound[contVar];
                    fprintf( stderr, "%s%s%g\n", GetVariableNames()[contVar] 
                            , " at right boundary is now ", bcBound[contVar]  );
                }
                if ( fContinType == kPolytrope ) {
                    Double	polExp = 1.33;

                    fPressure->vec[fPressure->len] *= 
                        pow( bcBound[contVar] / pConv, polExp / ( polExp - 1.0 ) );
                    fprintf( stderr, "pressure is now %g bar\n", GetPressure()*1.0e-5 );
                    // init pConv
                    pConv = -1.0;
                }
            }
        }
    }
    else {
        bt->SetLeaveContin();
        if ( toSpeciesIndex-firstSpecies >= 0 ) {
            int 	fromSpeciesIndex = GetFromSpeciesIndex() + firstSpecies;
            char	**names = GetSpecies()->GetNames();
            Double	*yBound;
            Double	*bcBound;

            if ( fromSpeciesIndex-firstSpecies < 0 ) {
                fromSpeciesIndex = GetFuelIndex() + firstSpecies;
            }

            if ( fContInc == 0.0 || fContinSide == kNoSide ) {
                fprintf( stderr, "%s%s%s\n", "#warning: increment of continuation is zero", NEWL, 
                        "          or no continuation side specified"  );
                return;
            }

            if ( fContinSide == kLeft ) {
                yBound = grid->GetFine()->GetYLeft()->vec;
                bcBound = grid->GetFine()->GetBcLeft()->vec;
            }
            else {
                yBound = grid->GetFine()->GetYRight()->vec;
                bcBound = grid->GetFine()->GetBcRight()->vec;
            }

            if ( fContinType == kPressure ) {
                fContInc = sqrt(fContInc);
            }
            else {
                fContInc *= 0.5;
            }

            if ( fContInc >= 0.005 ) {
                fSolver->ReInit();

                bcBound[toSpeciesIndex] = yBound[toSpeciesIndex] 
                    = MAX( MIN( 1.0, yBound[toSpeciesIndex] - fContInc ), 1.0e-60 );
                bcBound[fromSpeciesIndex] = yBound[fromSpeciesIndex] 
                    = MIN( MAX( 1.0e-60, yBound[fromSpeciesIndex] + fContInc ), 1.0 );
                fprintf( stderr, "%s%s%s%g%s%s%s%g\n"
                        , "Y_", names[toSpeciesIndex-firstSpecies], " is now "
                        , bcBound[toSpeciesIndex], " and Y_"
                        , names[fromSpeciesIndex-firstSpecies], " is now ", bcBound[fromSpeciesIndex] );
            }
        }
        else {
            if ( aVec && aVec->len < aVec->phys_len ) {
                Double	interStrainRate = aConv + ( GetStrainRate() - aConv ) * 0.5;
                if ( aConv >= 0.0 && fabs( interStrainRate - aConv ) / aConv >= 0.001 ) {
                    aNotConv = aVec->vec[aVec->len];
                    aVec->vec[aVec->len] = interStrainRate;
                    fSolver->ReInit();
                    fprintf( stderr, "%s%g\n", "strainRate is now ", GetStrainRate()  );
                }
            }
            if ( fPressure && fPressure->len < fPressure->phys_len ) {
                Double	interPressure = pConv + ( GetPressure() - pConv ) * 0.5;
                if ( pConv >= 0.0 && fabs( interPressure - pConv ) >= 0.1 ) {
                    pNotConv = fPressure->vec[fPressure->len];
                    fPressure->vec[fPressure->len] = interPressure;
                    fSolver->ReInit();
                    fprintf( stderr, "pressure is now %g bar\n", GetPressure()*1.0e-5 );
                }
            }
            if ( fContinType == kTemperature || fContinType == kVelocity || fContinType == kPressure
                    || fContinType == kMomentum || fContinType == kSeshaTemp ) {
                if ( fContInc == 0.0 || ( fContinSide == kNoSide && fContinType != kPressure )  ) {
                    fprintf( stderr, "%s%s%s\n", "#warning: increment of continuation is zero", NEWL, 
                            "          or no continuation side specified"  );
                    return;
                }

                switch( fContinType ) {
                    case kTemperature:
                        contVar = GetOffsetTemperature();
                        contSmall = 1.0;
                        break;
                    case kVelocity:
                        contVar = GetOffsetVVelocity();
                        contSmall = 0.001;
                        break;
                    case kMomentum:
                        contVar = GetOffsetVVelocity();
                        contSmall = 0.25 * ( fSolver->bt->GetRight() - fSolver->bt->GetLeft() ) // l/4
                            * sqrt( fProperties->GetDensity()->vec[kPrev] // sqrt( rhoLeft/rhoRight )
                                    / fProperties->GetDensity()->vec[bt->GetCurrentGridPoints()] ) 
                            * 0.5; // delta_aSeshaSmall
                        fContinSide = kBothSides;
                        break;
                    case kSeshaTemp:
                        contVar = GetOffsetTemperature();
                        contSmall = 1.0;
                        break;
                    case kPressure:
                        if ( fContInc < 1.0 ) {
                            fprintf( stderr, "not yet implemented feature kPressure && fContInc < 1.0\n" );
                        }
                        contVar = -1;
                        contSmall = 1.0;
                        fContinSide = kNoSide;
                        break;
                    default:
                        fprintf( stderr, "%s\n", "#something wrong in T1DFlame<Species>::PostIter"  );
                        exit(2);
                }

                Double	*yBound;
                Double	*bcBound;
                if ( fContinSide == kLeft || fContinSide == kBothSides ) {
                    yBound = grid->GetFine()->GetYLeft()->vec;
                    bcBound = grid->GetFine()->GetBcLeft()->vec;
                }
                else if ( fContinSide == kRight ) {
                    yBound = grid->GetFine()->GetYRight()->vec;
                    bcBound = grid->GetFine()->GetBcRight()->vec;
                }


                fContInc *= 0.5;
                Double	interCont;
                if ( fContinType == kPressure ) {
                    interCont = sqrt( GetPressure() * contConv );
                }
                else if ( fContinType == kMomentum ) {
                    interCont = contConv * ( 1.0 + fContInc );
                }
                else {
                    interCont = contConv + ( bcBound[contVar] - contConv ) * 0.5;
                }

                fprintf( stderr, "interCont = %g\tcontConv = %g\tfContInc = %g\tcontSmall = %g\n"
                        , interCont, contConv, fContInc, contSmall );

                if ( contConv >= 0.0 
                        && fabs( interCont - contConv ) >= contSmall ) {
                    if ( fContinType == kPressure ) {
                        contNotConv = GetPressure();
                        SetPressure( interCont );
                        fSolver->ReInit();
                        fprintf( stderr, "pressure is now %g bar\n", GetPressure()*1.0e-5 );
                    }
                    else if ( fContinType == kSeshaTemp ) {
                        Double oldTemp = yBound[contVar];
                        contNotConv = bcBound[contVar];
                        bcBound[contVar] = yBound[contVar] = interCont;
                        fSolver->ReInit();

                        Double	*yBoundLeft = grid->GetFine()->GetYLeft()->vec;
                        Double	*bcBoundLeft = grid->GetFine()->GetBcLeft()->vec;
                        Double	*yBoundRight = grid->GetFine()->GetYRight()->vec;
                        Double	*bcBoundRight = grid->GetFine()->GetBcRight()->vec;
                        int vVar = GetOffsetVVelocity();

                        Double  rhoLeft = fProperties->GetDensity()->vec[kPrev];
                        Double  rhoRight = fProperties->GetDensity()->vec[bt->GetCurrentGridPoints()];

                        bcBoundRight[vVar] *= oldTemp / yBound[contVar];
                        yBoundRight[vVar] = bcBoundRight[vVar];

                        yBoundLeft[vVar] = sqrt( yBoundRight[vVar] * yBoundRight[vVar] / rhoRight * rhoLeft );
                        bcBoundLeft[vVar] = yBoundLeft[vVar];

                        fprintf( stderr, "%s%s%s%s%g\n"
                                , GetVariableNames()[contVar], " at " 
                                , ( ( fContinSide == kRight ) ? "right" : "left" )
                                , " boundary is now ", bcBound[contVar]  );
                        fprintf( stderr, "MomentumLeft is now %g, momentum right is now %g\n", yBoundLeft[vVar] * yBoundLeft[vVar] / rhoLeft
                                , bcBoundRight[vVar] *bcBoundRight[vVar] / rhoRight );
                    }
                    else {
                        contNotConv = bcBound[contVar];
                        bcBound[contVar] = yBound[contVar] = interCont;
                        fSolver->ReInit();
                        fprintf( stderr, "%s%s%s%s%g\n", GetVariableNames()[contVar], " at " 
                                , ( ( fContinSide == kRight ) ? "right" : "left" )
                                , " boundary is now ", bcBound[contVar]  );
                    }
                    if ( fContinSide == kBothSides ) {
                        Double	*yRight = grid->GetFine()->GetYRight()->vec;
                        Double	*bcRight = grid->GetFine()->GetBcRight()->vec;
                        if ( fContinType != kMomentum ) {
                            yRight[contVar] = bcRight[contVar] = yBound[contVar];
                            fprintf( stderr, "%s%s%g\n", GetVariableNames()[contVar] 
                                    , " at boundaries is now ", bcBound[contVar]  );
                        }
                        else {
                            yRight[contVar] = bcRight[contVar] = -sqrt( bcBound[contVar] * bcBound[contVar]
                                    / fProperties->GetDensity()->vec[kPrev] 
                                    * fProperties->GetDensity()->vec[bt->GetCurrentGridPoints()] );
                            fprintf( stderr, "%s%s%g\n", GetVariableNames()[contVar] 
                                    , " at right boundary is now ", bcRight[contVar]  );
                        }
                    }
                }
                else {
                    fprintf( stderr, "leave continuation\n" );
                }
            }
        }
    }
}

    template<class Species>
void T1DFlame<Species>::SetInitialBC( TGridPtr grid, TInputDataPtr inp )
{
    int					i;
    Double				mixMolarMass;
    int					nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    //SpeciesPtr			species = inp->GetSpecies();	// ###
    BoundaryInputPtr	left = inp->leftBoundary;
    BoundaryInputPtr	right = inp->rightBoundary;
    int					inpVOffset = inp->fVVelocityOffset;
    int					inpUOffset = inp->fUVelocityOffset;
    int					inpTOffset = inp->fTemperatureOffset;
    int					fFirstSpecies = GetOffsetFirstSpecies();
    int					fUVelocity = GetOffsetUVelocity();
    int					fVVelocity = GetOffsetVVelocity();
    int					fTemperature = GetOffsetTemperature();
    int					*speciesIndexLeft = NULL;
    int					*speciesIndexRight = NULL;
    int					leftSpecifiedSpecies = left->fSpecifiedSpeciesBCs;
    int					rightSpecifiedSpecies = right->fSpecifiedSpeciesBCs;
    int					*bcFlagLeft = grid->GetBcFlagLeft();
    int					*bcFlagRight = grid->GetBcFlagRight();
    Double				*yleft = grid->GetYLeft()->vec;
    Double				*yright = grid->GetYRight()->vec;
    Double				*bcLeft = grid->GetBcLeft()->vec;
    Double				*bcRight = grid->GetBcRight()->vec;
    //	Double				*massFracsLeft = fSolMassFracs->mat[kPrev];
    //	Double				*massFracsRight = fSolMassFracs->mat[kPrev];

    //	allocate memory for speciesIndex
    speciesIndexLeft = new int[left->fSpecifiedSpeciesBCs];
    if ( !speciesIndexLeft ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );
    speciesIndexRight = new int[right->fSpecifiedSpeciesBCs];
    if ( !speciesIndexRight ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );

    //	set speciesIndex
    for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
        speciesIndexLeft[i] = inp->FindSpecies( left->speciesName[i] );
    }
    for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
        speciesIndexRight[i] = inp->FindSpecies( right->speciesName[i] );
    }

    // set fMixtureSpecification
    SetMixtureSpecificationLeft( left->fMixtureSpecification );
    SetMixtureSpecificationRight( right->fMixtureSpecification );

    // set BCFlags
    bcFlagLeft[fVVelocity] = left->fBcFlag[inpVOffset];
    bcFlagLeft[fUVelocity] = left->fBcFlag[inpUOffset];
    bcFlagLeft[fTemperature] = left->fBcFlag[inpTOffset];
    bcFlagRight[fVVelocity] = right->fBcFlag[inpVOffset];
    bcFlagRight[fUVelocity] = right->fBcFlag[inpUOffset];
    bcFlagRight[fTemperature] = right->fBcFlag[inpTOffset];
    for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
        bcFlagLeft[i] = left->fBcFlagSpecies;
        bcFlagRight[i] = right->fBcFlagSpecies;
    }

    // set value
    yleft[fVVelocity] = left->fValue[inpVOffset];
    yleft[fUVelocity] = left->fValue[inpUOffset];
    yleft[fTemperature] = left->fValue[inpTOffset];
    yright[fVVelocity] = right->fValue[inpVOffset];
    yright[fUVelocity] = right->fValue[inpUOffset];
    yright[fTemperature] = right->fValue[inpTOffset];

    bcLeft[fVVelocity] = left->fValue[inpVOffset];
    bcLeft[fUVelocity] = left->fValue[inpUOffset];
    bcLeft[fTemperature] = left->fValue[inpTOffset];
    bcRight[fVVelocity] = right->fValue[inpVOffset];
    bcRight[fUVelocity] = right->fValue[inpUOffset];
    bcRight[fTemperature] = right->fValue[inpTOffset];

    for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
        if ( speciesIndexLeft[i] >= nSpeciesInSystem ) {
            fprintf( stderr, "%s%s%s\n", "#warning: value at left boundary for steady state species '"
                    , fSpecies.GetNames()[i], "' specified, make no use of it"  );
        }
        else {
            yleft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
            bcLeft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
        }
    }

    if ( left->fMixtureSpecification == kMolarFraction ) {
        // first compute molar mass of mixture
        for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
            mixMolarMass += fSpecies.GetMolarMass()->vec[i] * yleft[i+fFirstSpecies];	// ###
        }
        // compute massfractions
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            yleft[i+fFirstSpecies] *= fSpecies.GetMolarMass()->vec[i] / mixMolarMass;	// ###
            bcLeft[i+fFirstSpecies] = yleft[i+fFirstSpecies];
        }
        for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
            bcFlagLeft[i] = kMassFraction;
        }
    }

    for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
        if ( speciesIndexRight[i] >= nSpeciesInSystem ) {
            fprintf( stderr, "%s%s%s\n", "#warning: value at right boundary for steady state species '"
                    , fSpecies.GetNames()[i], "' specified, make no use of it"  );
        }
        else {
            yright[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
            bcRight[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
        }
    }
    if ( right->fMixtureSpecification == kMolarFraction ) {
        // first compute molar mass of mixture
        for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
            mixMolarMass += fSpecies.GetMolarMass()->vec[i] * yright[i+fFirstSpecies];	// ###
        }
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            yright[i+fFirstSpecies] *= fSpecies.GetMolarMass()->vec[i] / mixMolarMass;	// ###
            bcRight[i+fFirstSpecies] = yright[i+fFirstSpecies];
        }
        for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
            bcFlagRight[i] = kMassFraction;
        }
    }
    delete speciesIndexRight;
    delete speciesIndexLeft;
}

    template<class Species>
void T1DFlame<Species>::CheckBC( void )
{
    int			i, iOff;
    int			firstSpecies = GetOffsetFirstSpecies();
    int			nOfSpecies = fSpecies.GetNOfSpecies();
    int			nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    Double		sumYLeft = 0.0;
    Double		sumYRight = 0.0;
    TNewtonPtr	bt = fSolver->bt;
    TGridPtr	grid = bt->GetGrid()->GetFine();
    Double		*yLeft = grid->GetYLeft()->vec;
    Double		*yRight = grid->GetYRight()->vec;
    int			nGridPoints = grid->GetNGridPoints();
    Double		*massFracLeft = fSolMassFracs->mat[kPrev];
    Double		*massFracRight = fSolMassFracs->mat[nGridPoints];

    if ( fClipNegativeConcs ) {
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            iOff = i + firstSpecies;
            yLeft[iOff] = MAX( 1.0e-60, yLeft[iOff] );
            yRight[iOff] = MAX( 1.0e-60, yRight[iOff] );

        }
    }
    for ( i = 0; i < nOfSpecies; ++i ) {
        if ( fClipNegativeConcs ) {
            massFracLeft[i] = MAX( 1.0e-60, massFracLeft[i] );
            massFracRight[i] = MAX( 1.0e-60, massFracRight[i] );
        }
        sumYLeft += massFracLeft[i];
        sumYRight += massFracRight[i];
    }

    if ( fabs( sumYLeft - 1.0 ) > 1.0e-3 ) {
        fprintf( stderr, "%s%g\n", "#warning: sum of Y_i at left boundary is ", sumYLeft  );
    }
    if ( fabs( sumYRight - 1.0 ) > 1.0e-3 ) {
        fprintf( stderr, "%s%g\n", "#warning: sum of Y_i at right boundary is ", sumYRight  );
    }
}

    template<class Species>
void T1DFlame<Species>::SetFlameNode( int k )
{
    int gridPoints = GetSolver()->bt->GetCurrentGridPoints();

    if ( k >= 0 && k < gridPoints ) {
        // std::cout << "fReaction.GetTBConcentrations()->mat[k][0]" << fReaction.GetTBConcentrations()->mat[k][0] << '\n';
        fFlameNode->tBodyConc = fReaction.GetTBConcentrations()->mat[k];
        // std::cout << "fReaction.GetRateCoefficients()->mat[k][0]" << fReaction.GetRateCoefficients()->mat[k][0] << '\n';
        fFlameNode->rateCoeff = fReaction.GetRateCoefficients()->mat[k];
        // std::cout << "fReaction.GetTempReaction()->vec[k]" << fReaction.GetTempReaction()->vec[k] << '\n';
        fFlameNode->tempReaction = &fReaction.GetTempReaction()->vec[k];
        // std::cout << "fReaction.GetPressureReaction()->vec[k]" << fReaction.GetPressureReaction()->vec[k] << '\n';
        fFlameNode->pressureReaction = &fReaction.GetPressureReaction()->vec[k];
        // std::cout << "fReaction.GetCurrRateCoeff()->mat[k][0]" << fReaction.GetCurrRateCoeff()->mat[k][0] << '\n';
        fFlameNode->currRateCoeff = fReaction.GetCurrRateCoeff()->mat[k];
        // std::cout << "&fReaction.GetKNewerThanW()[k]" << &fReaction.GetKNewerThanW()[k] << '\n';
        fFlameNode->kNewerThanW = &fReaction.GetKNewerThanW()[k];
        // std::cout << "fReaction.GetReactionRate()->mat[k][0]" << fReaction.GetReactionRate()->mat[k][0] << '\n';
        fFlameNode->reactionRate = fReaction.GetReactionRate()->mat[k];
        // std::cout << "fReaction.GetYReaction()->mat[k][0]" << fReaction.GetYReaction()->mat[k][0] << '\n';
        fFlameNode->YReaction = fReaction.GetYReaction()->mat[k];
        // std::cout << "fReaction.GetCurrReacRate()->mat[k][0]" << fReaction.GetCurrReacRate()->mat[k][0] << '\n';
        fFlameNode->currReacRate = fReaction.GetCurrReacRate()->mat[k];
        // std::cout << "fReaction.GetDmdy()->tensor[k][0][0]" << fReaction.GetDmdy()->tensor[k][0][0] << '\n';
        fFlameNode->dMdY = fReaction.GetDmdy()->tensor[k];
        // std::cout << "fSpecies.GetProductionRate()->mat[k][0]" << fSpecies.GetProductionRate()->mat[k][0] << '\n';
        fFlameNode->productionRate = fSpecies.GetProductionRate()->mat[k];
        if ( fProperties->GetRadiation() ) {
            fFlameNode->radiation = &fProperties->GetRadiation()->GetRadiation()->vec[k];
            // fFlameNode->kappa = &fProperties->GetRadiation()->GetKappa()->vec[k];
        }
        fFlameNode->diffusivityNext = fSpecies.GetDiffusivity()->mat[k+1];
        fFlameNode->diffusivityPrev = fSpecies.GetDiffusivity()->mat[k-1];
    }
    else if ( k >= 0 ) {
        fFlameNode->diffusivityPrev = fSpecies.GetDiffusivity()->mat[k-1];
    }
    else if ( k < fReaction.GetTBConcentrations()->phys_cols ) {
        fFlameNode->diffusivityNext = fSpecies.GetDiffusivity()->mat[k+1];
    }

    fFlameNode->tempProp = &fSpecies.GetTempProp()->vec[k];
    fFlameNode->pressureProp = &fSpecies.GetPressureProp()->vec[k];
    fFlameNode->diffCorr = &GetDiffCorr()->vec[k];
    fFlameNode->viscosityInf = fProperties->GetViscosity()->vec[gridPoints];
    fFlameNode->rhoInf = fProperties->GetDensity()->vec[gridPoints];

    fFlameNode->viscosity = fSpecies.GetViscosity()->mat[k];
    fFlameNode->heatCapacity = fSpecies.GetHeatCapacity()->mat[k];
    fFlameNode->conductivity = fSpecies.GetConductivity()->mat[k];
    fFlameNode->enthalpy = fSpecies.GetEnthalpy()->mat[k];
    fFlameNode->diffusivity = fSpecies.GetDiffusivity()->mat[k];
    fFlameNode->diffTherm = &fSpecies.GetDiffTherm()->mat[k];
    fFlameNode->savedY = fSpecies.GetSavedY()->mat[k];
    fFlameNode->savedDeltaiY = fSpecies.GetSavedDeltaiY()->mat[k];
    fFlameNode->sumDiff = fSpecies.GetSumDiff()->mat[k];
    fFlameNode->deltaI = fSpecies.GetDeltaI()->mat[k];
    fFlameNode->tempProp = &fSpecies.GetTempProp()->vec[k];
    fFlameNode->pressureProp = &fSpecies.GetPressureProp()->vec[k];
    fFlameNode->binDij = &fSpecies.GetBinDij()->tensor[k];
    fFlameNode->GijOverWj = fSpecies.GetGijOverWj()->tensor[k];
    fFlameNode->OmegaDOverDCoeff = fSpecies.GetOmegaDOverDCoeff()->tensor[k];
    if ( fThermoDiffusion ) {
        fFlameNode->DThermConst = fSpecies.GetDThermConst()->tensor[k];
    }

    fFlameNode->mixViscosity = &fProperties->GetViscosity()->vec[k];
    fFlameNode->mixDensity = &fProperties->GetDensity()->vec[k];
    fFlameNode->mixConductivity = &fProperties->GetConductivity()->vec[k];
    fFlameNode->mixHeatCapacity = &fProperties->GetHeatCapacity()->vec[k];
    fFlameNode->mixMolarMass = &fProperties->GetMolarMass()->vec[k];

    fFlameNode->Y = &fSolMassFracs->mat[k];
    fFlameNode->temp = &fSolTemp->vec[k];
    if ( fSoot ) {
        fFlameNode->Pij = fSoot->GetPij()->tensor[k];
        fFlameNode->sumPi = fSoot->GetSumPi()->mat[k];
        fFlameNode->pahMoments = fSoot->GetPAHMoments()->mat[k];
        fFlameNode->moments = fSoot->GetMoments()->mat[k];
        fFlameNode->pahReactionRate = fSoot->GetReactionRate()->mat[k];
        fFlameNode->diffSoot = &fSoot->GetSootDiff()->vec[k];
        if ( k >= 0 && k < gridPoints ) {
            fFlameNode->sootSource = fSoot->GetSource()->mat[k];
            fFlameNode->dMdx = fSoot->GetDMdx()->tensor[k];
        }
    }
}

    template<class Species>
void T1DFlame<Species>::SetBVPInput( TBVPInputPtr input )
{
    input->fWriteBT = fInputData->fWriteBT;
    input->fWriteResiduum = fInputData->fWriteResiduum;
    input->fWatchGridding = fInputData->fWatchGridding;
    input->fWriteEverySolution = fInputData->fWriteEverySolution;
    input->fOutputPath = fInputData->fOutputPath;

    input->fWriteFullRes = fInputData->fWriteFullRes;
    input->fUseModifiedNewton = fInputData->fUseModifiedNewton;
    input->fUseNumericalJac = fInputData->fUseNumericalJac;
    input->fUseSecOrdJac = fInputData->fUseSecOrdJac;

    // TNewton
    input->fNVariables = fInputData->fNVariables;;
    input->fInitialEquations = fInputData->fInitialEquations;
    input->fMaxGridPoints = fInputData->fMaxGridPoints;
    input->fInitialGridPoints = fInputData->fInitialGridPoints;
    input->fDampFlag = fInputData->fDampFlag;
    input->fTimeDepFlag = fInputData->fTimeDepFlag;
    input->fContinFlag = fInputData->fContinFlag;
    input->fDeltaNewGrid = fInputData->fDeltaNewGrid;
    input->fTolRes = fInputData->fTolRes;
    input->fTolDy = fInputData->fTolDy;
    input->fMaxIter = fInputData->fMaxIter;
    input->fLeft = fInputData->fLeft;
    input->fRight = fInputData->fRight;

    // TAdaptiveGrid
    input->fOneSolOneGrid = fInputData->fOneSolOneGrid;
    input->fR = fInputData->fR;
    input->fQ = fInputData->fQ;
    input->fStart = fInputData->fStart;
    input->fEnd = fInputData->fEnd;
    input->fAlpha = fInputData->fAlpha;

    // TDamp
    input->fLambdaMin = fInputData->fLambdaMin;

    // TTime
    if ( fInputData->fDeltaTStart > 0.0 ) {
        input->fDeltaTStart = fInputData->fDeltaTStart;
    }
    else {
        input->fDeltaTStart = 1.0;
    }

    // TContinuation
    input->fContSteps = fInputData->fContSteps;
}

    template<class Species>
void T1DFlame<Species>::ReadStartProfiles( TInputDataPtr inp )
{
    StartProfilePtr	sp = NULL;
    char 	*insp = inp->fStartProfileFile;
    FILE	*fpS = NULL;
    char	*fileName = GetFullPath( insp, kFileName );

    sp = new StartProfile;	

    if ( !insp || ( fpS = fopen( fileName, "r" ) ) == NULL ) {
        fprintf( stderr, "#error: can't open input file for startprofiles '%s'\n", fileName );
        exit( 2 );
    }
    else {
        fprintf( stderr, "use startprofiles file '%s'\n", fileName );
    }
    delete fileName;

    ::ReadStartProfiles( sp, fpS );
    SetInitialValues( fInputData, sp ); // initial values of coarse grid are set during gridgeneration
    CleanReadStartProfile();
    delete sp;
    fclose( fpS );
}

    template<class Species>
void T1DFlame<Species>::CheckInitialGuess( void )
{
    int			i, k;
    int 		gridPoints = fSolver->bt->GetCurrentGridPoints();
    int 		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    Double		sumY; 
    Double		*temp = fSolTemp->vec;
    Double		**Y = fSolMassFracs->mat;
    int			speciesEq;
    int			firstSpeciesOffset = GetOffsetFirstSpecies();
    int			temperatureOffset = GetOffsetTemperature();
    TGridPtr	grid = fSolver->bt->GetGrid()->GetCurrentGrid();
    Double		**y = grid->GetY()->mat;
    Double		error = 1.0e-3;
    Double		maxErrSumYi = 0.0;
    int			errPoint;

    for ( k = 0; k < gridPoints; ++k ) {
        temp[k] = MAX( temp[k], 10.0 );
        for ( i = 0, sumY = 0.0; i < nSpeciesInSystem; ++i ) {
            if ( fClipNegativeConcs ) {
                Y[k][i] = MAX( Y[k][i], 1.0e-60 );
            }
            Y[k][i] = MIN( Y[k][i], 1.0 );
            sumY += Y[k][i];
        }
        y[k][temperatureOffset] = MAX( y[k][temperatureOffset], 10.0 );
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            speciesEq = i + firstSpeciesOffset;
            if ( fClipNegativeConcs ) {
                y[k][speciesEq] = MAX( y[k][speciesEq], 1.0e-60 );
            }
            y[k][speciesEq] = MIN( y[k][speciesEq], 1.0 );
        }
        if ( fabs( sumY - 1.0 ) > error ) {
            error = fabs( sumY - 1.0 );
            maxErrSumYi = sumY;
            errPoint = k;
        }
    }
    if ( error > 1.0e-3 ) {
        fprintf( stderr, "#warning: maximum error in mass conservation at gridpoint no. %d: sumYi = %g\n", errPoint, maxErrSumYi );
    }
}

    template<class Species>
Double T1DFlame<Species>::GetZStoich( void )
{

    if ( GetNFuels() > 1 ) {
        return GetZStoich_mf();
    }

    int			oxIndex = fInputData->fOxIndex;
    int			fuelIndex = GetFuelIndex();
    int			firstSpecies = GetOffsetFirstSpecies();
    Double		nuOx = GetNu( GetInputData()->fGlobalReaction, fSpecies.GetNames()[oxIndex] );
    Double		nuFuel = GetNu( GetInputData()->fGlobalReaction, fSpecies.GetNames()[fuelIndex] );
    Double		molarMassOx = GetSpecies()->GetMolarMass()->vec[oxIndex];
    Double		molarMassFuel = GetSpecies()->GetMolarMass()->vec[fuelIndex];
    Double		nu;
    TGridPtr	grid = GetSolver()->bt->GetGrid()->GetFine();
    Double		Yf1, YO22;


    if ( grid->GetYRight()->vec[firstSpecies+fuelIndex] > grid->GetYLeft()->vec[firstSpecies+fuelIndex] ) {
        YO22 = grid->GetYLeft()->vec[firstSpecies+oxIndex];
        Yf1 = grid->GetYRight()->vec[firstSpecies+fuelIndex];
    }
    else {
        YO22 = grid->GetYRight()->vec[firstSpecies+oxIndex];
        Yf1 = grid->GetYLeft()->vec[firstSpecies+fuelIndex];
    }

    if ( nuOx == -1 ) {
        fprintf( stderr, "%s\n", "error: there is no oxidizer in reaction no. '0'"  );
        exit(2);
    }
    if ( nuFuel == -1 ) {
        fprintf( stderr, "%s\n", "error: there is no fuel in reaction no. '0'"  );
        exit(2);
    }

    nu = nuOx * molarMassOx / ( nuFuel * molarMassFuel );

    return 1.0 / ( 1.0 + nu * Yf1 / YO22 );
}

    template<class Species>
Double T1DFlame<Species>::GetZStoich_mf( void )
{//(gb)		
    int			n;
    int			*fuelIndex_i = GetFuelIndexVec()->vec;
    int			oxIndex = GetInputData()->FindSpecies( "O2" );
    int			firstSpecies = GetOffsetFirstSpecies();
    Double		nuOx = GetNu( GetInputData()->fGlobalReaction, "O2" );
    char		**names = fSpecies.GetNames();
    Double		nuFuel_i;
    VectorPtr	fMolarMass = fSpecies.GetMolarMass();
    TGridPtr	grid = GetSolver()->bt->GetGrid()->GetFine();
    Double 		sum1 = 0.0;
    Double 		sumY = 0.0;
    Double		nu;
    Double		Yf_i1, YO22;
    Double		molarMassOx, molarMassFuel;

    for (n=0; n < GetNFuels(); ++n) {		
        if ( grid->GetYRight()->vec[firstSpecies+fuelIndex_i[n]] > grid->GetYLeft()->vec[firstSpecies+fuelIndex_i[n]] ) {
            YO22 = grid->GetYLeft()->vec[firstSpecies+oxIndex];
            Yf_i1 = grid->GetYRight()->vec[firstSpecies+fuelIndex_i[n]];
            sumY += Yf_i1;			
        }
        else {
            YO22 = grid->GetYRight()->vec[firstSpecies+oxIndex];
            Yf_i1 = grid->GetYLeft()->vec[firstSpecies+fuelIndex_i[n]];
            sumY += Yf_i1;			
        }
    }

    if ( nuOx == -1 ) {
        fprintf( stderr, "%s\n", "error: there is no oxidizer in reaction no. '0'"  );
        exit(2);
    }
    molarMassOx = fMolarMass->vec[oxIndex];

    for (n=0; n < GetNFuels(); ++n) {
        nuFuel_i = GetNu( GetInputData()->fGlobalReaction, names[fuelIndex_i[n]] );
        molarMassFuel = fMolarMass->vec[fuelIndex_i[n]];

        if ( nuFuel_i == -1 ) {
            fprintf( stderr, "%s\n", "error: there is no fuel in reaction no. '0'"  );
            exit(2);
        }

        sum1 += nuFuel_i * molarMassFuel;			
    }

    nu = nuOx * molarMassOx / ( sum1 );

    return 1.0 / ( 1.0 + nu * sumY / YO22 );

}

    template<class Species>
Double T1DFlame<Species>::ComputeZBilger( Double *Y, Double *YFuelSide, Double *YOxSide )
{
    static const Double	molarMassC = 12.01, 
                 molarMassO = 16.0,
                 molarMassH = 1.008;
    Double				z, zC, zO, zH, zOO, zCF, zHF, zOF, zCO, zHO;
    TInputDataPtr		inp = GetInputData();
    int					CNum = inp->FindAtomIndex( "C" );
    int					HNum = inp->FindAtomIndex( "H" );
    int					ONum = inp->FindAtomIndex( "O" );
    //SpeciesPtr			species = inp->GetSpecies();		// ### composition can be determined by fSpecies.GetComposition()
    Double				nuC, nuH, nuO;
    int					oxIndex = inp->FindSpecies( "O2" );

    if ( oxIndex < 0 ) {
        return -1.0;
    }

    nuC = ( CNum >= 0 ) ? fSpecies.GetComposition()->mat[ GetFuelIndex() ][ CNum ] : 0;	// ###
    nuH = ( HNum >= 0 ) ? fSpecies.GetComposition()->mat[ GetFuelIndex() ][ HNum ] : 0;	// ###
    nuO = ( ONum >= 0 ) ? fSpecies.GetComposition()->mat[ oxIndex ][ ONum ] : 0;		// ###
    zC = GetElementMassFraction( Y, "C", molarMassC );
    zO = GetElementMassFraction( Y, "O", molarMassO );
    zH = GetElementMassFraction( Y, "H", molarMassH );
    zOO = GetElementMassFraction( YOxSide, "O", molarMassO );
    zCF = GetElementMassFraction( YFuelSide, "C", molarMassC );
    zHF = GetElementMassFraction( YFuelSide, "H", molarMassH );

    zOF = GetElementMassFraction( YFuelSide, "O", molarMassO );
    zCO = GetElementMassFraction( YOxSide, "C", molarMassC );
    zHO = GetElementMassFraction( YOxSide, "H", molarMassH );

    if ( nuO == 0 ) { return -1.0; }

    if ( nuC == 0 ) {
        z = ( ( zH - zHO ) / ( nuH * molarMassH ) 
                + 2.0 * ( zOO - zO ) / ( nuO * molarMassO ) )
            / ( ( zHF - zHO ) / ( nuH * molarMassH ) 
                    + 2.0 * ( zOO - zOF ) / ( nuO * molarMassO ) );
    }
    else if ( nuH == 0 ) {
        z = ( ( zC - zCO ) / ( nuC * molarMassC )
                + 2.0 * ( zOO - zO ) / ( nuO * molarMassO ) )
            / ( ( zCF - zCO ) / ( nuC * molarMassC )
                    + 2.0 * ( zOO - zOF ) / ( nuO * molarMassO ) );
    }
    else {
        z = ( ( zC - zCO ) / ( nuC * molarMassC ) + ( zH - zHO ) / ( nuH * molarMassH ) 
                + 2.0 * ( zOO - zO ) / ( nuO * molarMassO ) )
            / ( ( zCF - zCO ) / ( nuC * molarMassC ) + ( zHF - zHO ) / ( nuH * molarMassH ) 
                    + 2.0 * ( zOO - zOF ) / ( nuO * molarMassO ) );
    }

    return z;
}

    template<class Species>
Double T1DFlame<Species>::GetElementMassFraction( Double *Y, const char *const atomName, Double atomMolarMass )
{
    TInputDataPtr	inp = GetInputData();
    int				i, atomNumber = inp->FindAtomIndex( atomName );
    int 			nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    //SpeciesPtr		species = inp->GetSpecies();		// ### composition can be determined by fSpecies.GetComposition()
    Double			z = 0.0;

    if ( atomNumber > -1 ) {
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            z += fSpecies.GetComposition()->mat[i][atomNumber] * atomMolarMass / fSpecies.GetMolarMass()->vec[i] * Y[i];		// ###
        }
    }

    return z;
}

    template<class Species>
void T1DFlame<Species>::XToEta( TNewtonPtr bt, VectorPtr etaVec )
{
    int			k;
    int			gridPoints = bt->GetCurrentGridPoints();
    TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
    Double		factor = 0.5 * sqrt( GetStrainRate() / ( fFlameNode->rhoInf * fFlameNode->viscosityInf ) );
    Double		left = bt->GetLeft();
    Double		right = bt->GetRight();
    Double		*rho = GetProperties()->GetDensity()->vec;
    Double		*x = grid->GetX()->vec;
    Double		*eta = etaVec->vec;

    if ( etaVec->len < gridPoints + 2 ) {
        fprintf( stderr, "%s\n", "#warning: Vector etaVec too short, values for physical grid are not computed"  );
        return;
    }

    eta[0] = 0.0;
    eta[1] = eta[0] + factor * ( x[0] - left ) * ( rho[-1] + rho[0] );
    for ( k = 1; k < gridPoints; ++k ) {
        eta[k+1] = eta[k] + factor * ( x[k] - x[k-1] ) * ( rho[k-1] + rho[k] );
    }
    eta[gridPoints+1] = eta[gridPoints] + factor * ( right - x[gridPoints-1] ) * ( rho[gridPoints-1] + rho[gridPoints] );
}

    template<class Species>
void T1DFlame<Species>::OriginToZstoich( VectorPtr xVec, VectorPtr mixtureFraction, Double zStoich )
{
    // here xVec contains boundary points

    int			k;
    int			gridPoints = xVec->len - 2;
    int			kBefStoech = -1;
    Double		*Z = mixtureFraction->vec;
    Double		*x = xVec->vec;
    Double		xStoich;

    for ( k = 0; k < gridPoints-1; ++k ) {
        if ( ( Z[k] - zStoich ) * ( Z[k+1] - zStoich ) <= 0.0 ) {
            kBefStoech = k;
            break;
        }
    }

    if ( kBefStoech < 0 ) {
        fprintf( stderr, "%s%g\n", "##warning: can't find the point of stoichiometric mixture fraction, Zst = ", zStoich  );
        return;
    }

    xStoich = x[kBefStoech+1] + ( x[kBefStoech+2] - x[kBefStoech+1] ) 
        / ( Z[kBefStoech+1] - Z[kBefStoech] )
        * ( zStoich - Z[kBefStoech] );

    for ( k = 0; k < gridPoints+2; ++k ) {
        x[k] -= xStoich;
    }
}

    template<class Species>
Double T1DFlame<Species>::ComputeEmissionIndex( int speciesIndex, Double *x )
{
    int			k;
    TNewtonPtr	bt = fSolver->bt;
    TGridPtr 	grid = bt->GetGrid()->GetFine();
    int			nGridPoints = grid->GetNGridPoints();
    char		**names = fSpecies.GetNames();
    Double		**m = fSpecies.GetProductionRate()->mat;
    Double		sum;

    sum = m[0][speciesIndex] * ( x[1] - bt->GetLeft() );  // m[0] = 0, m[L] = 0
    for ( k = 1; k < nGridPoints-1; ++k ) {
        sum += m[k][speciesIndex] * ( x[k+1] - x[k-1] );
    }
    sum += m[nGridPoints-1][speciesIndex] * ( bt->GetRight() - nGridPoints-1 );  // m[0] = 0, m[L] = 0

    return 0.5 * sum;
}

    template<class Species>
Double T1DFlame<Species>::ThermoDiffusion( int speciesIndex, CoordType coordinate, NodeInfoPtr nodeInfo )
{
    // fills the jacobian with 
    //		d/dy( rho D_iT dln(T)/dy ) for similarity coordinate
    // and with
    //		d/dy( D_iT dln(T)/dy ) for physical coordinate

    int		tempOff = GetOffsetTemperature();
    Double	yPrev = nodeInfo->yPrev[tempOff];
    Double	y =  nodeInfo->y[tempOff];
    Double	yNext =  nodeInfo->yNext[tempOff];
    Double	**thermDiff = fFlameNode->diffTherm;
    Double	diffPlusHm, diffMinusH;

    if ( coordinate == kPhysical ) {
        diffPlusHm = nodeInfo->hm * ( thermDiff[kCurr][speciesIndex] / nodeInfo->y[tempOff]
                + thermDiff[kNext][speciesIndex] / nodeInfo->yNext[tempOff] );
        diffMinusH = nodeInfo->h * ( thermDiff[kCurr][speciesIndex] / nodeInfo->y[tempOff]
                + thermDiff[kPrev][speciesIndex] / nodeInfo->yPrev[tempOff] );
    }
    else {
        Double	*rho = fFlameNode->mixDensity;
        diffPlusHm = nodeInfo->hm * ( rho[kCurr] * thermDiff[kCurr][speciesIndex] / nodeInfo->y[tempOff]
                + rho[kNext] * thermDiff[kNext][speciesIndex] / nodeInfo->yNext[tempOff] );
        diffMinusH = nodeInfo->h * ( rho[kCurr] * thermDiff[kCurr][speciesIndex] / nodeInfo->y[tempOff]
                + rho[kPrev] * thermDiff[kPrev][speciesIndex] / nodeInfo->yPrev[tempOff] );
    }

    return ( ( diffPlusHm * ( yNext - y ) + diffMinusH * ( yPrev - y ) ) 
            / nodeInfo->hnenn );
}

    template<class Species>
void T1DFlame<Species>::FillJacThermoDiffusion( int nVariable, Double constCoeff, CoordType coordinate, NodeInfoPtr nodeInfo )
{
    // fills the jacobian with 
    //		constCoeff * d/dy( rho D_iT dln(T)/dy ) for similarity coordinate
    // and with
    //		constCoeff * d/dy( D_iT dln(T)/dy ) for physical coordinate

    int		tempOff = GetOffsetTemperature();
    int		specInd = nVariable - GetOffsetFirstSpecies();
    Double	*temp = fFlameNode->temp;
    Double	**thermDiff = fFlameNode->diffTherm;
    Double	diffPlusHm, diffMinusH;

    if ( coordinate == kPhysical ) {
        diffPlusHm = constCoeff * nodeInfo->hm * ( thermDiff[kCurr][specInd]
                + thermDiff[kNext][specInd] );
        diffMinusH = constCoeff * nodeInfo->h * ( thermDiff[kCurr][specInd]
                + thermDiff[kPrev][specInd] );
    }
    else {
        Double	*rho = fFlameNode->mixDensity;
        diffPlusHm = constCoeff * nodeInfo->hm * ( rho[kCurr] * thermDiff[kCurr][specInd]
                + rho[kNext] * thermDiff[kNext][specInd] );
        diffMinusH = constCoeff * nodeInfo->h * ( rho[kCurr] * thermDiff[kCurr][specInd]
                + rho[kPrev] * thermDiff[kPrev][specInd] );
    }

    nodeInfo->a[tempOff][nVariable] -= ( diffPlusHm + diffMinusH ) / temp[kCurr];
    if ( !nodeInfo->lastPoint ) {
        nodeInfo->b[tempOff][nVariable] += diffPlusHm / temp[kNext];
    }
    if ( !nodeInfo->firstPoint ) {
        nodeInfo->c[tempOff][nVariable] += diffMinusH / temp[kPrev];
    }
}

    template<class Species>
int T1DFlame<Species>::CheckFlameLocation( Double xMid, Double deltaXVar )
{
    //	returns number of gridpoints to move, where positive nPoints 
    //	means shift from left to right

    int			nPoints = 0;
    TNewtonPtr	bt = GetSolver()->bt;
    TGridPtr 	grid = bt->GetGrid()->GetFine();
    Double		*x = grid->GetX()->vec;
    int			nGridPoints = grid->GetNGridPoints();
    int			tempLoc = LocationOfMax( nGridPoints+2, &fSolTemp->vec[kPrev] ) - 1;
    Double		deltaX =  bt->GetRight() - bt->GetLeft();
    Double		lowerBound = ( xMid - 0.5 * deltaXVar ) * deltaX;
    Double		upperBound = ( xMid + 0.5 * deltaXVar ) * deltaX;

    if ( x[tempLoc] < lowerBound ) {
        Double	newRight = bt->GetRight() - 0.5 * deltaXVar * deltaX;
        nPoints = GetNOfX( newRight, nGridPoints, x ) - nGridPoints;
    }
    else if ( x[tempLoc] > upperBound ) {
        Double	newLeft = bt->GetLeft() + 0.5 * deltaXVar * deltaX;
        nPoints = GetNOfX( newLeft, nGridPoints, x );
    }

    return nPoints;
}

    template<class Species>
int T1DFlame<Species>::GetNOfX( Double theX, int nGridPoints, Double *x )
{
    //	returns the number of the smallest x greater than theX
    int	loc = 0;

    while ( loc < nGridPoints && x[loc] < theX ) ++loc;

    return loc;
}


    template<class Species>
Flag T1DFlame<Species>::AdjustGrid( PostIterFuncPtr PostIter )
{
    //	adjust flame location
    //	first compute number of points to move
    Double		xMid = 0.5;
    Double		deltaXVar = 0.25;
    Double		xToMove;
    int 		nPointsToMove = CheckFlameLocation( xMid, deltaXVar );

    if ( nPointsToMove ) {
        NodeMover::fromType	fromTo = NodeMover::fromLeftSide;
        if ( nPointsToMove < 0 ) {
            nPointsToMove = -nPointsToMove;
            fromTo = NodeMover::fromRightSide;
        }
        fprintf( stderr, "%s%d%s%s%s\n", "move grid by ", nPointsToMove, " points to the " 
                , ( ( fromTo == NodeMover::fromLeftSide ) ? "right" : "left" ) 
                , NEWL  );
        int			i, k;
        int			fVVelocity = GetOffsetVVelocity();
        TNewtonPtr 	bt = GetSolver()->bt;
        TGridPtr 	grid = bt->GetGrid()->GetFine();
        int			nOfVars = bt->GetNVariables();
        int			nGridPoints = grid->GetNGridPoints();
        Double		**y = grid->GetY()->mat;
        Double		*yLeft = grid->GetYLeft()->vec;
        Double		*yRight = grid->GetYRight()->vec;
        VectorPtr	xVec = grid->GetX();
        Double		*x = xVec->vec;
        MMDataBag	bag( nOfVars );

        //	init
        bag.Initialize();
        bag.SetOldInpedVar( xVec, "x" );
        for ( i = 0; i < nOfVars; ++i ) {
            bag.Insert( &y[0][i], nGridPoints, nOfVars, GetVariableNames()[i] );
        }
        bag[fVVelocity].SetXType( MMDataSet::lastGradient );

        //	set new grid
        VectorPtr	newXVec = NewVector( nGridPoints );
        NodeMover nm( xVec, newXVec, nPointsToMove, fromTo );
        nm.MoveIt();
        bag.SetNewInpedVar( newXVec, "xNew" );

#undef DEBUGMAPPING
#ifdef DEBUGMAPPING
        ofstream os( GetOutfileName( "DataBag", FileType::kText ), ios::out );
        os, bag;
#endif			

        //	map
        VectorPtr	newYVec = NewVector( nGridPoints );
        Double		*newY = newYVec->vec;
        for ( i = 0; i < bag.NumElems(); ++i ) {
            bag[i].Map( newYVec );
            for ( k = 0; k < nGridPoints; ++k ) {
                y[k][i] = newY[k];
            }
        }

        //	shift new grid
        Double	shift;
        Double	right;
        Double	*xNew = newXVec->vec;
        if ( fromTo == NodeMover::fromLeftSide ) {
            shift = x[nPointsToMove] - x[0];
            right = bt->GetRight() + nPointsToMove * ( x[nGridPoints-1] - x[nGridPoints-2] ) - shift;
        }
        else {
            shift = - nPointsToMove * ( x[1] - x[0] );
            right = x[nGridPoints-nPointsToMove] - shift;
        }
        for ( k = 0; k < nGridPoints; ++k ) {
            x[k] = xNew[k] - shift;
        }
        bt->SetRight( right );

        //	set bondary values for VVelocity
        int lp = nGridPoints-1;
        for ( i = 0; i < nOfVars; ++i ) {
            if ( bag[i].GetXType() == MMDataSet::lastGradient ) {
                yLeft[i] = y[0][i] - ( y[1][i] - y[0][i] ) / ( x[1] - x[0] ) 
                    * ( x[0] - bt->GetLeft() );
                yRight[i] = y[lp][i] + ( y[lp][i] - y[lp-1][i] ) / ( x[lp] - x[lp-1] ) 
                    * ( bt->GetRight() - x[lp] );
            }
        }

        PostIter( this );

        //	clean up
        DisposeVector( newYVec );
        DisposeVector( newXVec );
        fSolver->ReInit();

        return TRUE;
    }
    else if ( (xToMove = CheckGradientLeft()) != 0 ) {
        int			k;
        TNewtonPtr 	bt = GetSolver()->bt;
        TGridPtr 	grid = bt->GetGrid()->GetFine();
        Double		*x = grid->GetX()->vec;
        int			nGridPoints = grid->GetNGridPoints();

        fprintf( stderr, "%s%g%s\n", "enlarge left side of grid by ", xToMove, NEWL  );
        bt->SetRight( bt->GetRight() - xToMove );
        for ( k = 0; k < nGridPoints; ++k ) {
            x[k] -= xToMove;
        }

        //	set bondary values for VVelocity
        int 		lp = nGridPoints-1;
        int			vOff = GetOffsetVVelocity();
        Double		*yLeft = grid->GetYLeft()->vec;
        Double		*yRight = grid->GetYRight()->vec;
        Double		**y = grid->GetY()->mat;

        yLeft[vOff] = y[0][vOff] - ( y[1][vOff] - y[0][vOff] ) / ( x[1] - x[0] ) 
            * ( x[0] - bt->GetLeft() );
        yRight[vOff] = y[lp][vOff] + ( y[lp][vOff] - y[lp-1][vOff] ) / ( x[lp] - x[lp-1] ) 
            * ( bt->GetRight() - x[lp] );

        PostIter( this );

        fSolver->ReInit();
        return TRUE;
    } 
    else if ( (xToMove = CheckGradientRight()) != 0 ) {
        TBVPSolverPtr	solver = GetSolver();
        TNewtonPtr		bt = solver->bt;
        TGridPtr 	grid = bt->GetGrid()->GetFine();
        Double		*x = grid->GetX()->vec;
        int			nGridPoints = grid->GetNGridPoints();

        fprintf( stderr, "%s%s%g%s\n", ( ( xToMove > 0.0 ) ? "enlarge" : "cut" ) 
                , " right side of grid by ", xToMove, NEWL  );
        if ( xToMove < 0.0 ) {	// cut right
            nPointsToMove = (int) xToMove;
            bt->SetRight( x[nGridPoints+nPointsToMove] );
            grid->AdjustNGridPoints( nGridPoints+nPointsToMove );
            solver->UpdateAllDimensions( nGridPoints+nPointsToMove );
        }
        else {					// enlarge right
            bt->SetRight( bt->GetRight() + xToMove );
        }

        //	set bondary values for VVelocity
        int 		lp = grid->GetNGridPoints()-1;
        int			vOff = GetOffsetVVelocity();
        Double		*yLeft = grid->GetYLeft()->vec;
        Double		*yRight = grid->GetYRight()->vec;
        Double		**y = grid->GetY()->mat;

        yLeft[vOff] = y[0][vOff] - ( y[1][vOff] - y[0][vOff] ) / ( x[1] - x[0] ) 
            * ( x[0] - bt->GetLeft() );
        yRight[vOff] = y[lp][vOff] + ( y[lp][vOff] - y[lp-1][vOff] ) / ( x[lp] - x[lp-1] ) 
            * ( bt->GetRight() - x[lp] );

        PostIter( this );

        fSolver->ReInit();
        return TRUE;
    } 
    else {
        return FALSE;
    }
}

    template<class Species>
Double T1DFlame<Species>::CheckGradientLeft( void )
{
    TNewtonPtr	bt = GetSolver()->bt;
    TGridPtr 	grid = bt->GetGrid()->GetFine();
    Double		*x = grid->GetX()->vec;
    Double		*temp = fSolTemp->vec;
    int			nGridPoints = grid->GetNGridPoints();
    Double		left =  bt->GetLeft();
    Double		right =  bt->GetRight();
    Double		L =  right - left;
    Double		deltaX;
    Double		tempBound;

    // check left
    deltaX = x[0] - left;
    tempBound = 20.0 * deltaX / L;
    if ( fabs( temp[0] - temp[-1] ) > tempBound ) {
        return -2.0 * deltaX;
    }

    return 0.0;
}

    template<class Species>
Double T1DFlame<Species>::CheckGradientRight( void )
{
    TNewtonPtr	bt = GetSolver()->bt;
    TGridPtr 	grid = bt->GetGrid()->GetFine();
    Double		*x = grid->GetX()->vec;
    Double		*temp = fSolTemp->vec;
    int			nGridPoints = grid->GetNGridPoints();
    Double		left =  bt->GetLeft();
    Double		right =  bt->GetRight();
    Double		L =  right - left;
    Double		deltaX;
    Double		tempBound;

    // check right
    deltaX = right - x[nGridPoints-1];
    // upper
    tempBound = 20.0 * deltaX / L;
    if ( fabs( temp[nGridPoints] - temp[nGridPoints-1] ) > tempBound ) {
        return 2.0 * deltaX;
    }
    // lower
    tempBound = 2.0 * deltaX / L;
    if ( fabs( temp[nGridPoints] - temp[nGridPoints-1] ) < tempBound ) {
        return -2.0;
    }

    return 0.0;
}

    template<class Species>
int T1DFlame<Species>::SensitivityAnalysis( Double coeffSpecies, Double coeffTemp, CoordType coordinate )
{ 
    //cai: TUnstrPremFlamePhys.C:	flame->SensitivityAnalysis( 1.0, 1.0, kPhysical );

    if ( fNSensObj == 0 ) {
        return -1;
    }


    int 				i;				// variables
    int 				j;				// reactions
    int 				k;				// coordinate
    int					m;				// sensitivity objects
    int					n;				// species
    int					speciesOffset = GetOffsetFirstSpecies();
    int					tempOffset = GetOffsetTemperature();
    int					nSpeciesIn = fSpecies.GetNSpeciesInSystem();
    int					nReactions = fReaction.GetNOfReactions();
    VectorPtr			APtr = fReaction.GetA();
    Double				*A = APtr->vec;
    VectorPtr			*nu = fReaction.GetNu();
    Double				*nuvec;   //cai: v stochiometric coefficient
    IntVectorPtr		*speciesNumber = fReaction.GetSpeciesNumber();
    int					*specNumvec;
    Double				*molarMass = fSpecies.GetMolarMass()->vec;
    TNewtonPtr			bt = fSolver->bt;
    int					nGridPoints = bt->GetCurrentGridPoints();
    NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    MatrixPtr			dy = bt->GetDy();
    MatrixPtr			y = bt->GetGrid()->GetFine()->GetY();
    IntVectorPtr		objectsPtr = NewIntVector( fNSensObj );	// index pointer to sensitivity objects in solution vector
    int					*objects = objectsPtr->vec;
    int					&objectslen = objectsPtr->len;
    VectorPtr			maxObjPtr = NewVector( fNSensObj );
    Double				*maxObj = maxObjPtr->vec;
    int					&maxObjlen = maxObjPtr->len;
    ConstStringArray variableNames = GetVariableNames();
    int					nVariables = GetVariablesWithoutSpecies() + nSpeciesIn;
    Double				coeffCoord;
    char				name[128];
    TensorPtr			sensitivityCoefficients;
    Double				hnenn;
    Double				rateOverCoeff;


    objectslen = 0;//cai: int	&objectslen = objectsPtr->len;  objectsPtr->len is the adress of objetslen (int) 
    maxObjlen = 0;
    for ( m = 0; m < fNSensObj; ++m ) {
        for ( i = 0; i < nVariables; ++i ) {

            // search object cai:try to find species in nVariables
            strcpy( name, variableNames[i] );
            UpperString( name );
            if ( strcmp( fSensObj[m], name ) == 0 ) {
                objects[objectslen++] = i;

                // search maximum value and save it    cai:try to find the maximal value of  Y_speices    maxobjlen++ see below
                maxObj[maxObjlen] = 0.0;
                for ( k = 0; k < nGridPoints; ++k ) {
                    if ( y->mat[k][i] > maxObj[maxObjlen] ) maxObj[maxObjlen] = y->mat[k][i];
                }
                if ( maxObj[maxObjlen] == 0.0 ) FatalError( " # maxObj less or equal 0.0" );
                ++maxObjlen;
                break;
            }
        }
    }
    if ( objectslen == 0 ) return -1;

    fprintf( stderr, "%s", "sensitivity analysis .. " );

    sensitivityCoefficients = NewTensor( objectslen, nGridPoints, nReactions, kRowPointers );

    for ( j = 0; j < nReactions; ++j ) {
        specNumvec = speciesNumber[j]->vec;
        nuvec = nu[j]->vec;
        ClearMatrix( dy );
        for ( k = 0; k < nGridPoints; ++k ) {
            bt->SetNodeInfo( this, k );
            SetFlameNode( k );
            rateOverCoeff = fFlameNode->reactionRate[j] / A[j]; 
            hnenn = bt->GetNodeInfo()->hnenn;
            switch ( coordinate ) {
                case kPhysical:
                    coeffCoord = 1.0;
                    break;
                case kSimilarity:
                    coeffCoord = - 1.0 / (*fFlameNode->mixDensity);
                    break;
                default:
                    FatalError( "# unknown coordinate type in sensitivity analysis" );
                    break;
            }

            //	fill rhs for species
            for ( i = 0; i < speciesNumber[j]->len; ++i )
                dy->mat[k][specNumvec[i] + speciesOffset] =
                    coeffCoord * coeffSpecies * molarMass[specNumvec[i]] * (-nuvec[i])
                    * rateOverCoeff * hnenn;     //cai: coeffSpecies from funktion input

            // fill rhs for temperature
            for ( i = 0; i < speciesNumber[j]->len; ++i ) {
                dy->mat[k][tempOffset] -= fFlameNode->enthalpy[specNumvec[i]]
                    * molarMass[specNumvec[i]] * (-nuvec[i]);
            }
            dy->mat[k][tempOffset] *= coeffCoord * coeffTemp * rateOverCoeff / (*fFlameNode->mixHeatCapacity)
                * hnenn;
        }

        // backsolve
        bt->BackSolve( dy, NULL, FALSE );

        // save sensitivity coefficients
        for ( k = 0; k < nGridPoints; ++k ) {
            bt->SetNodeInfo( this, k );
            SetFlameNode( k );
            for ( m = 0; m < objectslen; ++m ) {
                sensitivityCoefficients->tensor[m][k][j] = A[j] / maxObj[m] * dy->mat[k][objects[m]];
                if ( fSpecies.FindSpecies( fSensObj[m] ) >= 0 ) {
                    for ( n = 0; n < nSpeciesIn; ++n ) sensitivityCoefficients->tensor[m][k][j] -=
                        A[j] * (*fFlameNode->mixMolarMass) /
                            molarMass[n] * dy->mat[k][n+speciesOffset];
                }
            }
        }
    }

    // output
    PrintSensitivityCoefficients( sensitivityCoefficients, objectsPtr );
    PrintSensMax( sensitivityCoefficients, objectsPtr );

    DisposeTensor( sensitivityCoefficients );
    DisposeVector( maxObjPtr );	
    DisposeIntVector( objectsPtr );	
    fprintf( stderr, "%s\n", "done." );
    return 0;
}

    template<class Species>
void T1DFlame<Species>::PrintSensMax( TensorPtr sensitivityCoefficients, IntVectorPtr objectsPtr )
{
    ConstStringArray	variableNames = GetVariableNames();
    ConstStringArray	reactionLabels = fReaction.GetLabels();
    char			fName[128];
    TNewtonPtr		bt = fSolver->bt;
    TGridPtr		fine = bt->GetGrid()->GetFine();
    int				fTemperature = GetOffsetTemperature();
    int				gridPoints = fine->GetX()->len;
    Double			*x = fine->GetX()->vec;
    int				*objects = objectsPtr->vec;
    int				&objectslen = objectsPtr->len;
    int				nReactions = fReaction.GetNOfReactions();
    int				i, j;

    for ( i = 0; i < objectslen; ++i ) {
        sprintf( fName, "Senm%.8s", variableNames[objects[i]] );  //right file for sensitivity analysis
        for ( j = 0; j < strlen(fName); ++j ) {
            if ( fName[j] == '/' ) {
                fName[j] = '_';
            }
        }
        PrintOneSensMax( sensitivityCoefficients->tensor[i], gridPoints, nReactions
                , fName );
    }
}

    template<class Species>
void T1DFlame<Species>::PrintOneSensMax( Double **sensCoeff, int gridPoints, int nReactions
        , const char *fileName )
{
    int		i;
    FILE	*fp = GetOutputFile( fileName, "", FileType::kText );
    Double	*max = new Double[nReactions];
    char	**label = fReaction.GetLabels();
    int		maxReaction;

    for ( i = 0; i < nReactions; ++i ) {
        max[i] = GetMaxSens( sensCoeff, gridPoints, i );
    }

    fprintf( fp, "maximum values of the sensitivity coefficients\n\n" );
    for ( i = 0; i < nReactions; ++i ) {
        maxReaction = LocationOfAbsMax( nReactions, max );
        fprintf( fp, "%-5s:\t%12.6g\t", label[maxReaction], max[maxReaction] );
        fReaction.PrintReactionEquation( maxReaction, &fSpecies, fp );
        fprintf( fp, "\n" );
        max[maxReaction] = 0.0;
    }

    delete[] max;
    fclose( fp );
}

    template<class Species>
Double T1DFlame<Species>::GetMaxSens( Double **sensCoeff, int gridPoints, int reaction )
{
    int		loc = 0;
    Double	high = sensCoeff[0][reaction];

    for ( int k = 1; k < gridPoints; ++k ) {
        if ( fabs( sensCoeff[k][reaction] ) > high ) {
            loc = k;
            high = fabs( sensCoeff[k][reaction] );
        }
    }

    return sensCoeff[loc][reaction];
}

    template<class Species>
void T1DFlame<Species>::PrintSensitivityCoefficients( TensorPtr sensitivityCoefficients, IntVectorPtr objectsPtr )
{
    ConstStringArray	variableNames = GetVariableNames();
    ConstStringArray	reactionLabels = fReaction.GetLabels();
    TNewtonPtr		bt = fSolver->bt;
    TGridPtr		fine = bt->GetGrid()->GetFine();
    int				fTemperature = GetOffsetTemperature();
    int				gridPoints = fine->GetX()->len;
    Double			*x = fine->GetX()->vec;
    int				*objects = objectsPtr->vec;
    int				&objectslen = objectsPtr->len;
    int				nReactions = fReaction.GetNOfReactions();
    char			name[128];
    unsigned int    i, j;

    for ( i = 0; i < objectslen; ++i ) {
        strcpy( name, variableNames[objects[i]] );
        for ( j = 0; j < strlen(name); ++j ) {
            if ( name[j] == '/' ) {
                name[j] = '_';
            }
        }
        sprintf( GetOutFileBuff(), "%sSens%.8s_p%.2dT%.4d"
                , GetOutputPath()
                , variableNames[objects[i]]
                , ( int )( GetPressure() * 1.0e-5 )
                , ( int )( fine->GetYLeft()->vec[fTemperature] ) );

        SaveArray( sensitivityCoefficients->tensor[i], gridPoints, nReactions,
                kRowPointers, x, reactionLabels, GetOutFileBuff() );
    }
}

    template<class Species>
void T1DFlame<Species>::ReactionFluxes( CoordType coordinate )
{
    int 				j;				// reactions
    int 				k;				// coordinate
    int					i;				// species
    int					nSpecies = fSpecies.GetNOfSpecies();
    int					nReactions = fReaction.GetNOfReactions();
    TNewtonPtr			bt = fSolver->bt;
    int					nGridPoints = bt->GetCurrentGridPoints();
    MatrixPtr			fluxesPtr = NewMatrix( nSpecies, nReactions, kRowPointers);
    Double				**fluxes = fluxesPtr->mat;
    MatrixPtr			fractionsPtr = NewMatrix( nSpecies, nReactions, kRowPointers);
    Double				**fractions = fractionsPtr->mat;
    VectorPtr			sumFormPtr = NewVector( nSpecies );
    Double				*sumForm = sumFormPtr->vec;
    VectorPtr			sumConsPtr = NewVector( nSpecies );
    Double				*sumCons = sumConsPtr->vec;
    VectorPtr			*nuPtr = fReaction.GetNu();
    Double				*nu; 
    IntVectorPtr		*speciesNumberPtr = fReaction.GetSpeciesNumber();
    int					*specNum;
    Double				reactionRate;
    Double				*molarMass = fSpecies.GetMolarMass()->vec;
    Double				dummy;
    const Double		kg2g = 1.0e3;	
    VectorPtr			xPtr = NewVector( nGridPoints + 2 );
    Double				*x = xPtr->vec + 1;		// now pointing to the second element

    Double				epsilon = 5.0;


    fprintf( stderr, "%s", "compute reaction fluxes .. " );

    if ( coordinate == kSimilarity ) {
        EtaToX( bt, xPtr );
    }
    else {
        TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
        Double		*theX = grid->GetX()->vec;
        copy( nGridPoints, theX, 1, x, 1 );
        x[-1] = bt->GetLeft();
        x[nGridPoints] = bt->GetRight();
    }

    for ( k = 0; k < nGridPoints; ++k ) {
        SetFlameNode( k );
        for ( j = 0; j < nReactions; ++j ) {
            reactionRate = fFlameNode->reactionRate[j];
            specNum = speciesNumberPtr[j]->vec;
            nu = nuPtr[j]->vec;
            for ( i = 0; i < speciesNumberPtr[j]->len; ++i ) {
                dummy = (-nu[i]) * reactionRate * (x[k+1] - x[k-1]) / 2.0 * kg2g
                    * molarMass[specNum[i]];
                fluxes[specNum[i]][j] += dummy;
                if ( dummy > 0.0 ) {
                    sumForm[specNum[i]] += dummy;
                }
                else {
                    sumCons[specNum[i]] += dummy;
                }
            }
        }
    }

    for ( j = 0; j < nReactions; ++j ) {
        specNum = speciesNumberPtr[j]->vec;
        for ( i = 0; i < speciesNumberPtr[j]->len; ++i ) {
            if ( fluxes[specNum[i]][j] > 0.0 ) {
                fractions[specNum[i]][j] = fluxes[specNum[i]][j] / sumForm[specNum[i]] * 100.0;
            }
            else if ( fluxes[specNum[i]][j] < 0.0 ) {
                fractions[specNum[i]][j] = fluxes[specNum[i]][j] / (-sumCons[specNum[i]]) * 100.0;
            }
        }
    }

    // Output relative values
    PrintReactionFluxesReac( fluxes, "Flux.reac.abs", "integrated production rates [g/m^3s]" );
    PrintReactionFluxesReac( fractions, "Flux.reac.rel", "production rate of one reaction over production rate [%]" );
    PrintReduceInfo( fractions, epsilon, "reactions with no contribution larger than epsilon" );
    PrintReactionFluxesSpec( fluxes, "Flux.spec.abs", "integrated production rates [g/m^3s]" );
    PrintReactionFluxesSpec( fractions, "Flux.spec.rel", "production rate of one reaction over production rate [%]" );

    // do the same for added forward and backward reactions
    ClearVector( sumConsPtr );
    ClearVector( sumFormPtr );
    ClearMatrix( fluxesPtr );
    ClearMatrix( fractionsPtr );

    int	*backw = fReaction.GetBackwardReacs()->vec;

    for ( k = 0; k < nGridPoints; ++k ) {
        SetFlameNode( k );
        for ( j = 0; j < nReactions; ++j ) {
            if ( fReaction.IsBackwardReaction( j ) ) {
                continue;
            }
            if ( fReaction.IsForwardReaction( j ) ) {
                reactionRate = ( fFlameNode->reactionRate[j] 
                        - fFlameNode->reactionRate[backw[j]] );
            }
            else {
                reactionRate = fFlameNode->reactionRate[j];
            }
            specNum = speciesNumberPtr[j]->vec;
            nu = nuPtr[j]->vec;
            for ( i = 0; i < speciesNumberPtr[j]->len; ++i ) {
                dummy = (-nu[i]) * reactionRate * (x[k+1] - x[k-1]) / 2.0 * kg2g
                    * molarMass[specNum[i]];
                fluxes[specNum[i]][j] += dummy;
            }
        }
    }
    for ( j = 0; j < nReactions; ++j ) {
        specNum = speciesNumberPtr[j]->vec;
        for ( i = 0; i < speciesNumberPtr[j]->len; ++i ) {
            if ( fluxes[specNum[i]][j] > 0.0 ) {
                sumForm[specNum[i]] += fluxes[specNum[i]][j];
            }
            else {
                sumCons[specNum[i]] += fluxes[specNum[i]][j];
            }
        }
    }

    for ( j = 0; j < nReactions; ++j ) {
        specNum = speciesNumberPtr[j]->vec;
        for ( i = 0; i < speciesNumberPtr[j]->len; ++i ) {
            if ( fluxes[specNum[i]][j] > 0.0 ) {
                fractions[specNum[i]][j] = fluxes[specNum[i]][j] / sumForm[specNum[i]] * 100.0;
            }
            else if ( fluxes[specNum[i]][j] < 0.0 ) {
                fractions[specNum[i]][j] = fluxes[specNum[i]][j] / (-sumCons[specNum[i]]) * 100.0;
            }
        }
    }

    // Output relative values
    PrintReactionFluxesReac( fluxes, "Fluxa.reac.abs", "integrated production rates [g/m^3s]" );
    PrintReactionFluxesReac( fractions, "Fluxa.reac.rel", "production rate of one reaction over production rate [%]" );
    PrintReactionFluxesSpec( fluxes, "Fluxa.spec.abs", "integrated production rates [g/m^3s]" );
    PrintReactionFluxesSpec( fractions, "Fluxa.spec.rel", "production rate of one reaction over production rate [%]" );

    // CleanUp
    DisposeVector( xPtr );
    DisposeVector( sumConsPtr );
    DisposeVector( sumFormPtr );
    DisposeMatrix( fractionsPtr );
    DisposeMatrix( fluxesPtr );

    fprintf( stderr, "%s\n", "done." );
}

    template<class Species>
void T1DFlame<Species>::PrintReactionFluxesReac( Double **fluxes, const char *name, const char *header )
{
    FILE				*fp = NULL;
    ConstStringArray reactionLabels = fReaction.GetLabels();
    ConstStringArray names = fSpecies.GetNames();
    IntVectorPtr		*speciesNumberPtr = fReaction.GetSpeciesNumber();
    int					*specNum;
    int					nReactions = fReaction.GetNOfReactions();
    int					nSpecies = fSpecies.GetNOfSpecies();
    TNewtonPtr			bt = fSolver->bt;
    TGridPtr			fine = bt->GetGrid()->GetFine();
    int					fTemperature = GetOffsetTemperature();
    int					j;		// reactions
    int					i;		// species

    fp = GetOutputFile( name, "", FileType::kText );

    fprintf( fp, "%s\n\n", header );
    for ( j = 0; j < nReactions; ++j ) {
        fprintf( fp, "%-5s\t", reactionLabels[j] );
        specNum = speciesNumberPtr[j]->vec;
        for ( i = 0; i < speciesNumberPtr[j]->len; ++i ) {
            fprintf( fp, "%9.5f %-6s\t", fluxes[specNum[i]][j],
                    names[specNum[i]]); 
        }
        fputs( "\n", fp );
    }
    fclose( fp );
}

    template<class Species>
void T1DFlame<Species>::PrintReduceInfo( Double **fluxes, Double epsilon, const char *header )
{
    FILE				*fp = NULL;
    ConstStringArray reactionLabels = fReaction.GetLabels();
    ConstStringArray names = fSpecies.GetNames();
    IntVectorPtr		*speciesNumberPtr = fReaction.GetSpeciesNumber();
    int					*specNum;
    int					nReactions = fReaction.GetNOfReactions();
    int					nSpecies = fSpecies.GetNOfSpecies();
    TNewtonPtr			bt = fSolver->bt;
    TGridPtr			fine = bt->GetGrid()->GetFine();
    int					fTemperature = GetOffsetTemperature();
    int					j;		// reactions
    int					i;		// species
    Flag				eliminate;

    sprintf( GetOutFileBuff(), "%s%s_eps%gp%.2dT%.4d"
            , GetOutputPath()
            , "RedInfo"
            , epsilon
            , ( int )( GetPressure() * 1.0e-5 )
            , ( int )( fine->GetYLeft()->vec[fTemperature] ) );
    if ( !( fp = fopen( GetOutFileBuff(), "w") ) ) { 
        fprintf( stderr, "%s%s\n", "#warning: unable to open file ", GetOutFileBuff()  );
        exit(2);
    }

    fprintf( fp, "%s = %g\n\n", header, epsilon );
    for ( j = 0; j < nReactions; ++j ) {
        eliminate = TRUE;
        specNum = speciesNumberPtr[j]->vec;
        for ( i = 0; i < speciesNumberPtr[j]->len; ++i ) {
            if ( fabs( fluxes[specNum[i]][j] ) > epsilon ) {
                eliminate = FALSE;
                break;
            }
        }
        if ( eliminate ) {
            fprintf( fp, "%-5s\t", reactionLabels[j] );
            fReaction.PrintReactionEquation( j, &fSpecies, fp );
            fputs( "\n", fp );
        }
    }
    fclose( fp );
}

extern Double gSortBuffer[5000];
    template<class Species>
void T1DFlame<Species>::PrintReactionFluxesSpec( Double **fluxes, const char *name, const char *header )
{
    FILE			*fp = NULL;
    ConstStringArray	reactionLabels = fReaction.GetLabels();
    ConstStringArray names = fSpecies.GetNames();
    Double			*flux;
    int			nReactions = fReaction.GetNOfReactions();
    int			nSpecies = fSpecies.GetNOfSpecies();
    TNewtonPtr		bt = fSolver->bt;
    TGridPtr		fine = bt->GetGrid()->GetFine();
    int			fTemperature = GetOffsetTemperature();
    IntVectorPtr		indexPtr = NewIntVector( nReactions );
    int			*index = indexPtr->vec;
    int			*nOfUsedReactions = fSpecies.GetNOfUsedReactions()->vec;
    int			j;		// reactions
    int			i;		// species

    char reac[128];

    fp = GetOutputFile( name, "", FileType::kText );

    fprintf( fp, "%s\n\n", header );
    for ( i = 0; i < nSpecies; ++i ) {
        if ( nOfUsedReactions[i] ){
            flux = fluxes[i];
            fprintf( fp, "%s\n", names[i]);
            for ( j = 0; j < nReactions; ++j ) {
                index[j] = j;
                gSortBuffer[j] = flux[j];
            }
            qsort( index, nReactions, sizeof( int ), myCompare );		
            for ( j = 0; j < nReactions; ++j ) {
                if( flux[index[j]] != 0.0 ) {
                    fReaction.PrintReactionEquation (index[j], &fSpecies, reac);
                    fprintf (fp, "\t%12g\t%s:\t%s\n", flux[index[j]],
                            reactionLabels[index[j]], reac);
                }
            }
            fputs( "\n", fp );
        }
    }
    fclose( fp );

    DisposeIntVector( indexPtr );
}

    template<class Species>
void T1DFlame<Species>::ComputePropertiesMod( TFlameNodePtr flameNode, Double temp
        , Double *Y, Double pressure )
{
    int	nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    Double	*molarMass = fSpecies.GetMolarMass()->vec;
    Flag	newTemp;

    //  First compute molar mass of mixture
    fProperties->ComputeMixtureMolarMass( *flameNode->mixMolarMass, Y, molarMass, nSpeciesInSystem );

    //  compute properties of Species
    newTemp = fSpecies.ComputeSpeciesProperties( flameNode, temp, pressure );

    //	compute Delta_i, which is used by CompMixtureProps and Compute_DTherm
#ifdef OPTIMIZEDELTAI
    fSpecies.ComputeDeltaIOpt( flameNode, Y, flameNode->GijOverWj, newTemp );
#else
    fSpecies.ComputeDeltaI( flameNode->deltaI, Y, flameNode->viscosity );
#endif

    //  compute properties of the mixture
    fProperties->CompMixtureProps( flameNode, Y, temp, pressure, &fSpecies );

    if ( fSpecies.IsConstantLewisNumber() ) {
        fSpecies.Compute_D( flameNode );
    }
    else {
        fSpecies.Compute_D( flameNode, temp, Y, pressure, newTemp );
        if ( fThermoDiffusion ) {
            fSpecies.Compute_DTherm( flameNode, newTemp );
        }
    }

    //  compute properties of soot
    if ( fSoot ) {
        fSoot->ComputePolymereConcs( Y, temp, flameNode->mixDensity[kCurr]
                , molarMass, flameNode->Pij, flameNode->sumPi, flameNode->pahMoments
                , flameNode->moments, &fReaction );
        fSoot->ComputeDiffusivity( this );
    }
}

    template<class Species>
void T1DFlame<Species>::UpdateThermoProps( TFlameNodePtr flameNode, NodeInfoPtr nodeInfo )
{
    Double				pressure = GetPressure();
    Double				*molarMass = fSpecies.GetMolarMass()->vec;
    T1DRadiationPtr		radiation = fProperties->GetRadiation();
    Double				density;
    Double				*rateCoeff = flameNode->rateCoeff;
    Double				*currRateCoeff = flameNode->currRateCoeff;
    Flag				&kNewerThanW = flameNode->kNewerThanW[kCurr];
    Double				&tempReaction = flameNode->tempReaction[kCurr];
    Double				&pressureReaction = flameNode->pressureReaction[kCurr];
    Double				*YReaction = flameNode->YReaction;
    Double 				*reactionRate = flameNode->reactionRate;
    Double				*currReacRate = flameNode->currReacRate;
    Double				*tBConc = flameNode->tBodyConc;
    Double				temp = flameNode->temp[kCurr];
    Double				*Y = flameNode->Y[kCurr];

    CheckSolution( temp, Y, fSpecies.GetNSpeciesInSystem() );

    ComputePropertiesMod( flameNode, temp, Y, pressure );

    density = flameNode->mixDensity[kCurr];
#ifdef PRODRATEFILE
    Double	*concs = fReaction.GetMolarConcs()->vec;
    Double	*prodRate = flameNode->productionRate;
    fReaction.ComputeConcs( concs, Y, molarMass, density );
    fSpecies.ComputeTheProductionRates( prodRate, reactionRate
            , temp, pressure, concs, rateCoeff, tBConc );
    for ( int i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
        prodRate[i] *= molarMass[i];
    }
#else
    fReaction.CompThirdBodyConcs( tBConc, Y, molarMass, density );
    fReaction.ComputeRateCoefficients( rateCoeff, currRateCoeff, kNewerThanW, temp
            , tempReaction, pressure, pressureReaction, tBConc, &fSpecies );
    fReaction.ComputeReactionRates( reactionRate, kNewerThanW , currReacRate, rateCoeff
            , tBConc, density, Y, YReaction, molarMass, &fSpecies );
    fSpecies.ComputeProductionRates( flameNode->productionRate, reactionRate );
#endif
    ComputeDiffusivityCorrection( &(flameNode->Y[kCurr]), nodeInfo );

    double HeatRelease = 0.0;
    for ( int i = 0; i < fSpecies.GetNSpeciesInSystem(); ++i ) {
        //Double *ent = fSpecies.GetEnthalpy()->vec;

        HeatRelease += flameNode->productionRate[i]*flameNode->enthalpy[i];
    }

    HeatRelease = -HeatRelease;

    //cout << "Calculate the diffisvity" << endl;
    //cout << fFlameNode->diffusivity[kCurr] << " " << nodeInfo->x[kCurr] << endl;

    if ( radiation ) {
        if ( fRadiationName == "Adiabatic"){

        }
        else if (fRadiationName == string("Thin")){
            fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, temp
                    , Y, molarMass, density);
        }
        else if (fRadiationName == string("RadiativeFrac")){
            fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, temp
                    , Y, molarMass, density, HeatRelease, fRadiativeFrac, fFlameNode->kappa);
        }
        else if (fRadiationName == string("WSGG")){
            //cout <<fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint] << " " << nodeInfo->gridPoint << endl;
            //fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
            fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
        }
        else if (fRadiationName == string("WSGGJohansson")){
            //cout <<fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint] << " " << nodeInfo->gridPoint << endl;
            //fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
            fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
        }
        else if (fRadiationName == string("WSGGBordbar")){
            //cout <<fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint] << " " << nodeInfo->gridPoint << endl;
            //fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
            fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
        }
        else if (fRadiationName == string("Grey")){
            fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
        }
        else if (fRadiationName == string("SNB")){
            fProperties->GetRadiation()->ComputeRadiationOnePoint( fFlameNode->radiation, fProperties->GetRadiation()->dqrF[nodeInfo->gridPoint]);
        }
        else{
            cout << "#2 Error Radiation Model " <<  fRadiationName << " not found " << endl;
        }
    }
    if ( fSoot ) {
        fSoot->UpdateProductionRates( &fSpecies, &fReaction, flameNode->productionRate, density
                , Y, temp, flameNode->sumPi, flameNode->moments, flameNode->pahReactionRate
                , flameNode->mixMolarMass[kCurr]);
        fSoot->FillSource( flameNode->sootSource, this );
    }
}

    template<class Species>
TFlameNodePtr T1DFlame<Species>::NewTFlameNodeSaved( void )
{
    TFlameNodePtr	flameNodeSaved = new TFlameNode();
    const int		nThirdBodies = fReaction.GetNOfThirdBodies();
    const int		nReactions = fReaction.GetNOfReactions();
    const int		nSpecies = fSpecies.GetNOfSpecies();
    const int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();

    // reaction
    flameNodeSaved->tBodyConc = New1DArray( nThirdBodies );
    flameNodeSaved->rateCoeff = New1DArray( nReactions );
    flameNodeSaved->tempReaction = New1DArray( 1 );
    flameNodeSaved->pressureReaction = New1DArray( 1 );
    flameNodeSaved->currRateCoeff = New1DArray( nReactions );
    flameNodeSaved->kNewerThanW = new Flag[1];
    flameNodeSaved->reactionRate = New1DArray( nReactions );
    flameNodeSaved->YReaction = New1DArray( nSpeciesInSystem );
    flameNodeSaved->currReacRate = New1DArray( nReactions );
    flameNodeSaved->dMdY = New2DArray( nSpeciesInSystem, (nSpeciesInSystem+1) );

    // species
    flameNodeSaved->viscosity = New1DArray( nSpecies );
    flameNodeSaved->heatCapacity = New1DArray( nSpecies );
    flameNodeSaved->conductivity = New1DArray( nSpecies );
    flameNodeSaved->enthalpy = New1DArray( nSpecies );
    flameNodeSaved->diffusivity = New1DArray( nSpecies );
    flameNodeSaved->diffTherm = New2DArray( 1, nSpecies );
    flameNodeSaved->productionRate = New1DArray( nSpecies );
    flameNodeSaved->tempProp = New1DArray( 1 );
    flameNodeSaved->pressureProp = New1DArray( 1 );

    // properties
    flameNodeSaved->mixViscosity = New1DArray( 1 );
    flameNodeSaved->mixDensity = New1DArray( 1 );
    flameNodeSaved->mixConductivity = New1DArray( 1 );
    flameNodeSaved->mixHeatCapacity = New1DArray( 1 );
    flameNodeSaved->mixMolarMass = New1DArray( 1 );

    // flame
    flameNodeSaved->Y = New2DArray( 1, nSpecies );
    flameNodeSaved->temp = New1DArray( 1 );
    if ( fProperties->GetRadiation() ) {
        flameNodeSaved->radiation = New1DArray( 1 );
        flameNodeSaved->kappa = New1DArray( 1 );
    }
    flameNodeSaved->diffCorr = New1DArray( 1 );

    return flameNodeSaved;
}

    template<class Species>
void T1DFlame<Species>::DisposeTFlameNodeSaved( TFlameNodePtr flameNodeSaved )
{
    // flame
    Free1DArray( flameNodeSaved->diffCorr );
    if ( fProperties->GetRadiation() ) {
        Free1DArray( flameNodeSaved->radiation );
        Free1DArray( flameNodeSaved->kappa );
    }
    Free1DArray( flameNodeSaved->temp );
    Free2DArray( flameNodeSaved->Y );

    // properties
    Free1DArray( flameNodeSaved->mixMolarMass );
    Free1DArray( flameNodeSaved->mixHeatCapacity );
    Free1DArray( flameNodeSaved->mixConductivity );
    Free1DArray( flameNodeSaved->mixDensity );
    Free1DArray( flameNodeSaved->mixViscosity );

    // species
    Free1DArray( flameNodeSaved->pressureProp );
    Free1DArray( flameNodeSaved->tempProp );
    Free1DArray( flameNodeSaved->productionRate );
    Free2DArray( flameNodeSaved->diffTherm );
    Free1DArray( flameNodeSaved->diffusivity );
    Free1DArray( flameNodeSaved->enthalpy );
    Free1DArray( flameNodeSaved->conductivity );
    Free1DArray( flameNodeSaved->heatCapacity );
    Free1DArray( flameNodeSaved->viscosity );

    // reaction
    Free2DArray( flameNodeSaved->dMdY );
    Free1DArray( flameNodeSaved->currReacRate );
    Free1DArray( flameNodeSaved->YReaction );
    Free1DArray( flameNodeSaved->reactionRate );
    delete flameNodeSaved->kNewerThanW;
    Free1DArray( flameNodeSaved->currRateCoeff );
    Free1DArray( flameNodeSaved->tempReaction );
    Free1DArray( flameNodeSaved->pressureReaction );
    Free1DArray( flameNodeSaved->rateCoeff );
    Free1DArray( flameNodeSaved->tBodyConc );


    delete flameNodeSaved;
}

    template<class Species>
void T1DFlame<Species>::CopyFlameNode( TFlameNodePtr flameNodeSource, TFlameNodePtr flameNodeDest )
{
    int				nThirdBodies = fReaction.GetNOfThirdBodies();
    int				nReactions = fReaction.GetNOfReactions();
    int				nSpecies = fSpecies.GetNOfSpecies();
    int				nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();

    // reaction
    memcpy( flameNodeDest->tBodyConc, flameNodeSource->tBodyConc, sizeof( Double ) * nThirdBodies ); 
    memcpy( flameNodeDest->rateCoeff, flameNodeSource->rateCoeff, sizeof( Double ) * nReactions ); 
    memcpy( flameNodeDest->tempReaction, flameNodeSource->tempReaction, sizeof( Double ) ); 
    memcpy( flameNodeDest->pressureReaction, flameNodeSource->pressureReaction, sizeof( Double ) ); 
    memcpy( flameNodeDest->currRateCoeff, flameNodeSource->currRateCoeff, sizeof( Double ) * nReactions ); 
    memcpy( flameNodeDest->kNewerThanW, flameNodeSource->kNewerThanW, sizeof( Flag ) ); 
    memcpy( flameNodeDest->reactionRate, flameNodeSource->reactionRate, sizeof( Double ) * nReactions ); 
    memcpy( flameNodeDest->YReaction, flameNodeSource->YReaction, sizeof( Double ) * nSpeciesInSystem ); 
    memcpy( flameNodeDest->currReacRate, flameNodeSource->currReacRate, sizeof( Double ) * nReactions ); 
    memcpy( flameNodeDest->dMdY, flameNodeSource->dMdY, sizeof( Double ) * nSpeciesInSystem * (nSpeciesInSystem+1) ); 

    // species
    memcpy( flameNodeDest->viscosity, flameNodeSource->viscosity, sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->heatCapacity, flameNodeSource->heatCapacity, sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->conductivity, flameNodeSource->conductivity, sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->enthalpy, flameNodeSource->enthalpy, sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->diffusivity, flameNodeSource->diffusivity, sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->diffTherm, flameNodeSource->diffTherm, sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->productionRate, flameNodeSource->productionRate, sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->tempProp, flameNodeSource->tempProp, sizeof( Double ) ); 
    memcpy( flameNodeDest->pressureProp, flameNodeSource->pressureProp, sizeof( Double ) ); 

    // properties
    memcpy( flameNodeDest->mixViscosity, flameNodeSource->mixViscosity, sizeof( Double ) ); 
    memcpy( flameNodeDest->mixDensity, flameNodeSource->mixDensity, sizeof( Double ) ); 
    memcpy( flameNodeDest->mixConductivity, flameNodeSource->mixConductivity, sizeof( Double ) ); 
    memcpy( flameNodeDest->mixHeatCapacity, flameNodeSource->mixHeatCapacity, sizeof( Double ) ); 
    memcpy( flameNodeDest->mixMolarMass, flameNodeSource->mixMolarMass, sizeof( Double ) ); 

    // flame
    memcpy( flameNodeDest->Y[kCurr], flameNodeSource->Y[kCurr], sizeof( Double ) * nSpecies ); 
    memcpy( flameNodeDest->temp, flameNodeSource->temp, sizeof( Double ) ); 
    if ( fProperties->GetRadiation() ) {
        memcpy( flameNodeDest->radiation, flameNodeSource->radiation, sizeof( Double ) );
        memcpy( flameNodeDest->kappa, flameNodeSource->kappa, sizeof( Double ) );
    }
    memcpy( flameNodeDest->diffCorr, flameNodeSource->diffCorr, sizeof( Double ) ); 
}

    template<class Species>
void T1DFlame<Species>::CompareFlameNode( TFlameNodePtr flameNodeSource, TFlameNodePtr flameNodeDest )
{
    int				nThirdBodies = fReaction.GetNOfThirdBodies();
    int				nReactions = fReaction.GetNOfReactions();
    int				nSpecies = fSpecies.GetNOfSpecies();
    int				nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    int				notEqual;

    // reaction
    notEqual = 0;
    notEqual += memcmp( flameNodeDest->tBodyConc, flameNodeSource->tBodyConc, sizeof( Double ) * nThirdBodies );
    notEqual += memcmp( flameNodeDest->rateCoeff, flameNodeSource->rateCoeff, sizeof( Double ) * nReactions ); 
    notEqual += memcmp( flameNodeDest->tempReaction, flameNodeSource->tempReaction, sizeof( Double ) ); 
    notEqual += memcmp( flameNodeDest->currRateCoeff, flameNodeSource->currRateCoeff, sizeof( Double ) * nReactions ); 
    notEqual += memcmp( flameNodeDest->kNewerThanW, flameNodeSource->kNewerThanW, sizeof( Flag ) ); 
    notEqual += memcmp( flameNodeDest->reactionRate, flameNodeSource->reactionRate, sizeof( Double ) * nReactions ); 
    for ( int j = 0; j < nReactions; ++j ) {
        if ( flameNodeDest->reactionRate[j] != flameNodeSource->reactionRate[j] ) {
            fprintf( stderr, "# reaction %d\t%g\t%g\n", j, flameNodeDest->reactionRate[j], flameNodeSource->reactionRate[j] );
        }
    }

    notEqual += memcmp( flameNodeDest->YReaction, flameNodeSource->YReaction, sizeof( Double ) * nSpeciesInSystem ); 
    notEqual += memcmp( flameNodeDest->currReacRate, flameNodeSource->currReacRate, sizeof( Double ) * nReactions ); 
    notEqual += memcmp( flameNodeDest->dMdY, flameNodeSource->dMdY, sizeof( Double ) * nSpeciesInSystem * (nSpeciesInSystem+1) ); 
    if ( notEqual ) {
        fprintf( stderr, "#error: something's fishy in T1DFlame<Species>::CompareFlameNode reaction\n" );
    }
    // species
    notEqual = 0;
    notEqual += memcmp( flameNodeDest->viscosity, flameNodeSource->viscosity, sizeof( Double ) * nSpecies ); 
    notEqual += memcmp( flameNodeDest->heatCapacity, flameNodeSource->heatCapacity, sizeof( Double ) * nSpecies ); 
    notEqual += memcmp( flameNodeDest->conductivity, flameNodeSource->conductivity, sizeof( Double ) * nSpecies ); 
    notEqual += memcmp( flameNodeDest->enthalpy, flameNodeSource->enthalpy, sizeof( Double ) * nSpecies ); 
    notEqual += memcmp( flameNodeDest->diffusivity, flameNodeSource->diffusivity, sizeof( Double ) * nSpecies ); 
    notEqual += memcmp( flameNodeDest->diffTherm, flameNodeSource->diffTherm, sizeof( Double ) * nSpecies ); 
    notEqual += memcmp( flameNodeDest->productionRate, flameNodeSource->productionRate, sizeof( Double ) * nSpecies ); 
    notEqual += memcmp( flameNodeDest->tempProp, flameNodeSource->tempProp, sizeof( Double ) ); 
    notEqual += memcmp( flameNodeDest->pressureProp, flameNodeSource->pressureProp, sizeof( Double ) ); 
    if ( notEqual ) {
        fprintf( stderr, "#error: something's fishy in T1DFlame<Species>::CompareFlameNode species\n" );
    }

    // properties
    notEqual = 0;
    notEqual += memcmp( flameNodeDest->mixViscosity, flameNodeSource->mixViscosity, sizeof( Double ) ); 
    notEqual += memcmp( flameNodeDest->mixDensity, flameNodeSource->mixDensity, sizeof( Double ) ); 
    notEqual += memcmp( flameNodeDest->mixConductivity, flameNodeSource->mixConductivity, sizeof( Double ) ); 
    notEqual += memcmp( flameNodeDest->mixHeatCapacity, flameNodeSource->mixHeatCapacity, sizeof( Double ) ); 
    notEqual += memcmp( flameNodeDest->mixMolarMass, flameNodeSource->mixMolarMass, sizeof( Double ) ); 
    if ( notEqual ) {
        fprintf( stderr, "#error: something's fishy in T1DFlame<Species>::CompareFlameNode properties\n" );
    }

    // flame
    notEqual = 0;
    notEqual += memcmp( flameNodeDest->Y[kCurr], flameNodeSource->Y[kCurr], sizeof( Double ) * nSpecies ); 
    for ( int i = 0; i < nSpecies; ++i ) {
        if ( flameNodeDest->Y[kCurr][i] != flameNodeSource->Y[kCurr][i] ) 
            fprintf( stderr, "Ydest[%d] = %g \t Ysource[%d] = %g\n", i, flameNodeDest->Y[kCurr][i], i , flameNodeSource->Y[kCurr][i] );
    }
    notEqual += memcmp( flameNodeDest->temp, flameNodeSource->temp, sizeof( Double ) ); 
    if ( fProperties->GetRadiation() ) {
        notEqual |= memcmp( flameNodeDest->radiation, flameNodeSource->radiation, sizeof( Double ) );
        notEqual |= memcmp( flameNodeDest->kappa, flameNodeSource->kappa, sizeof( Double ) );
    }
    notEqual += memcmp( flameNodeDest->diffCorr, flameNodeSource->diffCorr, sizeof( Double ) ); 
    if ( notEqual ) {
        fprintf( stderr, "#error: something's fishy in T1DFlame<Species>::CompareFlameNode flame\n" );
    }
}
    template<class Species>
Flag T1DFlame<Species>::RHSAction(void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
    switch ( rhsMode ) {
        case kUpdate:
            UpdateSolutionOnePoint( nodeInfo->y, nodeInfo->gridPoint );
            UpdateThermoProps( fFlameNode, nodeInfo );


            break;
        case kDoNothing:
            break;
        case kSave:
            CopyFlameNode( fFlameNode, fFlameNodeSaved );
            return FALSE;
        case kRestore:
            CopyFlameNode( fFlameNodeSaved, fFlameNode );
            return FALSE;
        case kTest:
            CompareFlameNode( fFlameNodeSaved, fFlameNode );
            return FALSE;
        default:
            cerr << "#error: unknown RHSMode in T1DFlame<Species>::RHSAction" << NEWL;
            exit(2);
    }
    Clear1DArray( nodeInfo->rhs, nodeInfo->nOfEquations );

    return TRUE;
}
    template<class Species>
Flag T1DFlame<Species>::RHSAction( NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
    switch ( rhsMode ) {
        case kUpdate:
            UpdateSolutionOnePoint( nodeInfo->y, nodeInfo->gridPoint );
            UpdateThermoProps( fFlameNode, nodeInfo );
            break;
        case kDoNothing:
            break;
        case kSave:
            CopyFlameNode( fFlameNode, fFlameNodeSaved );
            return FALSE;
        case kRestore:
            CopyFlameNode( fFlameNodeSaved, fFlameNode );
            return FALSE;
        case kTest:
            CompareFlameNode( fFlameNodeSaved, fFlameNode );
            return FALSE;
        default:
            cerr << "#error: unknown RHSMode in T1DFlame<Species>::RHSAction" << NEWL;
            exit(2);
    }
    Clear1DArray( nodeInfo->rhs, nodeInfo->nOfEquations );

    return TRUE;
}

    template<class Species>
void T1DFlame<Species>::FillJacDiffusion( int nVariable, int nEquation, Double constCoeff, Double *diffCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
    // fills the jacobian with     constCoeff * d/dy ( rho * diffCoeff * df/dy)

    Double	*density = fFlameNode->mixDensity;
    Double	diffPlus = constCoeff * ( density[kCurr] * diffCoeff[kCurr] + density[kNext] * diffCoeff[kNext] );
    Double	diffMinus = constCoeff * ( density[kPrev] * diffCoeff[kPrev] + density[kCurr] * diffCoeff[kCurr] );
    Double	hm = nodeInfo->hm;
    Double	h = nodeInfo->h;

    if ( sign == kPositive ) {
        nodeInfo->a[nVariable][nEquation] -= ( hm * diffPlus + h * diffMinus );
        if ( !nodeInfo->lastPoint ) {
            nodeInfo->b[nVariable][nEquation] += hm * diffPlus;
        }
        if ( !nodeInfo->firstPoint ) {
            nodeInfo->c[nVariable][nEquation] += h * diffMinus;
        }
    }
    else {
        nodeInfo->a[nVariable][nEquation] += ( hm * diffPlus + h * diffMinus );
        if ( !nodeInfo->lastPoint ) {
            nodeInfo->b[nVariable][nEquation] -= hm * diffPlus;
        }
        if ( !nodeInfo->firstPoint ) {
            nodeInfo->c[nVariable][nEquation] -= h * diffMinus;
        }
    }
}

    template<class Species>
Double T1DFlame<Species>::StandardDiffusion( int nVariable, Double *diffCoeff, NodeInfoPtr nodeInfo )
{
    // returns finite difference approximation of   d/dy( rho * diffCoeff * df/fy )

    Double	yPrev = nodeInfo->yPrev[nVariable];
    Double	y = nodeInfo->y[nVariable];
    Double	yNext = nodeInfo->yNext[nVariable];
    Double	*density = fFlameNode->mixDensity;
    Double	diffPlus = density[kCurr] * diffCoeff[kCurr] + density[kNext] * diffCoeff[kNext];
    Double	diffMinus = density[kPrev] * diffCoeff[kPrev] + density[kCurr] * diffCoeff[kCurr];
    Double	hm = nodeInfo->hm;
    Double	h = nodeInfo->h;

    return ( ( diffPlus * hm * ( yNext - y ) + diffMinus * h * ( yPrev - y ) ) 
            / nodeInfo->hnenn );
}

    template<class Species>
void T1DFlame<Species>::FillJacMixFracDiffusion( int nVariable, int nEquation, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
    //	fill jacobian with  constCoeff * d/dy( rho * lambda / cp * dZ/dy )

    Double	*lambda = fFlameNode->mixConductivity;
    Double	*rho = fFlameNode->mixDensity;
    Double	*cp = fFlameNode->mixHeatCapacity;
    Double	diffPlus = constCoeff * ( rho[kCurr] * lambda[kCurr] / cp[kCurr]
            + rho[kNext] * lambda[kNext] / cp[kNext] );
    Double	diffMinus = constCoeff * ( rho[kPrev] * lambda[kPrev] / cp[kPrev]
            + rho[kCurr] * lambda[kCurr] / cp[kCurr] );
    Double	hm = nodeInfo->hm;
    Double	h = nodeInfo->h;

    if ( sign == kPositive ) {
        nodeInfo->a[nVariable][nEquation] -= ( hm * diffPlus + h * diffMinus );
        if ( !nodeInfo->lastPoint ) {
            nodeInfo->b[nVariable][nEquation] += hm * diffPlus;
        }
        if ( !nodeInfo->firstPoint ) {
            nodeInfo->c[nVariable][nEquation] += h * diffMinus;
        }
    }
    else {
        nodeInfo->a[nVariable][nEquation] += ( hm * diffPlus + h * diffMinus );
        if ( !nodeInfo->lastPoint ) {
            nodeInfo->b[nVariable][nEquation] -= hm * diffPlus;
        }
        if ( !nodeInfo->firstPoint ) {
            nodeInfo->c[nVariable][nEquation] -= h * diffMinus;
        }
    }
}

    template<class Species>
Double T1DFlame<Species>::SecondDerivMixFracDiffusion( int nVariable, NodeInfoPtr nodeInfo )
{
    Double	*lambda = fFlameNode->mixConductivity;
    Double	*rho = fFlameNode->mixDensity;
    Double	*cp = fFlameNode->mixHeatCapacity;
    Double	diffPlus =  rho[kCurr] * lambda[kCurr] / cp[kCurr]
        + rho[kNext] * lambda[kNext] / cp[kNext];
    Double	diffMinus = rho[kPrev] * lambda[kPrev] / cp[kPrev]
        + rho[kCurr] * lambda[kCurr] / cp[kCurr];
    Double	hm = nodeInfo->hm;
    Double	h = nodeInfo->h;
    Double	yPrev = nodeInfo->yPrev[nVariable];
    Double	y = nodeInfo->y[nVariable];
    Double	yNext = nodeInfo->yNext[nVariable];

    return ( diffPlus * hm * yNext - ( diffPlus * hm + diffMinus * h ) * y + diffMinus * h * yPrev ) 
        / ( h * hm * ( h + hm ) );
}

#endif // TFLAME_HPP__
