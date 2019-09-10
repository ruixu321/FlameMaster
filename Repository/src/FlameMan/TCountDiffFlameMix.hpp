#ifndef TCOUNT_DIFF_FLAME_MIX_HPP__
#define TCOUNT_DIFF_FLAME_MIX_HPP__

#include "ListTool.h"
#include "Spline.h"

// For full transport model instead of constant Lewis numbers use the following three options
#undef MOLARDIFF
#undef CONVECTION
#undef COOLTRANSPORT

// Undefining CENTRALDIFF gives better convergence. However,
// small mass conservation errors occur due to different discretization
// of convection and convection term by the use of upwind differencing, which
// are in the order of 1e-3 and essentially don't influence the solution. These
// errors are first order discretization errors and can hence be linearly decreased
// with increasing number of gridpoints. However, for strongly varying Lewis numbers
// such as for H2 flames, undefining CENTRALDIFF can lead to severe errors.

// It seems to be the case that diffusivitycorrection can only run if mole fraction
// diffusion is used. This has been found for methane at unity Lewis numbers. Therefore
// diffcorr is disabled if MOLARDIFF is undefined

#define CENTRALDIFF

// scalar dissipation rate stuff
// define READCHI to read chi from a file. The file name is Chi.tout. Look at the format below
#undef READCHI

// define CHIKW for Kim, Williams expression for chi(Z)
#undef CHIKW

#define ENTFLUX
#define HEATCAPGRAD

// enables size dependent diffusion in soot moment equations, see paper Pitsch, Riessmeier, Peters, CST.
#undef SIZEDEPDIFFUSION


    template <typename Species>
void TCountDiffFlameMix<Species>::InitTCountDiffFlameMix( void )
{

    cout << "start IniTCountDiffFlameMix " << endl;
    int i;
    TBVPSolverPtr	solver = GetSolver();
    TNewtonPtr		bt = solver->bt;
    TGridPtr		fine = bt->GetGrid()->GetFine();
    TGridPtr		coarse = bt->GetGrid()->GetCoarse();
    int				nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();

    fDeltaArcLength = 0.0;
    fdTds = 0.0;
    fdlnChids = 0.0;
    fDeltaTRef = 20.0;
    fNFlameletsCount = 0;
    fMaxFlamelets = 100;
    fDeltaChiRef = log( 1.3 ); // will give chi/chiold = 1.3
    fArcLengthContin = fInputData->fArcLengthContin;
    fRadiativeFrac = fInputData->fRadiativeFrac;
    fRadiationName = fInputData->fRadiationName;
    fZlZr = fInputData->fZlZr;

    fCl = fInputData->fCl;
    fHv = fInputData->fHv;
    fT_B = fInputData->fT_B;

    if ( fSoot ) {
        fSoot->SetMomentsOffset( fSootMoments );
    }


    fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];


    // fTemperature is an index where to find information about temperature (like a dictionnary)

    fVariableNames[fTemperature] = new char[2];
    strcpy( fVariableNames[fTemperature], "T" );


    fPrintMolarFractions = fInputData->fPrintMolarFractions;


    // information about species

    for ( i = 0; i < nSpeciesInSystem; ++i ) {
        fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies.GetNames()[i] ) + 1];
        strcpy( fVariableNames[fFirstSpecies + i], fSpecies.GetNames()[i] );
    }

    // information about chi
    fVariableNames[fLnChi] = new char[6];
    strcpy( fVariableNames[fLnChi], "lnChi" );

    // information about soot
    if ( fSoot ) {
        int	offset = fSoot->GetOffsetSootMoments();
        for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
            fVariableNames[offset + i] = new char[8];
            sprintf( fVariableNames[offset + i], "M%dORHO", i );
        }
    }


    cout << "fSolLnChi " << endl;
    //	always solve for some extra equation now
    fSolLnChi = NewVector( bt->GetMaxGridPoints() + 2 );
    fSolLnChi->vec = &fSolLnChi->vec[kNext];
    fSolLnChi->len -= 2;
    fSavedLnChi = NewVector( bt->GetMaxGridPoints() + 2 );
    fSavedLnChi->vec = &fSavedLnChi->vec[kNext];
    fSavedLnChi->len -= 2;

  
#ifdef READCHI
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
#endif

    //	set dissipation rate vector from input file or command line;
    if ( fInputData->fParameterComm >= 0.0 ) {
        fDissRate = NewVector( 1 );
        fDissRate->vec[0] = fInputData->fParameterComm;
        fDissRate->len = 0;
        cerr << "initial scalar dissipation rate from command line is " << fInputData->fParameterComm << NEWL;
    }
    else {
        if ( fInputData->fDissRate ) {
            fDissRate = NewVector( fInputData->fDissRate->len );
            copy_vec( fDissRate, fInputData->fDissRate );
            fDissRate->len = 0;
            cerr << "initial scalar dissipation rate is from FM input " << fInputData->fDissRate->vec[0] << NEWL;
        }
        else {
            cout << "DissRate from start profiles " << endl;
            // set dissrate vector from start profiles file
            fDissRate = NULL;
        }
    }

    if ( GetArcLengthCont() ) {

        cout << "GetArcLengthCont" << endl;
        if ( fDissRate && fDissRate->phys_len > 1 ) {
            if (fDissRate->vec[1] >= fDissRate->vec[0] ) {
                fArcUp = TRUE;
            }
            else {
                fArcUp = FALSE;
            }
        }
        else {
            fArcUp = TRUE;
        }
    }

    if ( fUseNumericalJac ) {
        cout << "fUseNumericalJac" << endl;
        bt->SetUtFuncs( NULL, NULL, NULL
                , CountDiffMixRHSRest<TCountDiffFlameMix<Species> >, CountDiffMixRHSRest<TCountDiffFlameMix<Species> >, CountDiffMixRHSRest<TCountDiffFlameMix<Species> >
                , CountDiffMixOutput<TCountDiffFlameMix<Species> >, CountDiffMixPostIter<TCountDiffFlameMix<Species> >
                , SetCountDiffMixNodeInfo<TCountDiffFlameMix<Species> >, CountDiffMixPostConv<TCountDiffFlameMix<Species> >
                , GetCountDiffMixVarNames<TCountDiffFlameMix<Species> >
                , CountDiffMixUpdateBoundary<TCountDiffFlameMix<Species> >, CountDiffMixUpdateBoundary<TCountDiffFlameMix<Species> > );
    }
    else {
        bt->SetUtFuncs( CountDiffMixJacRest<TCountDiffFlameMix<Species> >, CountDiffMixJacRest<TCountDiffFlameMix<Species> >, CountDiffMixJacRest<TCountDiffFlameMix<Species> >
                , CountDiffMixRHSRest<TCountDiffFlameMix<Species> >, CountDiffMixRHSRest<TCountDiffFlameMix<Species> >, CountDiffMixRHSRest<TCountDiffFlameMix<Species> >
                , CountDiffMixOutput<TCountDiffFlameMix<Species> >, CountDiffMixPostIter<TCountDiffFlameMix<Species> >
                , SetCountDiffMixNodeInfo<TCountDiffFlameMix<Species> >, CountDiffMixPostConv<TCountDiffFlameMix<Species> >
                , GetCountDiffMixVarNames<TCountDiffFlameMix<Species> >
                , CountDiffMixUpdateBoundary<TCountDiffFlameMix<Species> >, CountDiffMixUpdateBoundary<TCountDiffFlameMix<Species> > );
    }
    cout << "InitialBC" << endl;
    SetInitialBC( fine, fInputData );
    SetInitialBC( coarse, fInputData );
    cout << "Read Start Profiles" << endl;
    ReadStartProfiles( fInputData );
    CheckBC();
    CheckInitialGuess();
    UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );

    fLnChiContStart = log( fDissRate->vec[0] );
    fTempContStart = fSolTemp->vec[fTmaxLoc];

    SaveSolution();

    fprintf( stderr, "ZRef = %g\n", GetZRef() );
    if ( !UseDiffCorr() ) {
        fprintf( stderr, "### WARNING: Diffusivity Correction disabled\n" );;
    }
    else {
#ifndef MOLARDIFF
#ifdef CONVECTION
        fprintf( stderr, "### WARNING: If MOLARDIFF is undefined, CONVECTION should also be undefined\n    -> CONVECTION undefined\n" );
#undef CONVECTION
#endif
        UnSetDiffCorr();
#endif
    }

    if ( fCl > 0.0 ) {   // means fuel is liquid
        Double	T_L = fine->GetYRight()->vec[fTemperature];
        Double	h_Ref = 0.0;

        Double	h2, cpSum, entSum, deltaT, dummyPressure = 0.0;
        Double	temp = MAX( 200.0, T_L - fHv / fCl );
        int		i, count = 0, nOfSpeciesIn = fSpecies.GetNSpeciesInSystem();
        int		gridPoints = fSolver->bt->GetGrid()->GetCurrentGrid()->GetNGridPoints();

        SetFlameNode( gridPoints );

        // entSum is vapor enthalpy at T_b
        for( i = 0, entSum = 0.0; i < nOfSpeciesIn; ++i ) {
            fSpecies.ComputeSpeciesProperties( fFlameNode, fT_B, dummyPressure, i );
            entSum += fFlameNode->Y[kCurr][i] * fFlameNode->enthalpy[i];
        }

        h2 = fCl * (T_L - fT_B) - fHv + entSum;

        do {
            for( i = 0, entSum = cpSum = 0; i < nOfSpeciesIn; ++i ) {
                fSpecies.ComputeSpeciesProperties( fFlameNode, temp, dummyPressure, i );
                cpSum += fFlameNode->Y[kCurr][i] * fFlameNode->heatCapacity[i];
                entSum += fFlameNode->Y[kCurr][i] * fFlameNode->enthalpy[i];
            }

            deltaT = -( entSum - h2 ) / cpSum;
            temp += deltaT;

            if ( ++count > 1000 ) {
                fprintf( stderr, "#Error: temperature iteration for h = %g, T = %g not converged\n", h2, temp );
                exit( 2 );
            }
        } while ( fabs( deltaT / temp ) > 1.0e-3 );

        fprintf(stderr, "Initial fuel temperature is T_f = %g K\n", temp );

        fine->GetYRight()->vec[fTemperature] = temp;
        UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );
    }

    cout << "End of Init" << endl;


    // try to create the spatial grid ?
    //

}

    template <typename Species>
TCountDiffFlameMix<Species>::~TCountDiffFlameMix( void )
{
    fSolLnChi->vec = &fSolLnChi->vec[kPrev];
    DisposeVector( fSolLnChi );
    fSavedLnChi->vec = &fSavedLnChi->vec[kPrev];
    DisposeVector( fSavedLnChi );

    if ( fDissRate ) {
        DisposeVector( fDissRate );
    }
    int	nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
        delete fVariableNames[i];
    }
    delete fVariableNames;
#ifdef READCHI
    DisposeVector( fChiCount );
    DisposeVector( fZCount );
#endif
}

#ifdef READCHI

    template <typename Species>
Double TCountDiffFlameMix<Species>::GetDissRate( Double z, Double rho, Double chiSt )
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

#else

    template <typename Species>
Double TCountDiffFlameMix<Species>::GetDissRate( Double z, Double rho, Double chiStoich )
{

    // just to test if we use the solution at the different orgin


    /*

      NEED TO FIND A BETTER WAY
       */

    if (fZlZr == "true"){
        double Zr = 0.566;
        double Zl = 0.0112;

        return chiStoich * DissRateFact( rho ) / DissRateFact( fRhoRef ) *  ExactChi( z ) / ExactChi( GetZRef() )/((Zr - Zl)*(Zr - Zl));

    }
    else{

        return chiStoich * DissRateFact( rho ) / DissRateFact( fRhoRef ) *  ExactChi( z ) / ExactChi( GetZRef() );

    }
}

#endif

    template <typename Species>
Double TCountDiffFlameMix<Species>::ExactChi( Double Z )
{

    // we are using this function
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
    if ( twoZ < 1.0e-7 ) {
        return 0.0;
    }

    i = 0;
    while ( erFunc[i] > twoZ ) ++i;

    twoErfc = Interpol( twoZ, erFunc[i-1], (i-1)*deltax, erFunc[i], i*deltax );

    return exp( -2.0 * twoErfc * twoErfc );
}

    template <typename Species>
Double TCountDiffFlameMix<Species>::Interpol( Double x, Double xOld, Double yOld, Double xNew, Double yNew )
{
    if ( yNew == yOld ) {
        return yNew;
    }
    return yOld + ( yNew - yOld ) / ( xNew - xOld ) * ( x - xOld );
}

    template <typename Species>
Double TCountDiffFlameMix<Species>::DissRateFact( Double rho )
{
#ifdef CHIKW
    Double	densRatio = sqrt( fRhoInf / rho );

    return 1.5 * ( densRatio + 1 ) * ( densRatio + 1 ) / ( 2.0 * densRatio + 1 );
#else
    return 1.0;
#endif
}

    template<typename Flame>
void CountDiffMixJacRest( void *object, NodeInfoPtr nodeInfo )
{
    Flame*	flame = ( Flame* )object;
    TFlameNodePtr	flameNode = flame->fFlameNode;
    int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
    int 	fTemperature = flame->GetOffsetTemperature();
    int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
    int 	fLnChi = flame->GetOffsetLnChi();
    int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
    int		lastSpeciesEq = fFirstSpecies + nSpeciesInSystem;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = nodeInfo->hnenn;
    Double	**a = nodeInfo->a;
    Double	*yPrev = nodeInfo->yPrev;
    Double	*y = nodeInfo->y;
    Double	*yNext = nodeInfo->yNext;
    Double	*enthalpy = flameNode->enthalpy;
    Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
    Double	**dMdY = flameNode->dMdY;
    Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
    Double	*productionRate = flameNode->productionRate;
    Double	*heatCapacity = flameNode->heatCapacity;
    Double	*lewis = flame->GetSpecies()->GetLewisNumber()->vec;
    Double	mixDensity = flameNode->mixDensity[kCurr];
    Double	speciesCoeff = 0.5 * mixDensity * flame->GetDissRate( nodeInfo->x[kCurr]
            , mixDensity, exp( y[fLnChi] ));
    Double	energyCoeff = speciesCoeff * flameNode->mixHeatCapacity[kCurr];

#ifdef COOLTRANSPORT
    Double	coeff = flameNode->mixConductivity[kCurr]
        / ( flameNode->mixHeatCapacity[kCurr] * flameNode->mixDensity[kCurr] );
    for ( speciesEq = 0; speciesEq < nSpeciesInSystem; ++speciesEq ) {
        lewis[speciesEq] = coeff / flameNode->diffusivity[speciesEq];
    }
#endif

    flame->FilldMdYOnePoint( flameNode );
    flame->FilldMdTOnePoint( flameNode );

    if ( flame->GetSoot() ) {
        if ( !flame->fUseNumericalDM )  {
            flame->GetSoot()->UpdateJacobian( flame );
        }
        flame->GetSoot()->FillJacobi( flame, nodeInfo, kMixtureFraction );

        int		nSootMoments = flame->GetSoot()->GetNSootMoments();
        int		sootOff = flame->fSootMoments;
        Double	sootCoeff = flame->GetDissRate() / ( 2.0 * flame->GetSoot()->GetLewis1() );

        for ( int i = 0; i < nSootMoments; ++i ) {
#ifdef SIZEDEPDIFFUSION
            int		ioff, joff;
            Double	fracIndex;
            Double	wMinus = 2.0 * h;
            Double	wCurr = - 2.0 * ( h + hm );
            Double	wPlus = 2.0 * hm;
            ioff = sootOff + i;
            fracIndex = i - 2.0 / 3.0;
            for ( int j = 0; j < nSootMoments; ++j ) {
                joff = sootOff + j;
                a[joff][ioff] += sootCoeff * wCurr / y[joff]
                    * flame->GetSoot()->GetAlphaI2( j, fracIndex )
                    * flame->GetSoot()->FracMom2( fracIndex, &y[sootOff] );
                if ( !nodeInfo->lastPoint ) {
                    b[joff][ioff] += sootCoeff * wPlus / yNext[joff]
                        * flame->GetSoot()->GetAlphaI2( j, fracIndex )
                        * flame->GetSoot()->FracMom2( fracIndex, &yNext[sootOff] );
                }
                if ( !nodeInfo->firstPoint ) {
                    c[joff][ioff] += sootCoeff * wMinus / yPrev[joff]
                        * flame->GetSoot()->GetAlphaI2( j, fracIndex )
                        * flame->GetSoot()->FracMom2( fracIndex, &yPrev[sootOff] );
                }
            }
#else
            FillJacSecondDerivCentral( sootOff+i, sootOff+i, sootCoeff, nodeInfo );
#endif
        }
    }

    //	species equations
    for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
        speciesIndexEq = speciesEq - fFirstSpecies;
        FillJacSecondDerivCentral( speciesEq, speciesEq, speciesCoeff / lewis[speciesIndexEq], nodeInfo );
        a[fTemperature][speciesEq] -= speciesCoeff / ( y[fTemperature] * lewis[speciesIndexEq] )
            * SecondDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) * hnenn;
        a[fTemperature][speciesEq] += dMdT[speciesIndexEq] * hnenn;
        for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq; ++speciesVar ) {
            speciesIndexVar = speciesVar - fFirstSpecies;
            a[speciesVar][speciesEq] += dMdY[speciesIndexVar][speciesIndexEq] * hnenn;
        }
    }

    //	energy equation
    FillJacSecondDerivCentral( fTemperature, fTemperature, energyCoeff, nodeInfo );
    a[fTemperature][fTemperature] -= energyCoeff / y[fTemperature]
        * SecondDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) * hnenn;

    for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
        a[fTemperature][fTemperature] -= dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
        a[fTemperature][fTemperature] -= productionRate[speciesIndexEq] * heatCapacity[speciesIndexEq] * hnenn;
        for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
            a[speciesIndexVar + fFirstSpecies][fTemperature] -= dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
        }
    }



    if ( flame->fProperties->GetRadiation() ) {
        flame->fProperties->GetRadiation()->FillJacRadiation( 1.0, flame, nodeInfo );
    }
}

    template<typename Flame>
void CountDiffMixRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
    Flame*	flame = ( Flame* )object;

   
    if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
        return;
    }

    TFlameNodePtr	flameNode = flame->fFlameNode;
    int		speciesEq, speciesInd;
    int		M = nodeInfo->nOfEquations;
    int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
    int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
    int 	fTemperature = flame->GetOffsetTemperature();
    int 	fLnChi = flame->GetOffsetLnChi();
    int		lastSpeciesEq = fFirstSpecies + nSpeciesInSystem;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = nodeInfo->hnenn;
    Double	*yPrev = nodeInfo->yPrev;
    Double	*y = nodeInfo->y;
    Double	*yNext = nodeInfo->yNext;
    Double	*rhs = nodeInfo->rhs;
    Double	*enthalpy = flameNode->enthalpy;
    Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
    Double	*productionRate = flameNode->productionRate;
    Double	*lewisIn = flame->GetSpecies()->GetLewisNumber()->vec;
    Double	*lewis = New1DArray( nSpeciesInSystem );
    Double	*lewisPrev = New1DArray( nSpeciesInSystem );
    Double	*lewisNext = New1DArray( nSpeciesInSystem );
    Double	MMPrev = flameNode->mixMolarMass[kPrev];
    Double	MM = flameNode->mixMolarMass[kCurr];
    Double	MMNext = flameNode->mixMolarMass[kNext];
    Double	speciesCoeff = 0.5 * flameNode->mixDensity[kCurr]
        * flame->GetDissRate( nodeInfo->x[kCurr]
                , flameNode->mixDensity[kCurr], exp( y[fLnChi] ) );
    Double	energyCoeff = speciesCoeff * flameNode->mixHeatCapacity[kCurr];
    Double	paramMDiff = 1.0; //flame->GetSolver()->bt->GetParameter();
    Double	paramDiffTerms = 1.0; //flame->GetSolver()->bt->GetParameter();
    Double	paramDiffCorr = 1.0; //flame->GetSolver()->bt->GetParameter();
    Double	paramNow = flame->GetSolver()->bt->GetParameter();
    Double	paramCool = 1.0;//flame->GetSolver()->bt->GetParameter();
    Double	Le_ZPrev = 1.0 / ( paramNow + ( 1.0 - paramNow ) );
    Double	Le_ZCurr = 1.0 / ( paramNow + ( 1.0 - paramNow ) );
    Double	Le_ZNext = 1.0 / ( paramNow + ( 1.0 - paramNow ) );

#ifdef COOLTRANSPORT
    Double	coeff = flameNode->mixConductivity[kCurr]
        / ( flameNode->mixHeatCapacity[kCurr] * flameNode->mixDensity[kCurr] );
    Double	coeffPrev = flameNode->mixConductivity[kPrev]
        / ( flameNode->mixHeatCapacity[kPrev] * flameNode->mixDensity[kPrev] );
    Double	coeffNext = flameNode->mixConductivity[kNext]
        / ( flameNode->mixHeatCapacity[kNext] * flameNode->mixDensity[kNext] );
#endif
    for ( speciesEq = 0; speciesEq < nSpeciesInSystem; ++speciesEq ) {
#ifdef COOLTRANSPORT
        lewis[speciesEq] = coeff / ( ( flameNode->diffusivity[speciesEq]==0.0)?0.001:flameNode->diffusivity[speciesEq] ) / Le_ZCurr;
        lewisPrev[speciesEq] = coeffPrev / ( ( flameNode->diffusivityPrev[speciesEq]==0.0)?0.001:flameNode->diffusivityPrev[speciesEq] ) / Le_ZPrev;
        lewisNext[speciesEq] = coeffNext / ( ( flameNode->diffusivityNext[speciesEq]==0.0)?0.001:flameNode->diffusivityNext[speciesEq] ) / Le_ZNext;
#else
        lewis[speciesEq] = lewisPrev[speciesEq] = lewisNext[speciesEq] = lewisIn[speciesEq]/Le_ZCurr;
#endif
    }

    //	species equations
    Double	*rho = flameNode->mixDensity;
    Double	chiPrev = ( nodeInfo->firstPoint )
        ? flame->GetDissRate( flame->GetSolver()->bt->GetLeft()
                , flameNode->mixDensity[kPrev], exp( y[fLnChi] ) )
        : flame->GetDissRate( nodeInfo->x[kPrev], flameNode->mixDensity[kPrev], exp( y[fLnChi] ) );
    Double	chiCurr = flame->GetDissRate( nodeInfo->x[kCurr], flameNode->mixDensity[kCurr], exp( y[fLnChi] ) );
    Double	chiNext = ( nodeInfo->lastPoint )
        ? flame->GetDissRate( flame->GetSolver()->bt->GetRight()
                , flameNode->mixDensity[kNext], exp( y[fLnChi] ) )
        : flame->GetDissRate( nodeInfo->x[kNext], flameNode->mixDensity[kNext], exp( y[fLnChi] ) );
    Double	rhoChiPrev = rho[kPrev] * chiPrev;
    Double	rhoChi = rho[kCurr] * chiCurr;
    Double	rhoChiNext = rho[kNext] * chiNext;
    Double	lambdaOverCpLeZCurr = flameNode->mixConductivity[kCurr]
        / ( flameNode->mixHeatCapacity[kCurr] * Le_ZCurr );
    Double	lambdaOverCpLeZPrev = flameNode->mixConductivity[kPrev]
        / ( flameNode->mixHeatCapacity[kPrev] * Le_ZPrev );
    Double	lambdaOverCpLeZNext = flameNode->mixConductivity[kNext]
        / ( flameNode->mixHeatCapacity[kNext] * Le_ZNext );

    Double	dY_kdZ;
    Double	dLe_kdZ;

    Double	sum1 = 0.0;
    Double	sum2 = 0.0;
    Double	sum3 = 0.0;
    Double	sum4 = 0.0;
#ifndef CENTRALDIFF
    Double	sum4back = 0.0;
    Double	sum4forw = 0.0;
#endif
    Double	sum5Check = 0.0;
    if ( flame->UseDiffCorr() ) {
        for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
            speciesInd = speciesEq-fFirstSpecies;
            dY_kdZ = FirstDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
            dLe_kdZ = FirstDeriv( lewisPrev[speciesInd]
                    , lewis[speciesInd]
                    , lewisNext[speciesInd], hm, h );
            sum1 += SecondDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h )
                / lewis[speciesInd];
            sum2 += y[speciesEq] / lewis[speciesInd];
            sum3 += dLe_kdZ * dY_kdZ / ( lewis[speciesInd] * lewis[speciesInd] );
            sum4 += dY_kdZ / lewis[speciesInd];
#ifndef CENTRALDIFF
            sum4back += ( y[speciesEq] - yPrev[speciesEq] ) / hm / lewis[speciesInd];
            sum4forw += ( yNext[speciesEq] - y[speciesEq] ) / h / lewis[speciesInd];
#endif
            sum5Check += speciesCoeff * dLe_kdZ * y[speciesEq] / MM / ( lewis[speciesInd] * lewis[speciesInd] );
        }
    }

    // start species equations loop
    for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
        speciesInd = speciesEq-fFirstSpecies;
        // diffusion
        rhs[speciesEq] += speciesCoeff / lewis[speciesInd]
            * SecondDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
#ifdef MOLARDIFF
        // X diffusion term
        rhs[speciesEq] += paramMDiff * speciesCoeff / lewis[speciesInd]
            * y[speciesEq] / MM * SecondDeriv( MMPrev, MM, MMNext, hm, h );
#endif
        if ( flame->UseDiffCorr() ) {
            // diffcorr from diffusion term
            rhs[speciesEq] -= /*paramDiffCorr **/ speciesCoeff * y[speciesEq] * sum1;
#ifdef MOLARDIFF
            // diffcorr from X diffusion term
            rhs[speciesEq] -= paramDiffCorr * paramMDiff * speciesCoeff * y[speciesEq] / MM * sum2
                * SecondDeriv( MMPrev, MM, MMNext, hm, h );
#endif
        }


#ifdef CONVECTION
        Double	convTermMassFrac = 0.0;
        Double	convTermXDiff = 0.0;
        Double	convTermDiffCorrY = 0.0;
        Double	convTermDiffCorrXDiff = 0.0;

        // mass fraction convection velocity
        Double	convVeloY = 0.25 * ( 1.0 - 1.0 / lewis[speciesInd] )
            * (
                    FirstDeriv( rhoChiPrev, rhoChi, rhoChiNext, hm, h )

                    + rhoChi / lambdaOverCpLeZCurr
                    * FirstDeriv( lambdaOverCpLeZPrev, lambdaOverCpLeZCurr, lambdaOverCpLeZNext, hm, h )
              )

#	ifdef COOLTRANSPORT
            + paramCool * 0.5 * rhoChi / ( lewis[speciesInd] * lewis[speciesInd] )
            * FirstDeriv( lewisPrev[speciesInd], lewis[speciesInd], lewisNext[speciesInd], hm, h )
#	endif
            ;

        Double	convVeloM = 0.25 / lewis[speciesInd]
            * (
                    // 2.
                    2.0 * rhoChi / MM * FirstDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h )
                    // 3.
                    + y[speciesEq] * FirstDeriv( rhoChiPrev / MMPrev, rhoChi / MM, rhoChiNext / MMNext, hm, h )
                    // 4.
                    + rhoChi / lambdaOverCpLeZCurr * y[speciesEq]
                    * FirstDeriv( lambdaOverCpLeZPrev / MMPrev
                        , lambdaOverCpLeZCurr / MM
                        , lambdaOverCpLeZNext / MMNext, hm, h )
              )

#	ifdef COOLTRANSPORT
            - paramCool * 0.5 * rhoChi / ( lewis[speciesInd] * lewis[speciesInd] ) * y[speciesEq] / MM
            * FirstDeriv( lewisPrev[speciesInd], lewis[speciesInd], lewisNext[speciesInd], hm, h )
#	endif
            ;

        // convection terms
#	ifdef CENTRALDIFF
        // mass fraction convection
        convTermMassFrac = convVeloY * FirstDeriv( yPrev[speciesEq]
                , y[speciesEq], yNext[speciesEq], hm, h );
#		ifdef MOLARDIFF
        // molar mass convection
        convTermXDiff = convVeloM * FirstDeriv( MMPrev, MM, MMNext, hm, h );
#		endif
#	else
        if ( convVeloY < 0.0 ) {
            convTermMassFrac = convVeloY
                * ( y[speciesEq] - yPrev[speciesEq] ) / hm;
        }
        else {
            convTermMassFrac = convVeloY
                * ( yNext[speciesEq] - y[speciesEq] ) / h;
        }
#		ifdef MOLARDIFF
        if ( convVeloM > 0.0 ) {
            convTermXDiff = convVeloM * ( MM - MMPrev ) / hm;
        }
        else {
            convTermXDiff = convVeloM * ( MMNext - MM ) / h;
        }
#		endif
#	endif

        // diffusivity correction convection terms
        if ( flame->UseDiffCorr() ) {
#	ifdef CENTRALDIFF
            // mass fraction convection
#		ifdef COOLTRANSPORT
            // 1.
            convTermDiffCorrY -= paramCool * speciesCoeff * y[speciesEq] * sum3;
#		endif
            // 2. + 3.
            convTermDiffCorrY += 0.25 * sum4 * (
                    FirstDeriv( yPrev[speciesEq] * rhoChiPrev
                        , y[speciesEq] * rhoChi
                        , yNext[speciesEq] * rhoChiNext, hm, h )
                    + rhoChi / lambdaOverCpLeZCurr
                    * FirstDeriv( yPrev[speciesEq] * lambdaOverCpLeZPrev
                        , y[speciesEq] * lambdaOverCpLeZCurr
                        , yNext[speciesEq] * lambdaOverCpLeZNext, hm, h )
                    );
#	else
            Double convVeloYCorr23 = 0.25 * (
                    FirstDeriv( yPrev[speciesEq] * rhoChiPrev
                        , y[speciesEq] * rhoChi
                        , yNext[speciesEq] * rhoChiNext, hm, h )
                    + rhoChi / lambdaOverCpLeZCurr
                    * FirstDeriv( yPrev[speciesEq] * lambdaOverCpLeZPrev
                        , y[speciesEq] * lambdaOverCpLeZCurr
                        , yNext[speciesEq] * lambdaOverCpLeZNext, hm, h )
                    );
            if ( convVeloYCorr23 > 0.0 ) {
                convTermDiffCorrY = convVeloYCorr23 * sum4back;
            }
            else {
                convTermDiffCorrY = convVeloYCorr23 * sum4forw;
            }
#		ifdef COOLTRANSPORT
            convTermDiffCorrY -= paramCool * speciesCoeff * y[speciesEq] * sum3;
#		endif
#	endif
            // molar mass convection
#	ifdef MOLARDIFF
#		ifdef CENTRALDIFF
            Double	dMdZ = FirstDeriv( MMPrev, MM, MMNext, hm, h );

            for ( int kc = fFirstSpecies; kc < lastSpeciesEq; ++kc ) {
#			ifdef COOLTRANSPORT
                convTermDiffCorrXDiff -= paramCool * speciesCoeff * y[speciesEq] / MM
                    * y[kc] / ( lewis[kc-fFirstSpecies] * lewis[kc-fFirstSpecies] )
                    * FirstDeriv( lewisPrev[kc-fFirstSpecies]
                            , lewis[kc-fFirstSpecies]
                            , lewisNext[kc-fFirstSpecies], hm, h );
#			endif
                convTermDiffCorrXDiff += 0.25 / lewis[kc-fFirstSpecies]
                    * ( FirstDeriv( yPrev[speciesEq] / MMPrev * yPrev[kc] * rhoChiPrev
                                , y[speciesEq] / MM * y[kc] * rhoChi
                                , yNext[speciesEq] / MMNext * yNext[kc] * rhoChiNext, hm, h )
                            + rhoChi / lambdaOverCpLeZCurr
                            * FirstDeriv( yPrev[speciesEq] / MMPrev * yPrev[kc] * lambdaOverCpLeZPrev
                                , y[speciesEq] / MM * y[kc] * lambdaOverCpLeZCurr
                                , yNext[speciesEq] / MMNext * yNext[kc] * lambdaOverCpLeZNext, hm, h ) );

            }
            convTermDiffCorrXDiff *= dMdZ;
#		else
            Double convVeloMCorr12345 =
                // 1.
#			ifdef COOLTRANSPORT
                - paramCool * y[speciesEq] * sum5Check
#			endif
                + 0.25 * (
                        // 2.
                        2.0 * y[speciesEq] / MM * rhoChi * sum4
                        // 3.
                        + sum2 * FirstDeriv( yPrev[speciesEq] / MMPrev * rhoChiPrev
                            , y[speciesEq] / MM * rhoChi
                            , yNext[speciesEq] / MMNext * rhoChiNext, hm, h )
                        // 4.
                        + sum2 * rhoChi / lambdaOverCpLeZCurr
                        * FirstDeriv( lambdaOverCpLeZPrev * yPrev[speciesEq] / MMPrev
                            , lambdaOverCpLeZCurr * y[speciesEq]  / MM
                            , lambdaOverCpLeZNext * yNext[speciesEq] / MMNext, hm, h )
                        );

            if ( convVeloM > 0.0 ) {
                convTermDiffCorrXDiff = convVeloMCorr12345 * ( MM - MMPrev ) / hm;
            }
            else {
                convTermDiffCorrXDiff = convVeloMCorr12345 * ( MMNext - MM ) / h;
            }
#		endif
#	endif
        }

        rhs[speciesEq] += paramDiffTerms * ( -convTermMassFrac
                + /*paramMDiff **/ convTermXDiff
                - /*paramDiffCorr **/ convTermDiffCorrY
                - /*paramDiffCorr * paramMDiff **/ convTermDiffCorrXDiff );


#endif
        rhs[speciesEq] += productionRate[speciesInd];
    }

    //	energy equation
    rhs[fTemperature] += energyCoeff * Le_ZCurr * SecondDeriv( yPrev[fTemperature]
            , y[fTemperature], yNext[fTemperature], hm, h );

#ifdef CONVECTION
    // mass fraction convection velocity
    Double	convVeloT = flameNode->mixHeatCapacity[kCurr] * ( 0.25 * ( Le_ZCurr - 1.0 )
            * (
                FirstDeriv( rhoChiPrev, rhoChi, rhoChiNext, hm, h )

                + rhoChi / lambdaOverCpLeZCurr
                * FirstDeriv( lambdaOverCpLeZPrev, lambdaOverCpLeZCurr, lambdaOverCpLeZNext, hm, h )
              )

#	ifdef COOLTRANSPORT
            + paramCool * 0.5 * rhoChi
            * FirstDeriv( Le_ZPrev, Le_ZCurr, Le_ZNext, hm, h )
#	endif
            );

#	ifdef CENTRALDIFF
    rhs[fTemperature] += convVeloT * FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#	else
    if ( convVeloT < 0.0 ) {
        rhs[fTemperature] += convVeloT
            * ( y[fTemperature] - yPrev[fTemperature] ) / hm;
    }
    else {
        rhs[fTemperature] += convVeloT
            * ( yNext[fTemperature] - y[fTemperature] ) / h;
    }
#	endif
#endif

#ifdef ENTFLUX
    Double	*cp = flameNode->heatCapacity;
    Double	entFluxSum = 0.0;
    Double	cpMix = flameNode->mixHeatCapacity[kCurr];
    if ( !flame->UseDiffCorr() ) {
        cpMix = 0.0;
    }
#endif
    for ( speciesEq = 0; speciesEq < nSpeciesInSystem; ++speciesEq ) {
        speciesInd = speciesEq+fFirstSpecies;
#ifdef ENTFLUX
        entFluxSum += ( cpMix - cp[speciesEq] ) / lewis[speciesEq]
            * ( FirstDeriv( yPrev[speciesInd], y[speciesInd], yNext[speciesInd], hm, h )
#	ifdef MOLARDIFF
                    + y[speciesInd] / MM * FirstDeriv( MMPrev, MM, MMNext, hm, h )
#	endif
              );
#endif
        rhs[fTemperature] -= productionRate[speciesEq] * enthalpy[speciesEq];
    }
#ifdef ENTFLUX
    rhs[fTemperature] -= speciesCoeff * entFluxSum
        * FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#endif
#ifdef HEATCAPGRAD
    rhs[fTemperature] += speciesCoeff * Le_ZCurr
        * FirstDeriv( flameNode->mixHeatCapacity[kPrev]
                , flameNode->mixHeatCapacity[kCurr]
                , flameNode->mixHeatCapacity[kNext], hm, h )
        * FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#endif
    if ( flame->fProperties->GetRadiation() ) {
        if(flame->GetRadiationName() == "WSGG"){
            rhs[fTemperature] += flameNode->radiation[kCurr];
        }
        else if(flame->GetRadiationName() == "Grey"){
            rhs[fTemperature] += flameNode->radiation[kCurr];
        }
        else if(flame->GetRadiationName() == "SNB"){
            rhs[fTemperature] += flameNode->radiation[kCurr];
        }
        else if(flame->GetRadiationName() == "WSGGJohansson"){
            rhs[fTemperature] += flameNode->radiation[kCurr];
        }
        else if(flame->GetRadiationName() == "WSGGBordbar"){
            rhs[fTemperature] += flameNode->radiation[kCurr];
        }
        else{
            rhs[fTemperature] += flameNode->radiation[kCurr];
        }
        if ( flame->GetSoot() ) {
            rhs[fTemperature] -= flame->GetSoot()->GetSootRadiation( y[fTemperature]
                    , flameNode->moments );
        }
    }

    if ( flame->GetArcLengthCont() ) {
        if ( nodeInfo->gridPoint < flame->fTmaxLoc ) {
            rhs[fLnChi] += FirstDerivUpwind( yNext[fLnChi], y[fLnChi], h );
        }
        else if ( nodeInfo->gridPoint == flame->fTmaxLoc ) {
            Double	dT = (y[fTemperature] - flame->GetTempContStart()) / flame->GetDeltaTref();
            Double	dlnChi = (y[fLnChi] - flame->GetChiContStart()) / flame->GetDeltaChiref();
            if ( flame->GetDeltaArcLength() == 0.0 ) {
                rhs[fLnChi] += dlnChi;
            }
            else {
                rhs[fLnChi] += dT * flame->GetdTds() + dlnChi * flame->GetdlnChids()
                    - flame->GetDeltaArcLength();
            }
        }
        else {
            rhs[fLnChi] += FirstDerivUpwind( y[fLnChi], yPrev[fLnChi], hm );
        }
    }
    else {
        rhs[fLnChi] += y[fLnChi] - flame->GetChiContStart();
    }

    //	soot equations
    if ( flame->GetSoot() ) {
        int		nSootMoments = flame->GetSoot()->GetNSootMoments();
        int		sootOff = flame->fSootMoments;
        Double	sootCoeff = flame->GetDissRate( nodeInfo->x[kCurr]
                , flameNode->mixDensity[kCurr], exp( y[fLnChi] ) ) / ( 2.0 * flame->GetSoot()->GetLewis1() );

        flame->GetSoot()->FillRHS( flame, nodeInfo, kMixtureFraction );
        for ( int i = 0; i < nSootMoments; ++i ) {
#ifdef SIZEDEPDIFFUSION
            Double	fracIndex;
            Double	wMinus = 2.0 * h / hnenn;
            Double	wCurr = - 2.0 * ( h + hm ) / hnenn;
            Double	wPlus = 2.0 * hm / hnenn;
            fracIndex = i - 2.0 / 3.0;
            rhs[sootOff+i] += sootCoeff
                * ( wMinus * flame->GetSoot()->FracMom2( fracIndex, &yPrev[sootOff] )
                        + wCurr * flame->GetSoot()->FracMom2( fracIndex, &y[sootOff] )
                        + wPlus * flame->GetSoot()->FracMom2( fracIndex, &yNext[sootOff] ) );
#else
            rhs[sootOff+i] += sootCoeff
                * SecondDeriv( yPrev[sootOff+i], y[sootOff+i], yNext[sootOff+i], hm, h );
#endif
        }
    }

    TTimePtr tim = flame->GetSolver()->time;
    for ( speciesEq = 0; speciesEq < M; ++speciesEq ) {
        if ( flame->GetSolver()->bt->GetTimedepFlag()
                && !flame->GetSolver()->time->GetTimeConverged() ) {
            rhs[speciesEq] -= ( y[speciesEq] - tim->GetYOld()->mat[nodeInfo->gridPoint][speciesEq] )
                / tim->GetDeltaT();
        }
        rhs[speciesEq] *= - hnenn;
    }

    Free1DArray( lewis );
    Free1DArray( lewisNext );
    Free1DArray( lewisPrev );
}

    template <typename Species>
void TCountDiffFlameMix<Species>::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
    // fills the jacobian with     constCoeff * sum_j ( Y_k/Le_j d2Y_j/dy2 )

    if ( sign == kNegative ) {
        constCoeff *= -1.0;
    }

    int		i, lVar;
    int		nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    Double	*Y = fFlameNode->Y[kCurr];
    Double	coeff = constCoeff * Y[nVariable-fFirstSpecies];
    Double	sumY = 0.0;
    Double	*lewis = fSpecies.GetLewisNumber()->vec;

    for ( i = 0; i < nSpeciesInSystem; ++i ) {
        sumY += Y[i];
    }
    coeff /= sumY;

    for ( i = 0; i < nSpeciesInSystem; ++i ) {
        lVar = fFirstSpecies + i;
        // d/dY_l
        FillJacSecondDerivCentral( lVar, nVariable, coeff/lewis[i], nodeInfo, sign );
    }
}

    template <typename Species>
void TCountDiffFlameMix<Species>::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
{
    int 				i, j, k;
    TBVPSolverPtr		solver = GetSolver();
    TNewtonPtr			bt = solver->bt;
    TAdaptiveGridPtr	adapGrid = bt->GetGrid();
    TGridPtr			grid = adapGrid->GetFine();
    int					variables = bt->GetNVariables();
    int					nGridPoints = grid->GetNGridPoints();
    int					maxGridPoints = bt->GetMaxGridPoints();
    int					initialGridPoints = bt->GetInitialGridpoints();
    int					nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    int					speciesIndex;
    MatrixPtr			yMat = grid->GetY();
    VectorPtr			yLeftVec = grid->GetYLeft();
    VectorPtr			yRightVec = grid->GetYRight();
    Double			 	*yLeft = yLeftVec->vec;
    Double 				*yRight = yRightVec->vec;
    Double				left = 0.0;
    Double				right = 1.0;
    Double				*locX = grid->GetX()->vec;
    int					gridPointsIn = sp->gridPoints;	// including boundaries
    Double				*yWork = adapGrid->GetWorkVector()->vec;
    Flag				ZSet = FALSE;
    Flag				chooseInputGrid = FALSE;
    Double				*xIn = new Double[gridPointsIn];
    if ( !xIn ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );
    Double				*yIn =  new Double[gridPointsIn];
    if ( !yIn ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );
    Double				*yInFloat = sp->data;
    Double				**y = grid->GetY()->mat;
    int					variable;
    char				*string = sp->labels;
    SplinePtr			theSpline = NULL;
    Double				leftSlope;
    Double				rightSlope;
    int					oxidizerSide; // this program assumes that oxidizerSide = kRight
    FILE				*fp;
    Double				dissRateIn;
    struct _parameter	*param = GetParameter( "chi" );

    if (!param) param = GetParameter( "chi_ref" );
    if (!param) param = GetParameter( "chi_st" );

    if ( !fDissRate || fDissRate->vec[0] < 0.0 ) {
        // get scalar dissipation rate from input file
        if ( param ) {
            dissRateIn = (Double)param->what.quantity.value;
        }
        else { // choose default
            cerr << "#warning: no value for 'scalar dissipation rate' in inputfile" << NEWL;
            dissRateIn = 10.0;
        }

        if ( !fDissRate ) {
            fDissRate = NewVector( 1 );
        }
        fDissRate->vec[0] = dissRateIn;
        fDissRate->len = 0;
        cerr << "initial scalar dissipation rate is " << GetDissRate() << NEWL;
    }


    //	choose grid from input or equidistant
    if ( gridPointsIn <= maxGridPoints+2 && gridPointsIn > initialGridPoints+2 && ( gridPointsIn % 2 ) != 0 ) {
        grid->AdjustNGridPoints( gridPointsIn-2 );
        solver->UpdateAllDimensions( gridPointsIn-2 );
        nGridPoints = grid->GetNGridPoints();
        chooseInputGrid = TRUE;
    }
    else {
        nGridPoints = initialGridPoints;
        grid->AdjustNGridPoints( nGridPoints );
        solver->UpdateAllDimensions( nGridPoints );
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
            bt->SetLeft( left );
            bt->SetRight( right );
            if ( chooseInputGrid ) {
                cerr << "choose inputGrid" << NEWL;
                if ( oxidizerSide == kRight ) { // turn vector
                    for ( j = 0; j < gridPointsIn-2; ++j ) {
                        locX[gridPointsIn-j-3] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
                    }
                }
                else { // copy vector
                    for ( j = 0; j < gridPointsIn-2; ++j ) {
                        locX[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
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
                grid->Make_equi_Grid();
                ZSet = TRUE;
            }
        }
        string += strlen(string) + 1;
    }

    // set default values
    for ( i = 0; i < nGridPoints; ++i ) {
        for ( j = 0; j < variables; ++j ) {
            y[i][j] = yLeft[j] + ( yRight[j] - yLeft[j] ) / ( right - left ) * locX[i];
        }
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
        else if ( GetSoot() && strncmp( string, "conc-soot", 9 ) == 0 ) {
            string += 9;
            int	num = atoi( string );
            if ( num < GetSoot()->GetNSootMoments() ) {
                variable = num + GetSoot()->GetOffsetSootMoments();
            }
            else {
                string += strlen(string) + 1;
                continue;
            }
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
        else {
            string += strlen(string) + 1;
            continue;
        }

        string += strlen(string) + 1;
        if ( chooseInputGrid ) {
            if ( oxidizerSide == kRight ) { // turn vector
                for ( k = 0; k < gridPointsIn-2; ++k ) {
                    y[k][variable] = yInFloat[(i+1)*gridPointsIn - k-2];	// copy workspace to vector of solution
                }
            }
            else {
                for ( k = 0; k < gridPointsIn-2; ++k ) {
                    y[k][variable] = yInFloat[i*gridPointsIn + k+1];	// copy workspace to vector of solution
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
            SplineInterpolate( theSpline, locX, yWork, nGridPoints );
            if ( oxidizerSide == kRight ) { // turn vector
                for ( k = 0; k < nGridPoints; ++k ) {
                    y[k][variable] = yWork[nGridPoints-k-1];	// copy workspace to vector of solution
                }
            }
            else {
                for ( k = 0; k < nGridPoints; ++k ) {
                    y[k][variable] = yWork[k];	// copy workspace to vector of solution
                }
            }
        }
        //	set bc
        if ( oxidizerSide == kLeft ) {
            if ( variable == fTemperature ) {
                if ( yLeft[fTemperature] <= 0.0 ) {
                    yLeft[fTemperature] = yInFloat[i*gridPointsIn];
                }
                if ( yRight[fTemperature] <= 0.0 ) {
                    yRight[fTemperature] = yInFloat[(i+1)*gridPointsIn-1];
                }
            }
        }
        else {
            if ( variable == fTemperature ) {
                if ( yLeft[fTemperature] <= 0.0 ) {
                    yLeft[fTemperature] = yInFloat[(i+1)*gridPointsIn-1];
                }
                if ( yRight[fTemperature] <= 0.0 ) {
                    yRight[fTemperature] = yInFloat[i*gridPointsIn];
                }
            }
        }
    }

    // set left to zero

    if ( GetPressure() <= 0.0 ) {
        param = GetParameter( "pressure" );
        Double	thePressure;
        if ( param ) {
            thePressure = (Double)param->what.quantity.value;
            if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
                thePressure *= 1.0e5;
            }
            SetPressure( thePressure );
            fprintf( stderr, "%s%g%s\n", "initial pressure is ", GetPressure()/1.0e5, " bar"  );
        }
        else { // exit
            cerr << "#error: no value for 'pressure' in inputfile" << NEWL;
            exit(2);
        }
    }


    // set initial Boundary values
    //	UpdateThermoProps();
    UpdateSolution( yMat, yLeftVec, yRightVec );
    CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );

    // convert M_i to M_i/rho
    if ( GetSoot() ) {
        Double	mixMolarMass;

        for ( k = 0; k < nGridPoints; ++k ) {
            fProperties->ComputeMixtureMolarMass( mixMolarMass, &y[k][fFirstSpecies]
                    , fSpecies.GetMolarMass()->vec, fSpecies.GetNSpeciesInSystem() );
            for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
                y[k][fSootMoments+i] *=  RGAS * y[k][fTemperature]
                    / ( GetPressure() * mixMolarMass );
            }
        }
    }

    for ( i = 0; i < nGridPoints; ++i ) {
        y[i][fLnChi] = log( fDissRate->vec[0] );
    }
    yLeft[fLnChi] = yRight[fLnChi] = log( fDissRate->vec[0] );

    CountDiffMixPostIter<TCountDiffFlameMix<Species> >( this );

    FreeSpline( theSpline );
    delete[] yIn;
    delete[] xIn;

    adapGrid->SetSolutionScaler();

    fp = GetOutfile( "initialguess", FileType::kData );
    bt->PrintSolution( locX, y, GetVariableNames(), fp );
    fclose(fp);

}

    template <typename Flame>
int CountDiffMixPostIter( void *object )
{
    Flame*	                flame = ( Flame* )object;
    int						nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
    int						fTemperature = flame->GetOffsetTemperature();
    int						fFirstSpecies = flame->GetOffsetFirstSpecies();
    TNewtonPtr				bt = flame->GetSolver()->bt;
    TGridPtr 				currGrid = bt->GetGrid()->GetCurrentGrid();
    int						nGridPoints = currGrid->GetNGridPoints();
    MatrixPtr				yMat = currGrid->GetY();
    VectorPtr				yLeftVec = currGrid->GetYLeft();
    VectorPtr				yRightVec = currGrid->GetYRight();
    Double					**y = yMat->mat;
    Double	                *yLeft = yLeftVec->vec;
    Double	                *yRight = yRightVec->vec;
    Double			        *x = currGrid->GetX()->vec;
    MatrixPtr               test = flame->GetSpecies()->GetDeltaI();

    string                  radiationName = flame->GetRadiationName();

    for ( int i = 0; i < nGridPoints; ++i ) {
        if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
            return 1;
        }
    }
    if ( flame->GetSoot() ) {
        flame->GetSoot()->PostIter( flame );
    }

    flame->UpdateSolution( yMat, yLeftVec, yRightVec );

    flame->SetTMaxLoc( flame->GetZRefLoc( currGrid->GetX() ) );

    //cout << " chi " << flame->GetDissRate() << endl;

    flame->UpdateThermoProps();    // The radiation calculation of the thin model will be done in this step

    if(radiationName == "WSGG" || radiationName == "Grey" || radiationName == "SNB" || 
        radiationName == "WSGGJohansson" || radiationName == "WSGGBordbar"){

        /* NEED TO ADD COMMENT HERE

        */
        // compute the physical grid

        vector<double> physicalGrid(nGridPoints+2);

        double DeltaX = 0.01;

       
        Double	**ent = flame->GetSpecies()->GetEnthalpy()->mat;
        Double	**prod = flame->GetSpecies()->GetProductionRate()->mat;
       
        vector<double> heatrel(nGridPoints);

        for (int k = 0; k < nGridPoints; ++k ) {
        heatrel[k] = 0.0;
          for (int i = 0; i < flame->GetSpecies()->GetNSpeciesInSystem(); ++i ) {
            heatrel[k] += ent[k][i] * prod[k][i];
          }
          heatrel[k] = - heatrel[k];
        }

        // bc
        physicalGrid[0] = DeltaX;

        double chiRef = flame->GetDissRate();

        for (int k = 0; k < nGridPoints; ++k){

            double chi = flame->GetDissRate( x[k], flame->GetProperties()->GetDensity()->vec[k], exp( y[k][flame->GetLnChi()] ) );

            double D = (flame->GetProperties()->GetConductivity()->vec[k])/
                (flame->GetProperties()->GetDensity()->vec[k]*flame->GetProperties()->GetHeatCapacity()->vec[k]);

            physicalGrid[k+1] = (x[k] - x[k-1])*sqrt(2*D/chi) + physicalGrid[k]; 

        }

        physicalGrid[nGridPoints+1] = physicalGrid[nGridPoints] + DeltaX;

        // checking the physical grid
        //for (int k=0; k< physicalGrid.size(); k++){

        //    cout << physicalGrid[k] << endl;
        //}
        flame->UpdateThermoProps(physicalGrid, DeltaX, chiRef, heatrel, radiationName);
        flame->UpdateThermoPropsRadiation();
    }

    flame->fRhoInf = flame->GetProperties()->GetDensity()->vec[kPrev];
    flame->fRhoRef = InterpolOne( flame->GetZRef(), currGrid->GetX()->vec
            , flame->GetProperties()->GetDensity()->vec, nGridPoints );

    return 0;
}

    template <typename Species>
Double TCountDiffFlameMix<Species>::GetZRef( void )
{
    return GetZStoich();
}

    template <typename Species>
void TCountDiffFlameMix<Species>::UpdateDimensions( int len )
{
    T1DFlame<Species>::UpdateDimensions( len );
    fSolLnChi->len = len;
}

    template <typename Species>
int TCountDiffFlameMix<Species>::GetZRefLoc( VectorPtr xVec )
{
    int		k = 0, nGridPoints = xVec->len;
    Double	*x = xVec->vec;
    Double	zRef = GetZRef();

    while ( x[k] <= zRef ) {
        ++k;
    };

    return k - 1; // return the location equal or left of zstoich
}

    template <typename Species>
void TCountDiffFlameMix<Species>::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
    int		nGridPoints = yMat->cols;
    Double	*lnChi = fSolLnChi->vec;
    Double	**y = yMat->mat;
    Double	*yLeft = yLeftVec->vec;
    Double	*yRight = yRightVec->vec;

    UpdateDimensions( nGridPoints );

    T1DFlame<Species>::UpdateSolution( yMat, yLeftVec, yRightVec );

    lnChi[kPrev] = yLeft[fLnChi];
    for ( int k = 0; k < nGridPoints; ++k ) {
        lnChi[k] = y[k][fLnChi];
    }
    lnChi[nGridPoints] = yRight[fLnChi];
}

    template <typename Species>
int	TCountDiffFlameMix<Species>::GetOffsetVVelocity( void )
{
    cerr << "#error: class has no member fVVelocity" << NEWL;
    exit( 2 );
    return 0;
}

    template <typename Species>
int	TCountDiffFlameMix<Species>::GetOffsetUVelocity( void )
{
    cerr << "#error: class has no member fUVelocity" << NEWL;
    exit( 2 );
    return 0;
}

    template <typename Species>
int TCountDiffFlameMix<Species>::GetOffsetMixFrac( void )
{
    cerr << "#error: class has no member fMixFrac" << NEWL;
    exit( 2 );
    return 0;
}

    template <typename Species>
int	TCountDiffFlameMix<Species>::GetOffsetFirstSpecies( void )
{
    return fFirstSpecies;
}

    template <typename Species>
int	TCountDiffFlameMix<Species>::GetOffsetTemperature( void )
{
    return fTemperature;
}

    template <typename Species>
int	TCountDiffFlameMix<Species>::GetOffsetLnChi( void )
{
    return fLnChi;
}

    template <typename Species>
ConstStringArray TCountDiffFlameMix<Species>::GetVariableNames( void )
{
    return (ConstStringArray)fVariableNames;
}

    template <typename Species>
int TCountDiffFlameMix<Species>::GetVariablesWithoutSpecies( void )
{
    return fVariablesWithoutSpecies;
}

    template <typename Species>
FILE *TCountDiffFlameMix<Species>::GetOutputFile( const char *head, const char *tail, const FileType type )
{
    int				fuelIndex = GetFuelIndex();
    char			*name = new char[64];
    FILE			*fp;
    char			**speciesNames = fSpecies.GetNames();
    int				tFuel = ( int ) fSolTemp->vec[fSolTemp->len];
    int				tOxidizer = ( int ) fSolTemp->vec[kPrev];
    Double			press = GetPressure() * 1.0e-5;

    sprintf( name, "%s%s%.8s_p%.2d_%.1dchi%05gtf%.4dto%.4d%s"
            , ( head ) ? head : "", ( head ) ? "_" : ""
            , speciesNames[fuelIndex]
            , ( int ) floor( press )	// in [bar]
            , ( int ) ( ( press - ( floor( press ) ) ) * 10 + 0.5 )
            , ( GetDissRate() )					// in [1/s]
            , ( int )( tFuel )							// in [K]
            , ( int )( tOxidizer ) 						// in [K]
            , ( tail ) ? tail : "" );

    fp = GetOutfile( name, type );
    delete[] name;

    return fp;
}

    template <typename Species>
void TCountDiffFlameMix<Species>::SetInitialBC( TGridPtr grid, TInputDataPtr inp )
{
    int					i;
    Double				mixMolarMass;
    int					nSpeciesInSystem = fSpecies.GetNSpeciesInSystem();
    //SpeciesPtr			species = inp->GetSpecies();					// ###
    // change left and right boundary conditions
    BoundaryInputPtr	right = inp->leftBoundary;
    BoundaryInputPtr	left = inp->rightBoundary;
    int					inpTOffset = inp->fTemperatureOffset;
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
    bcFlagLeft[fTemperature] = left->fBcFlag[inpTOffset];
    bcFlagRight[fTemperature] = right->fBcFlag[inpTOffset];
    for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
        bcFlagLeft[i] = left->fBcFlagSpecies;
        bcFlagRight[i] = right->fBcFlagSpecies;
    }

    // set value
    yleft[fTemperature] = left->fValue[inpTOffset];
    yright[fTemperature] = right->fValue[inpTOffset];

    bcLeft[fTemperature] = left->fValue[inpTOffset];
    bcRight[fTemperature] = right->fValue[inpTOffset];


    for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
        yleft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
        bcLeft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
    }
    if ( left->fMixtureSpecification == kMolarFraction ) {
        // first compute molar mass of mixture
        for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
            mixMolarMass += fSpecies.GetMolarMass()->vec[i] * yleft[i+fFirstSpecies];	// ###
        }
        // compute massfractions
        for ( i = 0; i < nSpeciesInSystem; ++i ) {
            yleft[i+fFirstSpecies] *= fSpecies.GetMolarMass()->vec[i] / mixMolarMass; 	// ###
            bcLeft[i+fFirstSpecies] = yleft[i+fFirstSpecies];
        }
        for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
            bcFlagLeft[i] = kMassFraction;
        }
    }

    for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
        yright[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
        bcRight[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
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

    template <typename Flame>
ConstStringArray GetCountDiffMixVarNames( void *object )
{
    Flame*	flame = ( Flame* )object;

    return flame->GetVariableNames();
}

    template <typename Flame>
void SetCountDiffMixNodeInfo( int k, void *object )
{
    Flame*	flame = ( Flame* )object;

    flame->SetFlameNode( k );
}

    template <typename Flame>
void CountDiffMixPostConv( void *object )
{
    Flame*	flame = ( Flame* )object;
    VectorPtr				dissRate = flame->GetDissRateVector();
    static Double			chiConv = -1.0;	//means not set
    static Double			chiNotConv = -1.0;	//means not set
    TNewtonPtr				bt = flame->GetSolver()->bt;
    TAdaptiveGridPtr		grid = bt->GetGrid();
    int						isConverged = bt->GetConvergeNewton();

    Double			*temp = flame->GetTemperature()->vec;
    TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
    Double			*x = currentGrid->GetX()->vec;
    int				gridPoints = currentGrid->GetNGridPoints();

    fprintf( stderr, "Tmax = %g @ Z = %g and Chi_ref = %g\n"
            , temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1]
            , x[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1], flame->GetDissRate() );
    if ( isConverged ) {
        flame->SaveSolution();
        if ( flame->fReactionFluxes ) {
            flame->ReactionFluxes( kPhysical );
            flame->GetReaction()->PrintReactionRates( flame );
            flame->fReaction.PrintRateCoeffs( flame );
            flame->fReaction.PrintDetailedHeatRelease( flame );
        }
        for ( int i = 0; i < flame->fNSensObj; ++i ) {
            if ( flame->fSpecies.FindSpecies( flame->fSensObj[i] ) >= 0 ) {
                flame->GetSpecies()->PrintProdRateTerms( flame->fSensObj[i], flame );
            }
        }
        if ( flame->fSensAnal ) {
            flame->SensitivityAnalysis( -1.0, -1.0, kPhysical );
        }

        if ( flame->GetArcLengthCont() ) {
            char	tail[28];
            int		izrefloc = flame->GetZRefLoc( currentGrid->GetX() );
            Double	tst = temp[izrefloc] + ( temp[izrefloc+1] - temp[izrefloc] )
                / (currentGrid->GetX()->vec[izrefloc+1] - currentGrid->GetX()->vec[izrefloc])
                * (flame->GetZRef() - currentGrid->GetX()->vec[izrefloc]);
            sprintf( tail, "Tst%04.0f", tst );
            bt->WriteOutput( object, NULL, tail );

            flame->IncNFlameletCount();
            if ( flame->GetDeltaArcLength() == 0.0 ) {
                flame->SetdlnChids( 1.0 );
                flame->SetdTds( 0.0 );
                flame->SetDeltaArcLength( flame->GetArcUp() ? flame->GetdlnChids() : -flame->GetdlnChids() );
            }
            else {
                Double dlnChids = ( log( flame->GetDissRate() ) - flame->GetChiContStart() )
                    / ( flame->GetDeltaChiref() * flame->GetDeltaArcLength() );
                Double dTds = ( temp[flame->GetZRefLoc( currentGrid->GetX() )] - flame->GetTempContStart() )
                    / ( flame->GetDeltaTref() * flame->GetDeltaArcLength() );
                flame->SetdlnChids( dlnChids );
                flame->SetdTds( dTds );
                Double ds = 0;
                Double absdlnChids = abs( dlnChids );
                Double absdTds = abs( dTds );
                ds = absdTds + absdlnChids;
                flame->SetDeltaArcLength( (flame->GetArcUp()) ? ds : -ds );
            }
            flame->SetTempContStart( temp[flame->GetZRefLoc( currentGrid->GetX() )] );
            flame->SetChiContStart( log( flame->GetDissRate() ) );
            if ( flame->GetNFlameletsCount() < flame->GetMaxFlamelets() ) {
                flame->fSolver->ReInit();
                fprintf( stderr, "\n\nFlamelet #%d: Tst is now %g chi_st is now %g\n", flame->GetNFlameletsCount()+1, flame->GetTempContStart(), exp(flame->GetChiContStart()) );
            }
            else{
                fprintf( stderr, "\n\n%d flamelets have been computed.\n", flame->GetMaxFlamelets() );
                fprintf( stderr, "The computation is stopped here. Increase the value of the\n" );
                fprintf( stderr, "variable fMaxFlamelets in file $FlameManSource/TCountDiffFlameMix.C\n" );
                fprintf( stderr, "if you want to compute more flamelets. You can also\n" );
                fprintf( stderr, "just start from one of the last solutions if it is not at a turning point in chi_st\n" );
            }
        }
        else {
            bt->WriteOutput( object, NULL, "" );
            if ( dissRate && ( dissRate->len < dissRate->phys_len - 1 || chiNotConv > 0.0 ) ) {
                flame->fSolver->ReInit();
                chiConv = flame->GetDissRate();
                if ( chiNotConv < 0.0 ) {
                    ++dissRate->len;
                }
                else {
                    dissRate->vec[dissRate->len] = chiNotConv;
                    chiNotConv = -1.0;
                }
                flame->SetChiContStart( log( dissRate->vec[dissRate->len] ) );
                fprintf( stderr, "%s%g\n", "scalar dissipation rate is now ", flame->GetDissRate()  );
            }
            else {
                flame->PostConvergence( object );
                CountDiffMixPostIter<Flame>( flame );
            }
        }
    }
    else {
        flame->RestoreSolution();
        if ( flame->GetArcLengthCont() ) {
            flame->SetDeltaArcLength( 0.5 * flame->GetDeltaArcLength() );

            if ( abs( flame->GetDeltaArcLength() ) > 0.1 * ( abs( flame->GetdlnChids() ) + abs(flame->GetdTds() ) ) ) {
                fprintf( stderr, "ds = %g\n", flame->GetDeltaArcLength() );
                flame->fSolver->ReInit();
            }
            else{
                flame->ReInitArcLengthCont();
                flame->PostConvergence( object );
                CountDiffMixPostIter<Flame>( flame );
            }
        }
        else if ( dissRate && dissRate->len < dissRate->phys_len ) {
            Double	interDissRate = chiConv + ( flame->GetDissRate() - chiConv ) * 0.25;
            if ( chiConv >= 0.0 && fabs( interDissRate - chiConv ) / chiConv >= 0.000001 ) {
                chiNotConv = dissRate->vec[dissRate->len];
                dissRate->vec[dissRate->len] = interDissRate;
                flame->SetChiContStart( log( dissRate->vec[dissRate->len] ) );
                flame->fSolver->ReInit();
                fprintf( stderr, "%s%g\n", "scalar dissipation rate is now ", flame->GetDissRate()  );
            }
        }
        else {
            flame->PostConvergence( object );
            CountDiffMixPostIter<Flame>( flame );
        }
    }
}

    template <typename Species>
void TCountDiffFlameMix<Species>::SaveSolution( void )
{
    int		k;
    int		len = fSavedTemp->len;
    Double	*lnChi = fSolLnChi->vec;
    Double	*saveLnChi = fSavedLnChi->vec;

    T1DFlame<Species>::SaveSolution();
    fSavedLnChi->len = fSolLnChi->len;

    for ( k = -1; k <= len; ++k ) {
        saveLnChi[k] = lnChi[k];
    }

    fSavedTempContStart = fTempContStart;
    fSavedLnChiContStart = fLnChiContStart;
}

    template <typename Species>
void TCountDiffFlameMix<Species>::RestoreSolution( void )
{
    int		len = fSavedTemp->len;
    int		k;
    Double	*lnChi = fSolLnChi->vec;
    Double	*saveLnChi = fSavedLnChi->vec;

    UpdateDimensions( len );

    T1DFlame<Species>::RestoreSolution();

    for ( k = -1; k <= len; ++k ) {
        lnChi[k] = saveLnChi[k];
    }

    fTempContStart = fSavedTempContStart;
    fLnChiContStart = fSavedLnChiContStart;

    SolutionToSolver();
}

    template <typename Species>
void TCountDiffFlameMix<Species>::SolutionToSolver( void )
{
    TNewtonPtr	bt = fSolver->bt;
    TGridPtr	grid = bt->GetGrid()->GetFine();
    int		nGridPoints = fSolLnChi->len;
    Double	*lnChi = fSolLnChi->vec;
    Double	**y = grid->GetY()->mat;

    T1DFlame<Species>::SolutionToSolver();

    for ( int k = 0; k < nGridPoints; ++k ) {
        y[k][fLnChi] = lnChi[k];
    }

    CountDiffMixPostIter<TCountDiffFlameMix<Species> >( this );
}

    template <typename Flame>
void CountDiffMixOutput( void *object, FILE *fp, const char* tail )
{
    Flame*	flame = ( Flame* )object;
    TNewtonPtr		bt = flame->GetSolver()->bt;
    T1DPropertiesPtr	props = flame->GetProperties();
    NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
    Double			*rho = props->GetDensity()->vec;
    Double			*mixMolarMass = props->GetMolarMass()->vec;
    Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
    TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
    Double			*x = currentGrid->GetX()->vec;
    Double			**massFracs = flame->GetMassFracs()->mat;
    Double			*temp = flame->GetTemperature()->vec;
    Double			**y = currentGrid->GetY()->mat;
    Double			*yLeft = currentGrid->GetYLeft()->vec,
                    *yRight = currentGrid->GetYRight()->vec;
    int				i, k;
    int				gridPoints = currentGrid->GetNGridPoints();
    int				nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
    int				nOfVariables = bt->GetNVariables();
    int				nOfEquations = bt->GetNEquations();
    int				firstSpecies = flame->GetOffsetFirstSpecies();
    int				tempOffset = flame->GetOffsetTemperature();
    time_t			theDate;
    char			buffer[80];
    ConstStringArray	varNames = flame->GetVariableNames();
    char			**names = flame->GetSpecies()->GetNames();
    Flag			fpOpen = FALSE;
    VectorPtr 		ZBilgerVec = NewVector( gridPoints + 2 );
    Double			*ZBilger = &ZBilgerVec->vec[kNext];

    if ( !fp ) {
        fpOpen = TRUE;
        fp = flame->GetOutputFile( NULL, tail, FileType::kNone );
    }

    ZBilger[kPrev] = bt->GetLeft();
    ZBilger[gridPoints] = bt->GetRight();
    for ( k = 0; k < gridPoints; ++k ) {
        ZBilger[k] = flame->ComputeZBilger( massFracs[k], massFracs[gridPoints], massFracs[-1] );
    }

    // write header
    fprintf( fp, "header\n\n" );

    fprintf( fp, "title = \"planar counterflow diffusion flame\"\n" );
    fprintf( fp, "mechanism = \"%s\"\n", flame->GetInputData()->fReactionFile );
    fprintf( fp, "author = \"%s\"\n", flame->GetAuthor() );
    time( &theDate );
    strcpy( buffer, ctime( &theDate ) );
    if ( buffer[strlen(buffer)-1] == '\n' )
        buffer[strlen(buffer)-1] = '\0';
    fprintf( fp, "date = \"%s\"\n\n", buffer );
    for ( i = 0; i < flame->GetNFuels(); ++i ) {
        fprintf( fp, "fuel = \"%s\"\n", varNames[firstSpecies+flame->GetFuelIndex( i )] );
    }
    fprintf( fp, "pressure = %g [bar]\n", flame->GetPressure() / 1.0e5 );
    fprintf( fp, "Z_st = %g [1/s]\n", flame->GetZRef() );
    fprintf( fp, "chi_st = %g [1/s]\n", flame->GetDissRate() );

    if ( flame->GetSpecies()->IsConstantLewisNumber() ) {
        fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
    }

    if ( flame->GetSolver()->bt->GetParameter() != 1.0 ) {
        fprintf( fp, "Parameter = %g\n", flame->GetSolver()->bt->GetParameter() );
    }

    fprintf( fp, "FlameLoc = %g\n", x[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
    fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );

    Double	EIFuel = 0.0;
    for ( i = 0; i < flame->GetNFuels(); ++i ) {
        EIFuel += flame->ComputeEmissionIndex( flame->GetFuelIndex( i ), x );
    }

    int	indNO = flame->GetSpecies()->FindSpecies( "NO" );
    if ( indNO > -1 ) {
        fprintf( fp, "EmissionIndexNO = %g [g/kg]\n"
                , -1000.0 * flame->ComputeEmissionIndex( indNO, x )
                / EIFuel );
    }

    int	indNO2 = flame->GetSpecies()->FindSpecies( "NO2" );
    if ( indNO2 > -1 ) {
        fprintf( fp, "EmissionIndexNO2 = %g [g/kg]\n"
                , -1000.0 * flame->ComputeEmissionIndex( indNO2, x )
                / EIFuel );
    }

    int	indN2O = flame->GetSpecies()->FindSpecies( "N2O" );
    if ( indN2O > -1 ) {
        fprintf( fp, "EmissionIndexN2O = %g [g/kg]\n"
                , -1000.0 * flame->ComputeEmissionIndex( indN2O, x )
                / EIFuel );
    }

    fprintf( fp, "\nFuelSide\n" );
    fprintf( fp, "begin\n" );
    fprintf( fp, "\tTemperature = %g [K]\n", yRight[tempOffset] );
    for ( i = 0; i < nOfSpecies; ++i ) {
        if ( fabs( massFracs[gridPoints][i] ) > 1.0e-20 ) {
            fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[gridPoints][i] );
        }
    }
    fprintf( fp, "end\n\n" );

    fprintf( fp, "OxidizerSide\n" );
    fprintf( fp, "begin\n" );
    fprintf( fp, "\tTemperature = %g [K]\n", yLeft[tempOffset] );
    for ( i = 0; i < nOfSpecies; ++i ) {
        if ( fabs( massFracs[kPrev][i] ) > 1.0e-20 ) {
            fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[kPrev][i] );
        }
    }
    fprintf( fp, "end\n\n" );

    fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
    fprintf( fp, "gridPoints = %d\n\n", gridPoints+2 );

    fprintf( fp, "body\n" );

    // write independent coordinate
    fprintf( fp, "Z\n" );
    fprintf( fp, "\t%-.6e", bt->GetLeft() );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", x[k] );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", bt->GetRight() );

    // write solution
    // write temperature
    flame->PrintFlameletVector( gridPoints+2, &temp[kPrev], "temperature [K]", fp );

    // write massfractions of species
    for ( i = 0; i < nOfSpecies; ++i ) {
        fprintf( fp, "massfraction-%s\n", names[i] );
        for ( k = 0; k < gridPoints+2; ++k ) {
            fprintf( fp, "\t%-.6e", massFracs[k-1][i] );
            if ( (k+1) % 5 == 0 ) {
                fprintf( fp, "\n" );
            }
        }
        if ( k % 5 ) {
            fprintf( fp, "\n" );
        }
    }

    if ( flame->fSoot ) {
        flame->GetSoot()->PrintFlameletFile( gridPoints, flame, fp );
    }

    if ( flame->fPrintMolarFractions ) {
        Double	locMolarMass;
        for ( i = 0; i < nOfSpecies; ++i ) {
            locMolarMass = molarMass[i];
            fprintf( fp, "molarfraction-%s\n", names[i] );
            fprintf( fp, "\t%-.6e", massFracs[kPrev][i] * mixMolarMass[-1] / locMolarMass );
            for ( k = 0; k < gridPoints; ++k ) {
                fprintf( fp, "\t%-.6e", massFracs[k][i] * mixMolarMass[k] / locMolarMass );
                if ( (k+2) % 5 == 0 ) {
                    fprintf( fp, "\n" );
                }
            }
            fprintf( fp, "\t%-.6e\n", massFracs[gridPoints][i] * mixMolarMass[gridPoints] / locMolarMass );
        }
    }

    flame->PrintFlameletVector( gridPoints+2, &mixMolarMass[kPrev], "W", fp );
    flame->PrintFlameletVector( gridPoints+2, &ZBilger[kPrev], "ZBilger", fp );


    //	write chi
    fprintf( fp, "chi [1/s]\n" );
    fprintf( fp, "\t%-.6e", flame->GetDissRate( bt->GetLeft(), rho[kPrev], exp( yLeft[flame->GetLnChi()] ) ) );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", flame->GetDissRate( x[k], rho[k], exp( y[k][flame->GetLnChi()] ) ) );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", flame->GetDissRate( bt->GetRight(), rho[gridPoints], exp( yRight[flame->GetLnChi()] ) ) );


    //	write density
    fprintf( fp, "density\n" );
    for ( k = 0; k < gridPoints+2; ++k ) {
        fprintf( fp, "\t%-.6e", rho[k-1] );
        if ( (k+1) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    //	if ( (k+1) % 5 ) {
    if ( k % 5 ) {
        fprintf( fp, "\n" );
    }

    //	write heat conductivity
    flame->PrintFlameletVector( gridPoints+2, &props->GetConductivity()->vec[kPrev], "lambda [W/m K]", fp );

    //	write heat capacity
    flame->PrintFlameletVector( gridPoints+2, &props->GetHeatCapacity()->vec[kPrev], "cp [J/kg K]", fp );

    //      write lambda/cp
    Double *lambda = props->GetConductivity()->vec;
    Double *cp = props->GetHeatCapacity()->vec;
    fprintf( fp, "lambdaOverCp [kg/ms]\n" );
    fprintf( fp, "\t%-.6e", lambda[kPrev]/cp[kPrev] );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", lambda[k]/cp[k] );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", lambda[gridPoints]/cp[gridPoints] );

    //  write viscosity
    flame->PrintFlameletVector( gridPoints+2, &props->GetViscosity()->vec[kPrev], "mu [kg/sm]", fp );

    //  write source term CO2
    int	indCO2 = flame->GetSpecies()->FindSpecies( "CO2" );
    fprintf( fp, "ProdRateCO2 [kg/m^3s]\n" );
    Double	**prod = flame->GetSpecies()->GetProductionRate()->mat;
    fprintf( fp, "\t%-.6e", 0.0 );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", prod[k][indCO2] );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", 0.0 );

    //  write source term H2O
    int indH2O = flame->GetSpecies()->FindSpecies( "H2O" );
    fprintf( fp, "ProdRateH2O [kg/m^3s]\n" );
    fprintf( fp, "\t%-.6e", 0.0 );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", prod[k][indH2O] );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", 0.0 );

    //  write source term CO
    int indCO = flame->GetSpecies()->FindSpecies( "CO" );
    fprintf( fp, "ProdRateCO [kg/m^3s]\n" );
    fprintf( fp, "\t%-.6e", 0.0 );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", prod[k][indCO] );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", 0.0 );

    //  write source term H2
    int indH2 = flame->GetSpecies()->FindSpecies( "H2" );
    fprintf( fp, "ProdRateH2 [kg/m^3s]\n" );
    fprintf( fp, "\t%-.6e", 0.0 );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", prod[k][indH2] );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", 0.0 );

    //  write ProgVar
    fprintf( fp, "ProgVar\n" );
    for ( k = 0; k < gridPoints+2; ++k ) {
        fprintf( fp, "\t%-.6e", massFracs[k-1][indCO2]+massFracs[k-1][indH2O]
                +massFracs[k-1][indCO]+massFracs[k-1][indH2] );
        if ( (k+1) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    if ( k % 5 ) {
        fprintf( fp, "\n" );
    }

    //  write source term ProdRateProgVar
    fprintf( fp, "ProdRateProgVar [kg/m^3s]\n" );
    fprintf( fp, "\t%-.6e", 0.0 );
    for ( k = 0; k < gridPoints; ++k ) {
        fprintf( fp, "\t%-.6e", prod[k][indCO2]+prod[k][indH2O]+prod[k][indCO]
                +prod[k][indH2] );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", 0.0 );

    //  write source term NO
    int indspec = flame->GetSpecies()->FindSpecies( "NO" );
    if ( indspec > -1 ) {
        flame->PrintProdRate( indspec, fp );
    }

    //	write enthalpy
    fprintf( fp, "TotalEnthalpy [J/kg]\n" );
    int		nSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
    Double	hTot = 0.0;
    Double	**ent = flame->GetSpecies()->GetEnthalpy()->mat;
    for ( k = 0; k < gridPoints+2; ++k ) {
        hTot = 0.0;
        for ( i = 0; i < nSpeciesIn; ++i ) {
            hTot += massFracs[k-1][i] * ent[k-1][i];
        }
        fprintf( fp, "\t%-.6e", hTot );
        if ( (k+1) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    if ( k % 5 ) {
        fprintf( fp, "\n" );
    }

    //	write heat release
    fprintf( fp, "HeatRelease [J/m^3 s]\n" );
    Double	heatrel = 0.0;

    fprintf( fp, "\t%-.6e", 0.0 );

    for ( k = 0; k < gridPoints; ++k ) {
        heatrel = 0.0;
        for ( i = 0; i < nSpeciesIn; ++i ) {
            heatrel += ent[k][i] * prod[k][i];
        }
        fprintf( fp, "\t%-.6e", -heatrel );
        if ( (k+2) % 5 == 0 ) {
            fprintf( fp, "\n" );
        }
    }
    fprintf( fp, "\t%-.6e\n", 0.0 );

    if ( flame->fProperties->GetRadiation() ) {
        //    write radiation source
        Double          *rad = props->GetRadiation()->GetRadiation()->vec;

        fprintf( fp, "RadiationSource [J/m^3 s]\n" );
        fprintf( fp, "\t%-.6e", 0.0 );

        for ( k = 0; k < gridPoints; ++k ) {
            fprintf( fp, "\t%-.6e", -rad[k] );

            // if ( (props->GetRadiation()->dqrF[k]+rad[k]) == 0){ cout << "equal" << endl;}
            // else{cout << "Nooooooooo" << endl;}

    	    if ( (k+2) % 5 == 0 ) {
    		    fprintf( fp, "\n" );
    	    }
        }

        fprintf( fp, "\t%-.6e\n", 0.0 );


        if ( flame->GetRadiationName() != "Thin" && flame->GetRadiationName() != "SNB"){
    		// kappa
    		fprintf( fp, "kappa [1/m]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );

    		for ( k = 0; k < gridPoints; ++k ) {
    			fprintf( fp, "\t%-.6e", props->GetRadiation()->kappa[k]);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}

            fprintf( fp, "\t%-.6e\n", 0.0 );

    		// Emission term
    		fprintf( fp, "Em [J/s/m^3]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );

    		for ( k = 0; k < gridPoints; ++k ) {
    			double em = props->GetRadiation()->kappa[k]*5.670367E-8*pow(temp[k],4)/M_PI;
    			fprintf( fp, "\t%-.6e", em);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );
        }

        if ( flame->GetRadiationName() == "WSGGJohansson" || flame->GetRadiationName() == "WSGGBordbar"){
        	// kappa1
        	fprintf( fp, "kappa1 [1/m]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			fprintf( fp, "\t%-.6e", props->GetRadiation()->kappa1[0][k]);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );
    		// Emission term1
    		fprintf( fp, "Em1 [J/s/m^3]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			double em = props->GetRadiation()->kappa1[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
    			fprintf( fp, "\t%-.6e", em);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );

        	// kappa2
        	fprintf( fp, "kappa2 [1/m]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			fprintf( fp, "\t%-.6e", props->GetRadiation()->kappa2[0][k]);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );
    		// Emission term2
    		fprintf( fp, "Em2 [J/s/m^3]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			double em = props->GetRadiation()->kappa2[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
    			fprintf( fp, "\t%-.6e", em);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );

        	// kappa3
        	fprintf( fp, "kappa3 [1/m]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			fprintf( fp, "\t%-.6e", props->GetRadiation()->kappa3[0][k]);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );
    		// Emission term3
    		fprintf( fp, "Em3 [J/s/m^3]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			double em = props->GetRadiation()->kappa3[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
    			fprintf( fp, "\t%-.6e", em);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );        

        	// kappa4
        	fprintf( fp, "kappa4 [1/m]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			fprintf( fp, "\t%-.6e", props->GetRadiation()->kappa4[0][k]);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );
    		// Emission term4
    		fprintf( fp, "Em4 [J/s/m^3]\n" );
    		fprintf( fp, "\t%-.6e", 0.0 );
    		for ( k = 0; k < gridPoints; ++k ) {
    			double em = props->GetRadiation()->kappa4[1][k]*5.670367E-8*pow(temp[k],4)/M_PI;
    			fprintf( fp, "\t%-.6e", em);
    			if ( (k+2) % 5 == 0 ) {
    				fprintf( fp, "\n" );
    			}
    		}
            fprintf( fp, "\t%-.6e\n", 0.0 );
        }

      //   else if( flame->GetRadiationName() == "RadiativeFrac"){

    		// Double          *rad = props->GetRadiation()->GetRadiation()->vec;
    		// Double          *kappa = props->GetRadiation()->GetKappa()->vec;

    		// // kappa
    		// fprintf( fp, "kappa [1/m]\n" );
    		// fprintf( fp, "\t%-.6e", 0.0 );

    		// for ( k = 0; k < gridPoints; ++k ) {
    		// 	fprintf( fp, "\t%-.6e", kappa[k]);
    		// 	if ( (k+2) % 5 == 0 ) {
    		// 		fprintf( fp, "\n" );
    		// 	}
    		// }
    		// fprintf( fp, "\t%-.6e\n", 0.0 );

    		// // Emission term
    		// fprintf( fp, "Em [J/s/m^3/sr]\n" );
    		// fprintf( fp, "\t%-.6e", 0.0 );

    		// for ( k = 0; k < gridPoints; ++k ) {
    		// 	double em = kappa[k]*5.670367E-8*pow(temp[k],4)/M_PI;
    		// 	fprintf( fp, "\t%-.6e", em);
    		// 	if ( (k+2) % 5 == 0 ) {
    		// 		fprintf( fp, "\n" );
    		// 	}
    		// }
    		// fprintf( fp, "\t%-.6e\n", 0.0 );
      //   }
    }

    if ( flame->fReactionFluxes ) {
	    //  write all reaction rates
	    T1DReactionPtr		reaction = flame->GetReaction();
	    int					nReactions = reaction->GetNOfReactions();
	    Double				**reactionRate = reaction->GetReactionRate()->mat;
	    char**				labels = reaction->GetLabels();

	    for ( i = 0; i < nReactions; ++i ) {
		    fprintf( fp, "ReactionRate_%s\n", labels[i] );
		    fprintf( fp, "\t%-.6e", 0.0 );
		    for ( k = 0; k < gridPoints; ++k ) {
			    fprintf( fp, "\t%-.6e", reactionRate[k][i] );
			    if ( (k+2) % 5 == 0 ) {
				    fprintf( fp, "\n" );
			    }
		    }
		    fprintf( fp, "\t%-.6e\n", 0.0 );
	    }
    }


    fprintf( fp, "trailer\n" );

    Double	*Le = flame->GetSpecies()->GetLewisNumber()->vec;
    for ( i = 0; i < nOfSpecies; ++i ) {
	    fprintf( fp, "%s\t%g\n", names[i], Le[i] );
    }
    if ( nOfEquations < nOfVariables) {
	    fprintf( fp, "number of converged equations is %d\n", nOfEquations );
    }

    if ( fpOpen ) {
	    fclose( fp );
    }
    DisposeVector( ZBilgerVec );
}

	template <typename Flame>
void CountDiffMixUpdateBoundary( void  * /*object*/ )
{
	/* nothing to do */
}

	template <typename Species>
void TCountDiffFlameMix<Species>::UpdateSolutionOnePoint( Double *y, int gridPoint )
{
	T1DFlame<Species>::UpdateSolution( y, gridPoint );

	fSolLnChi->vec[gridPoint] = y[fLnChi];
}


	template <typename Species>
void TCountDiffFlameMix<Species>::PrintProdRate( int speciesIndex, FILE *fp )
{
	int             k;
	TGridPtr        currentGrid = fSolver->bt->GetGrid()->GetCurrentGrid();
	int             fNGridPoints = currentGrid->GetNGridPoints();
	Double          *prodOfZ = New1DArray( fNGridPoints+2 );
	Double          *prodPlusOfZ = New1DArray( fNGridPoints+2 );
	Double          *prodMinusOfZ = New1DArray( fNGridPoints+2 );
	Double          *prodMinusOverYOfZ = New1DArray( fNGridPoints+2 );
	Double          **massFracs = GetMassFracs()->mat;
	Double          **reactionRate = fReaction.GetReactionRate()->mat;
	Double          **prodRate = fSpecies.GetProductionRate()->mat;
	char            buffer[128];
	sprintf( buffer, "ProdRate-%s [kg/m^3s]", fSpecies.GetNames()[speciesIndex] );

	prodOfZ[0] = prodOfZ[fNGridPoints+1] = 0.0;
	prodPlusOfZ[0] = prodPlusOfZ[fNGridPoints+1] = 0.0;
	prodMinusOfZ[0] = prodMinusOfZ[fNGridPoints+1] = 0.0;
	prodMinusOverYOfZ[0] = prodMinusOverYOfZ[fNGridPoints+1] = 0.0;
	for ( k = 0; k < fNGridPoints; ++k )
	{
		prodOfZ[k+1] = prodRate[k][speciesIndex];
		prodPlusOfZ[k+1] = fSpecies.GetPosNegProductionRate( speciesIndex, reactionRate[k], TRUE );
		prodMinusOfZ[k+1] = fSpecies.GetPosNegProductionRate( speciesIndex, reactionRate[k], FALSE );
		prodMinusOverYOfZ[k+1] = prodMinusOfZ[k+1] / MAX(1.0e-10, massFracs[k][speciesIndex] );
	}
	PrintFlameletVector( fNGridPoints+2, prodOfZ, buffer, fp );
	sprintf( buffer, "ProdRatePos-%s [kg/m^3s]", fSpecies.GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodPlusOfZ, buffer, fp );
	sprintf( buffer, "ProdRateNeg-%s [kg/m^3s]", fSpecies.GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodMinusOfZ, buffer, fp );
	sprintf( buffer, "ProdRateNegOverYNO-%s [kg/m^3s]", fSpecies.GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodMinusOverYOfZ, buffer, fp );

	Free1DArray( prodMinusOverYOfZ );
	Free1DArray( prodPlusOfZ );
	Free1DArray( prodMinusOfZ );
	Free1DArray( prodOfZ );
}

#endif // TCOUNT_DIFF_FLAME_MIX_HPP__
