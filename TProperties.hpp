#ifndef TPROPERTIES_HPP__
#define TPROPERTIES_HPP__

template<class Species>
void TProperties::InitProperties( TInputDataPtr /*input*/, Species* species )
{
	int nSpecies = species->GetNOfSpecies();
}

template<class Species>
void TProperties::ComputeDCpDT( Double &dCpdT, Double *Y, Double temp, Species* species )
{
	int 	nSpeciesInSystem = species->GetNSpeciesInSystem(); 
	
	dCpdT = 0;
	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
		dCpdT += Y[i] * species->ComputeDCpiDT( i, temp );
	}
}

template<class Species>
void T0DProperties::CompMixtureProps( Double *heatCapacity, Double *conductivity, Double *mu, Double *Y, Double temp
			, Double &pressure, Double &density, EqOfState what, int nSpeciesIn, Species* species )
{
//  calculate c_p and either pressure or density, depending on 'EqOfState what'

	int 	i;
	Double	*deltaI = species->GetDeltaI()->vec;
	
// ATTENTION hp
// the following is a very dirty fix of the problem that
// in homogeneous combustion transport data is not needed, hence deltaI
// not computed. In 1D unsteady, however deltaI is computed and transport
// data is needed
	if ( deltaI[0] ) {
		fMixHeatCapacity = Y[0] * heatCapacity[0];
		fMixConductivity = Y[0] / deltaI[0] * conductivity[0];
		fMixViscosity = Y[0] / deltaI[0] * mu[0];
		for ( i = 1; i < nSpeciesIn; ++i ){
			fMixHeatCapacity += Y[i] * heatCapacity[i];
			fMixConductivity += Y[i] / deltaI[i] * conductivity[i];
			fMixViscosity += Y[i] / deltaI[i] * mu[i];
		}
	}
	else {
		fMixHeatCapacity = Y[0] * heatCapacity[0];
		for ( i = 1; i < nSpeciesIn; ++i ){
			fMixHeatCapacity += Y[i] * heatCapacity[i];
		}
	}

	switch ( what ) {
		case kDensFromPress:
			density = pressure * fMixMolarMass / ( RGAS * temp );
			break;
		case kPressFromDens:
			pressure = density * RGAS * temp / fMixMolarMass;
			break;
		default:
			cerr << "#error in function T0DProperties::CompMixtureProps" << NEWL;
			exit(2);
	}
//	fprintf( stderr, "cp = %g\trho = %g\n", fMixHeatCapacity, density );
}

template<class Species>
void T1DProperties::CompMixtureProps( TFlameNodePtr flameNode, Double *Y, Double temp
					, Double pressure, Species* species )
{
//  calculate c_p, lambda and mu of the mixture

	int 	i;
	int		nSpeciesInSystem = species->GetNSpeciesInSystem();
	Double	y_over_delta = 0.0;
	Double	*c_p = flameNode->heatCapacity;
	Double	*lambda = flameNode->conductivity;
	Double	*mu = flameNode->viscosity;
	Double	*molarMass = species->GetMolarMass()->vec;
	Double	&mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	&mixConductivity = *flameNode->mixConductivity;
	Double	&mixViscosity = *flameNode->mixViscosity;
	Double	&mixDensity = *flameNode->mixDensity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*deltaI = flameNode->deltaI;
	

	mixHeatCapacity = Y[0] * c_p[0];
	
//	y_over_delta = Y[0] / CompDeltaI( 0, nSpeciesInSystem, molarMass, mu, Y );
	y_over_delta = Y[0] / deltaI[0];

	mixConductivity = y_over_delta * lambda[0];
	mixViscosity = y_over_delta * mu[0];

	for ( i = 1; i < nSpeciesInSystem; ++i ){
		mixHeatCapacity += Y[i] * c_p[i];
	
//		y_over_delta = Y[i] / CompDeltaI( i, nSpeciesInSystem, molarMass, mu, Y );
		y_over_delta = Y[i] / deltaI[i];
		mixConductivity += y_over_delta * lambda[i];
		mixViscosity += y_over_delta * mu[i];
	}
	mixDensity = pressure * mixMolarMass / ( RGAS * temp );
}

#endif // TPROPERTIES_HPP__

