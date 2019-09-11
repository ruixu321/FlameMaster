#include "TProperties.h"
#include "RadiationModel.h"
#include <string>
#include <vector>

// multiplies soot source term in energy equation, used to adjust number of optically thin directions
#define RADCALCONST 1.0 

TProperties::~TProperties( void )
{
}

void TProperties::ComputeMixtureMolarMass( Double &mixMolarMass, Double *Y, Double *molarMass, int nSpeciesInSystem )
{
	int i;	
   
   	for ( i = 0, mixMolarMass = 0.0; i < nSpeciesInSystem; ++i ) {
		mixMolarMass += Y[i] / molarMass[i];
	}

	mixMolarMass = 1.0 / mixMolarMass;
}

void T0DProperties::InitT0DProperties( TInputDataPtr input )
{
	fRadiationName = input->fRadiationName;

	if ( input->fPressure ) {
		fPressure = input->fPressure->vec[0];
/*		if ( input->fPressure->len > 1 ) {
			cerr << "#warning: more than one pressure specified, use" 
					<< fPressure << NEWL;
		}*/
		cerr << "initial pressure is " << fPressure / 1.0e5 << " bar" << NEWL;
	}
	else {
		fPressure = -1.0;
//		cerr << "#error: no pressure specified" << NEWL;
//		exit( 2 );
	}
	fRadiation = NULL;
	if ( input->fWithRadiation ) {
 	  	if ( !( fRadiation = new T0DRadiation( input ) ) ) {
			FatalError( "memory allocation of TRadiation failed" );
		}
	}
}

T0DProperties::~T0DProperties( void )
{
	if ( fRadiation ) {
		delete fRadiation;
	}
}

void TRadiation::InitRadiation( TInputDataPtr input )
{
	fH2OIndex = input->fH2OIndex;
	fCO2Index = input->fCO2Index;
	fCH4Index = input->FindSpecies( "CH4" );
	fCOIndex = input->FindSpecies( "CO" );
}

TRadiation::~TRadiation( void )
{
    //delete WSGG;
}

void TRadiation::ComputeRadiationOnePoint( Double *radiation, Double temp, Double *Y, Double *molarMass, Double density, 
	Double HeatRelease, Double RadiativeFrac, Double *Kappa)
{
	*radiation = -RadiativeFrac*HeatRelease;

	//cout << "Calculation of Radiation similar to RadCal" << endl;

	Double			alphaCO2, alphaH2O, alphaCH4, alphaCO;
	Double			C_CO2 = ( fCO2Index == -1 ) ? 0.0 : Y[fCO2Index] / molarMass[fCO2Index];
	Double			C_H2O = ( fH2OIndex == -1 ) ? 0.0 : Y[fH2OIndex] / molarMass[fH2OIndex];
	Double			C_CH4 = ( fCH4Index == -1 ) ? 0.0 : Y[fCH4Index] / molarMass[fCH4Index];
	Double			C_CO = ( fCOIndex == -1 ) ? 0.0 : Y[fCOIndex] / molarMass[fCOIndex];
	Double			Tm1 = 1000.0/temp;
	Double			T2 = temp * temp;
	Double			c0, c1, c2, c3, c4, c5;
	Double			atmm1ToPam1 = 1.0 / 101330.0; // 1/atm to 1/Pa

	// radcal
	/* corrected version 06/15/00*/	
	// first H2O
	c0 = -0.23093;
	c1 = -1.12390;
	c2 =  9.41530;
	c3 = -2.99880;
	c4 =  0.51382;
	c5 = -1.86840E-05;
	alphaH2O = c0+Tm1*(c1+Tm1*(c2+Tm1*(c3+Tm1*(c4+Tm1*c5)))); // (atm m)^-1

	// then CO2
	c0 =  18.741;
	c1 = -121.310;
	c2 =  273.500;
	c3 = -194.050;
	c4 =  56.310;
	c5 = -5.8169;
	alphaCO2 = c0+Tm1*(c1+Tm1*(c2+Tm1*(c3+Tm1*(c4+Tm1*c5)))); // (atm m)^-1


	// then CH4
	alphaCH4 = 6.6334 - 0.0035686*temp + 1.6682e-08*T2 + 2.5611e-10*T2*temp - 2.6558e-14*T2*T2; // (atm m)^-1

	// and CO
	if ( temp < 750.0 ) {
		c0 = 4.7869;
		c1 = -0.06953;
		c2 = 2.95775e-4;
		c3 = -4.25732e-7;
		c4 = 2.02894e-10;
	}
	else{
		c0 = 10.09;
		c1 = -0.01183;
		c2 = 4.7753e-6;
		c3 = -5.87209e-10;
		c4 = -2.5334e-14;
	}
	alphaCO =  c0 + temp * ( c1 + temp * ( c2 + temp * ( c3 + temp * c4) ) ); // (atm m)^-1

	*Kappa = RGAS * density * temp * atmm1ToPam1 * 
		( alphaCO2 * C_CO2 + alphaH2O * C_H2O + alphaCH4 * C_CH4 + alphaCO * C_CO );

}


void TRadiation::ComputeRadiationOnePoint( Double *radiation, double dqr)
{
	*radiation = -dqr;
}

void TRadiation::ComputeRadiationRadiation(Double *temp, vector<double> XH2OF, vector<double>XCO2F, int nSpeciesInSystem,
                                                  vector<double> physicalGrid, double Delta, double chiRef, 
                                                  vector<double> heatrel, string radiationName) 
{
	// declaration 

	vector <double> T(physicalGrid.size()+3, 200.);
	vector <double> XH2O(physicalGrid.size()+1,0.);
	vector <double> XCO2(physicalGrid.size()+1,0.);
	vector <double> HRR(physicalGrid.size()+1,0.);
	vector <double> P(physicalGrid.size()+1,1.);
	vector <double> Grid(physicalGrid.size()+2,0.);

	// interpolation function

	auto interpolation = [](double x, double xa, double xb, double ya, double yb){
		return ya + (x - xa)*(yb-ya)/(xb-xa);
	};

	// calculate the molar fraction

	// oxy boundary

	T[0] = temp[-1]; //ok because of size of nGridPoints
	T[1] = temp[-1]; //ok because of size of nGridPoints


	XCO2[0] = XCO2F[0]; //based on physicalGrid
	XH2O[0] = XH2OF[0];


	for (int i = 0; i<physicalGrid.size();i++){
		Grid[i+1] = physicalGrid[i];
	}

	// add two points

	Grid[0] = 0.;

	Grid[Grid.size()-1] = physicalGrid[physicalGrid.size()-1] + Delta;

	// fuel

	T[T.size()-1] = temp[physicalGrid.size()-3];
	T[T.size()-2] = temp[physicalGrid.size()-3];


	XCO2[XCO2.size()-1] = XCO2F[XCO2F.size()-3];
	XH2O[XH2O.size()-1] = XH2OF[XCO2F.size()-3];

	HRR[HRR.size()-1] = heatrel[heatrel.size()-1];

	for (int i = 0; i<physicalGrid.size()-1; ++i){


		double x = (physicalGrid[i+1] - physicalGrid[i])/2. + physicalGrid[i];
		XCO2[i+1] = interpolation(x, physicalGrid[i], physicalGrid[i+1], XCO2F[i], XCO2F[i+1]);
		XH2O[i+1] = interpolation(x, physicalGrid[i], physicalGrid[i+1], XH2OF[i], XH2OF[i+1]);
		HRR[i+1] = interpolation(x, physicalGrid[i], physicalGrid[i+1], heatrel[i], heatrel[i+1]);
		T[i+2] = interpolation(x, physicalGrid[i], physicalGrid[i+1], temp[i-1], temp[i]);
	}

	T[T.size()-3] = temp[physicalGrid.size()-3];

	string chiString = std::to_string(chiRef);
	string filedqr =  "./res/dqr_" + radiationName + "_" + chiString + ".res";
	string fileqr =  "./res/qr_" + radiationName + "_" + chiString + ".res";
	string model = radiationName;


	RadiationModel Rad(model, Grid, T, XH2O, XCO2, P, filedqr, fileqr, HRR);

	// cout << model << endl;

	Rad.Solve();

	// cout << "After  Rad.Solve" << endl;

	// reconstruction of radiative source term
	vector<double> dqrOutput = Rad.GetRadSource();

	vector<double> kappaOutput;

	vector<vector<double>> kappaOutput1;
	vector<vector<double>> kappaOutput2;
	vector<vector<double>> kappaOutput3;
	vector<vector<double>> kappaOutput4;

	if ( model != "SNB" && model != "Thin") {
		
		kappaOutput = Rad.GetKappa();

		if (model == "WSGGJohansson" || model == "WSGGBordbar"){
			kappaOutput1 = Rad.GetKappa(1);
			kappaOutput2 = Rad.GetKappa(2);
			kappaOutput3 = Rad.GetKappa(3);
			kappaOutput4 = Rad.GetKappa(4);
		}

	}
	// interpolation of dqr from WSGG -> FlameMaster

	dqrF = vector <double> (physicalGrid.size()-2,0);

	if ( model != "SNB" && model != "Thin") {

		kappa = vector <double> (physicalGrid.size()-2,0);
	
		if (model == "WSGGJohansson" || model == "WSGGBordbar") {

			kappa1 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );
			kappa2 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );
			kappa3 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );
			kappa4 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );

		}
	}

	for (int i =0; i<dqrF.size(); ++i){

		double xm = (Grid[i+2] - Grid[i+1])/2. + Grid[i+1]; //ok
		double xp = (Grid[i+3] - Grid[i+2])/2. + Grid[i+2]; //ok

		dqrF[i] = interpolation(physicalGrid[i+1], xm, xp, dqrOutput[i+1], dqrOutput[i+2]);

		if ( model != "SNB" && model != "Thin") {

			kappa[i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput[i+1], kappaOutput[i+2]);

			if ( model == "WSGGJohansson" || model == "WSGGBordbar" ){    

				kappa1[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput1[0][i+1], kappaOutput1[0][i+2]);
				kappa1[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput1[1][i+1], kappaOutput1[1][i+2]);
				kappa2[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput2[0][i+1], kappaOutput2[0][i+2]);
				kappa2[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput2[1][i+1], kappaOutput2[1][i+2]);
				kappa3[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput3[0][i+1], kappaOutput3[0][i+2]);
				kappa3[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput3[1][i+1], kappaOutput3[1][i+2]);
				kappa4[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput4[0][i+1], kappaOutput4[0][i+2]);
				kappa4[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput4[1][i+1], kappaOutput4[1][i+2]);

			}
		}
		//dqrF[i] = 0.; 
		//cout << physicalGrid[i+1] << " " << dqrF[i] <<" "<< xm << " " << xp << " " << dqrOutput[i+1] << " " << dqrOutput[i+2] << " " << temp[i] << " " << dqrOutput.size() << " " << physicalGrid.size() <<  endl;
	}
}

void TRadiation::ComputeRadiationUnsteady(vector<double> temp, vector<double> XH2OF, vector<double>XCO2F, int nSpeciesInSystem,
                                                  vector<double> physicalGrid, double Delta, double chiRef, vector<double> heatrel, string radiationName) 
{
	// declaration 

	vector <double> T(physicalGrid.size()+3, 200.);
	vector <double> XH2O(physicalGrid.size()+1,0.);
	vector <double> XCO2(physicalGrid.size()+1,0.);
	vector <double> HRR(physicalGrid.size()+1,0.);
	vector <double> P(physicalGrid.size()+1,1.);
	vector <double> Grid(physicalGrid.size()+2,0.);

	// interpolation function

	auto interpolation = [](double x, double xa, double xb, double ya, double yb)
		{ return ya + (x - xa)*(yb-ya)/(xb-xa); };

	// calculate the molar fraction

	// oxy boundary

	// T[0] = temp[-1]; //ok because of size of nGridPoints
	T[0] = temp[0];

	// T[1] = temp[-1]; //ok because of size of nGridPoints
	T[1] = temp[0];

	XCO2[0] = XCO2F[0]; //based on physicalGrid
	XH2O[0] = XH2OF[0];


	for (int i = 0; i<physicalGrid.size();i++){
		Grid[i+1] = physicalGrid[i];
	}

	// add two points

	Grid[0] = 0.;

	Grid[Grid.size()-1] = physicalGrid[physicalGrid.size()-1] + Delta;

	// fuel

	T[T.size()-1] = temp[physicalGrid.size()-3];
	T[T.size()-2] = temp[physicalGrid.size()-3];


	XCO2[XCO2.size()-1] = XCO2F[XCO2F.size()-3];
	XH2O[XH2O.size()-1] = XH2OF[XCO2F.size()-3];

	HRR[HRR.size()-1] = heatrel[heatrel.size()-1];

	for (int i = 0; i<physicalGrid.size()-1; ++i){

		double x = (physicalGrid[i+1] - physicalGrid[i])/2. + physicalGrid[i];
		XCO2[i+1] = interpolation(x, physicalGrid[i], physicalGrid[i+1], XCO2F[i], XCO2F[i+1]);
		XH2O[i+1] = interpolation(x, physicalGrid[i], physicalGrid[i+1], XH2OF[i], XH2OF[i+1]);
		HRR[i+1] = interpolation(x, physicalGrid[i], physicalGrid[i+1], heatrel[i], heatrel[i+1]);
		// T[i+2] = interpolation(x, physicalGrid[i], physicalGrid[i+1], temp[i-1], temp[i]);
		T[i+2] = interpolation(x, physicalGrid[i], physicalGrid[i+1], temp[i], temp[i+1]);
	}

	T[T.size()-3] = temp[physicalGrid.size()-3];

	string chiString = std::to_string(chiRef);
	string filedqr =  "./res/dqr_" + radiationName + "_" + chiString + ".res";
	string fileqr =  "./res/qr_" + radiationName + "_" + chiString + ".res";
	string model = radiationName;

	RadiationModel Rad(model, Grid, T, XH2O, XCO2, P, filedqr, fileqr, HRR);

	Rad.Solve();
	// reconstruction of radiative source term

	vector<double> dqrOutput = Rad.GetRadSource();

	vector<double> kappaOutput;

	vector<vector<double>> kappaOutput1;
	vector<vector<double>> kappaOutput2;
	vector<vector<double>> kappaOutput3;
	vector<vector<double>> kappaOutput4;

	if ( model != "SNB" && model != "Thin") {
		
		kappaOutput = Rad.GetKappa();

		if (model == "WSGGJohansson" || model == "WSGGBordbar"){
			kappaOutput1 = Rad.GetKappa(1);
			kappaOutput2 = Rad.GetKappa(2);
			kappaOutput3 = Rad.GetKappa(3);
			kappaOutput4 = Rad.GetKappa(4);
		}

	}

	// interpolation of dqr from WSGG -> FlameMaster

	dqrF = vector <double> (physicalGrid.size()-2,0);

	if ( model != "SNB" && model != "Thin") {

		kappa = vector <double> (physicalGrid.size()-2,0);
	
		if (model == "WSGGJohansson" || model == "WSGGBordbar") {
			kappa1 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );
			kappa2 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );
			kappa3 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );
			kappa4 = vector<vector<double>> (2, vector <double> (physicalGrid.size()-2,0) );
		}
	}
	
	for (int i =0; i<dqrF.size(); ++i){

		double xm = (Grid[i+2] - Grid[i+1])/2. + Grid[i+1]; //ok
		double xp = (Grid[i+3] - Grid[i+2])/2. + Grid[i+2]; //ok

		dqrF[i] = interpolation(physicalGrid[i+1], xm, xp, dqrOutput[i+1], dqrOutput[i+2]);

		if ( model != "SNB" && model != "Thin") {

			kappa[i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput[i+1], kappaOutput[i+2]);

			if ( model == "WSGGJohansson" || model == "WSGGBordbar" ){    
				kappa1[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput1[0][i+1], kappaOutput1[0][i+2]);
				kappa1[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput1[1][i+1], kappaOutput1[1][i+2]);
				kappa2[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput2[0][i+1], kappaOutput2[0][i+2]);
				kappa2[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput2[1][i+1], kappaOutput2[1][i+2]);
				kappa3[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput3[0][i+1], kappaOutput3[0][i+2]);
				kappa3[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput3[1][i+1], kappaOutput3[1][i+2]);
				kappa4[0][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput4[0][i+1], kappaOutput4[0][i+2]);
				kappa4[1][i] = interpolation(physicalGrid[i+1], xm, xp, kappaOutput4[1][i+1], kappaOutput4[1][i+2]);
			}
		}
		//dqrF[i] = 0.; 
		//cout << physicalGrid[i+1] << " " << dqrF[i] <<" "<< xm << " " << xp << " " << dqrOutput[i+1] << " " << dqrOutput[i+2] << " " << temp[i] << " " << dqrOutput.size() << " " << physicalGrid.size() <<  endl;
	}
}


void TRadiation::ComputeRadiationOnePoint( Double *radiation, Double temp, Double *Y, Double *molarMass, Double density )
{
	//cout << "Calculation of Radiation similar to RadCal" << endl;
	static Double	constant = RADCALCONST * 4.0 * STEFBOLTZ * RGAS;

	Double			alphaCO2, alphaH2O, alphaCH4, alphaCO;
	Double			C_CO2 = ( fCO2Index == -1 ) ? 0.0 : Y[fCO2Index] / molarMass[fCO2Index];
	Double			C_H2O = ( fH2OIndex == -1 ) ? 0.0 : Y[fH2OIndex] / molarMass[fH2OIndex];
	Double			C_CH4 = ( fCH4Index == -1 ) ? 0.0 : Y[fCH4Index] / molarMass[fCH4Index];
	Double			C_CO = ( fCOIndex == -1 ) ? 0.0 : Y[fCOIndex] / molarMass[fCOIndex];
	Double			Tm1 = 1000.0/temp;
	Double			T2 = temp * temp;
	Double			c0, c1, c2, c3, c4, c5;
	Double			atmm1ToPam1 = 1.0 / 101330.0; // 1/atm to 1/Pa

	// radcal
	/* corrected version 06/15/00*/	
	// first H2O
	c0 = -0.23093;
	c1 = -1.12390;
	c2 =  9.41530;
	c3 = -2.99880;
	c4 =  0.51382;
	c5 = -1.86840E-05;
	alphaH2O = c0+Tm1*(c1+Tm1*(c2+Tm1*(c3+Tm1*(c4+Tm1*c5)))); // (atm m)^-1

	// then CO2
	c0 =  18.741;
	c1 = -121.310;
	c2 =  273.500;
	c3 = -194.050;
	c4 =  56.310;
	c5 = -5.8169;
	alphaCO2 = c0+Tm1*(c1+Tm1*(c2+Tm1*(c3+Tm1*(c4+Tm1*c5)))); // (atm m)^-1


	// then CH4
	alphaCH4 = 6.6334 - 0.0035686*temp + 1.6682e-08*T2 + 2.5611e-10*T2*temp - 2.6558e-14*T2*T2; // (atm m)^-1

	// and CO
	if ( temp < 750.0 ) {
		c0 = 4.7869;
		c1 = -0.06953;
		c2 = 2.95775e-4;
		c3 = -4.25732e-7;
		c4 = 2.02894e-10;
	}
	else{
		c0 = 10.09;
		c1 = -0.01183;
		c2 = 4.7753e-6;
		c3 = -5.87209e-10;
		c4 = -2.5334e-14;
	}
	alphaCO =  c0 + temp * ( c1 + temp * ( c2 + temp * ( c3 + temp * c4) ) ); // (atm m)^-1

	//	*radiation = -constant * density * temp * ( pow( temp, 4.0 ) - pow( 300.0, 4.0 ) )
	*radiation = -constant * density * pow( temp, 5.0 ) 
		* atmm1ToPam1 * ( alphaCO2 * C_CO2 + alphaH2O * C_H2O + alphaCH4 * C_CH4 + alphaCO * C_CO );
}


#ifndef ZEROD
void T1DRadiation::InitT1DRadiation( TInputDataPtr input )
{
	fRadiation = NewVector( input->fMaxGridPoints );
	fKappa = NewVector( input->fMaxGridPoints );
}

T1DRadiation::~T1DRadiation( void )
{
	DisposeVector( fRadiation );
}

void T1DProperties::InitT1DProperties( TInputDataPtr input )
{
	int	maxGridPoints = input->fMaxGridPoints;

	fViscosity = NewVector( maxGridPoints+2 );
	fViscosity->vec = &fViscosity->vec[kNext];
	fConductivity = NewVector( maxGridPoints+2 );
	fConductivity->vec = &fConductivity->vec[kNext];
	fDensity = NewVector( maxGridPoints+2 );
	fDensity->vec = &fDensity->vec[kNext];
	fHeatCapacity = NewVector( maxGridPoints+2 );
	fHeatCapacity->vec = &fHeatCapacity->vec[kNext];
	fMolarMass = NewVector( maxGridPoints+2 );
	fMolarMass->vec = &fMolarMass->vec[kNext];
	fRadiation = NULL;
	if ( input->fWithRadiation ) {
		if ( !( fRadiation = new T1DRadiation( input ) ) ) FatalError( "memory allocation of TRadiation failed" );
	}
	fdCpdT = NewVector( maxGridPoints );
}

T1DProperties::~T1DProperties( void )
{
	DisposeVector( fdCpdT );
	if ( fRadiation ) {
		delete fRadiation;
	}
	fMolarMass->vec = &fMolarMass->vec[kPrev];
	DisposeVector( fMolarMass );
	fHeatCapacity->vec = &fHeatCapacity->vec[kPrev];
	DisposeVector( fHeatCapacity );
	fDensity->vec = &fDensity->vec[kPrev];
	DisposeVector( fDensity );
	fConductivity->vec = &fConductivity->vec[kPrev];
	DisposeVector( fConductivity );
	fViscosity->vec = &fViscosity->vec[kPrev];
	DisposeVector( fViscosity );
}

void T1DProperties::PrintProperties( int k )
{
	static Flag 	init = FALSE;
	FILE		*fp = NULL;

	if ( !init ) {
		if ( !( fp = fopen( "properties.tout", "w" ) ) ) { 
			cerr << "#warning: unable to open file 'properties.tout'" << NEWL;
			return;
		}
		init = TRUE;
	}
	else {
		if ( !( fp = fopen( "properties.tout", "a" ) ) ) { 
			cerr << "#warning: unable to open file 'properties.tout'" << NEWL;
			return;
		}
	}
	fprintf( fp, "Viscosity = %g\n", fViscosity->vec[k] );
	fprintf( fp, "Density = %g\n", fDensity->vec[k] );
	fprintf( fp, "HeatCapacity = %g\n", fHeatCapacity->vec[k] );
	fprintf( fp, "Conductivity = %g\n", fConductivity->vec[k] );
	fprintf( fp, "fMolarMass = %g\n", fMolarMass->vec[k] );

	fprintf( fp, "\n\n\n");
}

void T1DProperties::PrintProperties( TNewtonPtr bt )
{
	char		fName[32];
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	Double		*mu = fViscosity->vec;
	Double		*rho = fDensity->vec;
	Double		*cp = fHeatCapacity->vec;
	Double		*lambda = fConductivity->vec;
	Double		*molarMass = fMolarMass->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	static int	counter = 0;

	sprintf( fName, "props%d.dout", ++counter );
	if ( counter >= 10 ) counter = 0;
	if ( !( fp = fopen( fName, "w" ) ) ) { 
		cerr << "#warning: unable to open file " << fName << NEWL;
		return;
	}
	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n", "eta", "mu", "rho", "cp", "lambda", "M" );
	fprintf( fp, "%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\n", left, mu[-1], rho[-1], cp[-1], lambda[-1], molarMass[-1] );
	for ( int i = 0; i < nGridPoints; ++i ) {
		fprintf( fp, "%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\n", x[i], mu[i], rho[i], cp[i], lambda[i], molarMass[i] );
	}
	fprintf( fp, "%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\n", right, mu[nGridPoints], rho[nGridPoints], cp[nGridPoints], lambda[nGridPoints], molarMass[nGridPoints] );
	fclose( fp );
}
#endif // ZEROD
