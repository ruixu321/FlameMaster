#ifndef WSGGBordbar_H
#define WSGGBordbar_H

#include <string>
#include <iostream>
#include <vector>

#include "AbsorptionModel.h"

using namespace std;

class WSGGBordbar: public AbsorptionModel
{
    protected:

        vector<double> getKappa();
        vector<vector<double>> getKappa(int);


    private:

        int NBands = 5;
        
	//vector<vector<vector<double>>> c;
	//vector<vector<double>> d;

		vector<double> c;
		vector<double> d;

		vector <double> kappaOutput;

		vector<vector<double>> kappaOutput1;
		vector<vector<double>> kappaOutput2;
		vector<vector<double>> kappaOutput3;
		vector<vector<double>> kappaOutput4;

		
    public:

	// compute the coefficients for WSGGBordbar model

	double ComputeCoeffs() {return 0;};
	double ComputeCoeffs( int n, double T);

	// read WSGGBordbar parameters
	void ReadData();

	// return coeffs

	double transmissivity(const vector <double>& T, const vector <double>& P, const vector <double>& XH2O, const vector <double>& XCO2, const vector <double>& delta, int iBand, double mu, int i, int j){ return 0;};


	double kappa(int iBand, double MolarFracH2O, double MolarFracCO2);
	double kappa(int iBand, double T){ return 0;};
	double kappa(double T, double pH2O, double pCO2){return 0;};
	// vector<double> getKappa() {return vector<double> (0,0.);};

	double a();
	double a(int n, double T); // return a coeff for band n at temperature T
	double a(int n, double T, double MR); // return a coeff for band n at temperature T


	double b(int i, int j, double MR);

	vector<vector<double>> solve(vector<double> mu, int Npoints, vector <double> T, vector <double> P, vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign);


	// wavenumber
	double nu(int iBand) {return 0.;};

	// return number of bands

	int NbBands();

	// write data

	WSGGBordbar();

};

#endif
