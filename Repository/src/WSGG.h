#ifndef WSGG_H
#define WSGG_H

#include <string>
#include <iostream>
#include <vector>

#include "AbsorptionModel.h"

using namespace std;

class WSGG: public AbsorptionModel
{
    protected:

        vector<double> getKappa();
        vector < vector<double> > getKappa(int);
    private:
        int NBands = 25;
        
        vector<vector<double>> b;

        vector <double> kappaOutput;

    public:

        // compute the coefficients for WSGG model

        double ComputeCoeffs() {return 0;};
        double ComputeCoeffs( int n, double T);

        // read WSGG parameters
        void ReadData();

        // return coeffs

        double transmissivity(const vector <double>& T, const vector <double>& P, const vector <double>& XH2O, const vector <double>& XCO2, const vector <double>& delta, int iBand, double mu, int i, int j){ return 0;};
        
        
        double kappa(int iBand, double pH2O, double pCO2);
        double kappa(int iBand, double T){ return 0;};
        double kappa(double T, double pH2O, double pCO2){return 0;};

        // vector<double> getKappa() {return vector<double> (0,0.);};

        double a();
        double a(int n, double T); // return a coeff for band n at temperature T


        vector<vector<double>> solve(vector<double> mu, int Npoints, vector <double> T, vector <double> P, vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign);
    

        // wavenumber
        double nu(int iBand) {return 0.;};

        // return number of bands

        int NbBands();

        // write data

        WSGG();

};

#endif
