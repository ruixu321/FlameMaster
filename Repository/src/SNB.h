#ifndef SNB_H
#define SNB_H

#include <string>
#include <iostream>
#include <vector>

#include "AbsorptionModel.h"

using namespace std;

#define H2OSize1 449
#define H2OSize2 44002
#define Ncols 49
#define CO2Size1 323
#define CO2Size2 31654



class SNB: public AbsorptionModel
{
    private:

        int NBands = 449;

        double dataH2O [H2OSize2];
        double dataCO2 [CO2Size2];

    public:

        // compute the coefficients for SNB model

        double ComputeCoeffs() {return 0;};
        double ComputeCoeffs( int n, double T);

        // read SNB parameters
        void ReadData();
        void  readCO2();
        void  readH2O();
        //tuple<double, int> dataInterpolation(double T);
        void dataInterpolation(double T, int &IT, double &RT);

        // return coeffs

        double transmissivity(const vector <double>& T, const vector <double>& P, const vector <double>& XH2O, const vector <double>& XCO2, const vector <double>& delta, int iBand, double mu, int i, int j);
        
        double kappa(int iBand, double pH2O, double pCO2);
        double kappa(double T, double pH2O, double pCO2){return 0;};
        double kappa(int iBand, double T){ return 0;};

        vector<double> getKappa() {return vector<double> (0,0.);};
        vector <vector<double>> getKappa(int i) {return vector<vector<double>> (0, vector <double> (0, 0.) );};
        
        double a();
        double a(int n, double T); // return a coeff for band n at temperature T

        vector<vector<double>> solve(vector<double> mu, int Npoints, vector <double> T, vector <double> P, vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign);
    

        // wavenumber
        double nu(int iBand);

        // return number of bands

        int NbBands();

        // constructor 

        SNB();

};

#endif
