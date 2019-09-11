#ifndef ABSORPTIONMODEL_H
#define ABSORPTIONMODEL_H

#include <string>
#include <vector>

using namespace std;

class AbsorptionModel
{
    public:

        // compute coeffs

        virtual double ComputeCoeffs() = 0;
        
        // return absorption coeff

        virtual double a(int iBand, double T) = 0;

        // transmissivity 
        
        
        virtual double transmissivity(const vector <double>& T, const vector <double>& P, const vector <double>& XH2O, 
                                                 const vector <double>& XCO2, const vector <double>& delta, int iBand, double mu, int i, int j) = 0;
        
        
        virtual double kappa(int iBand, double T) = 0;
        virtual double kappa(int iBand, double pH2O, double pCO2) = 0;
        virtual double kappa(double T, double pH2O, double pCO2) = 0;

        virtual vector<double> getKappa() = 0;
        virtual vector<vector<double>> getKappa(int i) = 0;


        virtual vector<vector<double>> solve(vector<double> mu, int Npoints, vector <double> T, vector <double> P, 
                                                 vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign) = 0;

        // wavenumber
        virtual double nu(int iBand) = 0;

        // read data

        virtual void ReadData() = 0;

        // return number of bands

        virtual int NbBands() = 0;

        // Factory method to choose an absorption Model at run time

        static AbsorptionModel *make_absorptionModel( string absorptionName);

    private:

        // Absorption Model

        AbsorptionModel* absorption;
};


#endif
