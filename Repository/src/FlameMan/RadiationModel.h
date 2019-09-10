#ifndef RADIATIONMODEL_H
#define RADIATIONMODEL_H

#include <string>
#include <iostream>
#include <vector>

#include "AbsorptionModel.h"

using namespace std;

class RadiationModel
{

    private:

        string radiationModelName = "DOM";

        // Absorption model name
        string absorptionName;

        // Grid

        vector<double> gridPoint;
        vector<double> gridMidPoint;

        vector<double> DeltaX;

       
        // Parameters of DOM method
        vector<double> mu = {0.1422555, 0.5773503,0.8040087,0.9795543};
        vector<double> w = {2.1637144, 2.6406988,0.7938272,0.6849436};

        // Intensities

        vector<vector<vector<double>>> Iminus;
        vector<vector<vector<double>>> Iplus;

        vector<vector<double>> IminusTot;
        vector<vector<double>> IplusTot;

        // radiative heat flux

        vector <double> qr;

        // radiative source term

        vector<double> dqr;

        // Mole fraction of species

        vector<double> XH2O; 
        vector<double> XCO2; 

        // Temperature

        vector<double> T;

        // Pressure

        vector<double> P;

        // Heat Release Rate

        vector<double> HRR;

        // output
        string qrFileName;
        string dqrFileName;

    protected:
        
        // Nbands
        int NBands;

        // NPoints
        int NPoints = 0;

        // Absorption Model
        AbsorptionModel* absorption;

    public:


        // constructor
        RadiationModel(string absorptionName, vector<double> gridMiddlePoint, vector<double> T, vector<double> XCO2, vector<double> XH2O, vector<double> P,
                string dqrFileName, string qrFileName, vector<double> HRR);
        ~RadiationModel( void ){ delete absorption;}
        // empty constructor
        RadiationModel(){};
        RadiationModel(double){};
        RadiationModel(vector<double>){};
        RadiationModel(string);

        // solve the radiative transfert equation
        double Solve();

        // compute grid for DOM
        void ComputeDeltaX();

        // Calculat the middle of the cells
        void ComputeMidGrid();

        // compute intensities

        void SolveIPlus();
        void SolveIMinus();

        // radiative heat flux
        void ComputeRadiativeHeatFlux();

        // radiative source term
        void ComputeRadiativeSource();

        //output
        void WriteData();

        // access function

        vector<double> GetRadSource(){return dqr;}
        vector<double> GetHeatFlux(){return qr;}
        vector<vector<double>> GetKappa(int i);
        vector<double> GetKappa();
};


#endif
