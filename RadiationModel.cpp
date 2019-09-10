#include <cmath>
#include <fstream>
#include <algorithm>

#include "RadiationModel.h"
#include "RadiativeFunctions.h"

using namespace std;

RadiationModel::RadiationModel(string name){
    
}
RadiationModel::RadiationModel(string name, vector<double> gridInput, vector<double> Tinput, vector<double> XH2Oinput, 
        vector<double> XCO2input, vector<double> Pinput, string dqrFileNameInput, string qrFileNameInput, vector<double> HRRinput){

    //cout << "--> New Radiation Model " << endl;

    // grid 

    gridPoint = gridInput;
    NPoints = gridPoint.size();
    gridMidPoint = vector<double> (NPoints-1);
    DeltaX = vector<double> (NPoints-1);

    ComputeDeltaX();
    ComputeMidGrid();

    // absorption model
    absorptionName = name;
    absorption = AbsorptionModel::make_absorptionModel(absorptionName);

    if(radiationModelName == "DOM"){

        if(name == "Grey" || name == "WSGG" || name == "WSGGJohansson" || name == "WSGGBordbar"){
            // cout << "S8" << endl;
            mu = {0.1422555, 0.5773503,0.8040087,0.9795543};
            w = {2.1637144, 2.6406988,0.7938272,0.6849436};
        }
        else{
            // cout << "S1" << endl;    // This is for SNB model
            mu = {0.5000000};
            w = {6.2831853};
        }
    }

    // intensities

    Iminus = vector<vector<vector<double>>> (mu.size(), vector<vector<double>> (absorption->NbBands(), vector<double> (NPoints,0)));
    Iplus = vector<vector<vector<double>>> (mu.size(), vector<vector<double>> (absorption->NbBands(), vector<double> (NPoints,0)));

    IminusTot = vector<vector<double>> (mu.size(), vector<double> (NPoints,0)); 
    IplusTot = vector<vector<double>> (mu.size(), vector<double> (NPoints,0)); 

    // radiative heat flux

    qr = vector<double> (NPoints,0.);

    // radiative source term
    //
    dqr = vector<double> (NPoints-1, 0.);

    // input parameters

    XH2O = XH2Oinput; 
    XCO2 = XCO2input; 

    P = Pinput;
    T = Tinput; // need to add boundary condition for temperature

    // output

    dqrFileName = dqrFileNameInput;
    qrFileName = qrFileNameInput;

    HRR = HRRinput;

}

/**
 * Calculate the width of each cell based on a grid
 */

void RadiationModel::ComputeDeltaX(){

    for (int i = 0; i<DeltaX.size(); i++){

        DeltaX[i] = (gridPoint[i+1] - gridPoint[i]);

    }
}

void RadiationModel::ComputeMidGrid(){

    for (int i = 0; i<gridMidPoint.size(); i++){

        gridMidPoint[i] = (i==0) ? (gridPoint[i+1] - gridPoint[i])/2. : (gridPoint[i+1] - gridPoint[i])/2. + gridPoint[i];

    }
}

double RadiationModel::Solve(){

    // cout << "--> Solving RTE " << endl;

    //DOM 

    // cout << "----> Number of bands to solve " << absorption->NbBands() << endl;

    // compute the wsgg coefficients based on temperature


    absorption->ComputeCoeffs(); 

    // cout << " ----> I+ " << endl;

    SolveIPlus();

    // cout << " ----> I- " << endl;

    SolveIMinus();

    // cout << " Compute radiative heat flux " << endl;

    ComputeRadiativeHeatFlux();

    // cout << " Compute radiation source term "  << endl;

    ComputeRadiativeSource();

    // cout << " Writting data  " <<  endl;

    // WriteData();

    // cout << " End " << endl;

    return 0;
}

void RadiationModel::SolveIPlus(){

    IplusTot = absorption->solve(mu, NPoints, T, P, XH2O, XCO2, DeltaX, "+");

}

void RadiationModel::SolveIMinus(){

    reverse(P.begin(), P.end());
    reverse(T.begin(), T.end());
    reverse(XH2O.begin(), XH2O.end());
    reverse(XCO2.begin(), XCO2.end());
    reverse(DeltaX.begin(), DeltaX.end());


    IminusTot = absorption->solve(mu, NPoints, T, P, XH2O, XCO2, DeltaX, "-");


    reverse(P.begin(), P.end());
    reverse(T.begin(), T.end());
    reverse(XH2O.begin(), XH2O.end());
    reverse(XCO2.begin(), XCO2.end());
    reverse(DeltaX.begin(), DeltaX.end());

}

/*
 * Calculate the radiative heat flux based on Iminus and Iplus
 */

void RadiationModel::ComputeRadiativeHeatFlux(){

    for(int i = 0; i<NPoints; i++){

        qr[i] = 0.;

        for(int muj = 0; muj<mu.size(); muj++){

            qr[i] += mu[muj]*w[muj]*(IplusTot[muj][i] - IminusTot[muj][i]); 
        }
    }
}

/*
 * Calculate the Radiative source term based on the radiative heat flux qr
 */

void RadiationModel::ComputeRadiativeSource(){

    for (int i = 0; i<NPoints-1; i++){

        dqr[i] = (qr[i+1] - qr[i])/DeltaX[i];
    }

}

/* * Write data to files
 * qr
 * div qr
 */

// void RadiationModel::WriteData(){

//     ofstream qrfile (qrFileName);
//     ofstream dqrfile (dqrFileName);

//     vector<double> kappaOutput = GetKappa();

//     if (qrfile.is_open()){

//         qrfile << "x [m]\tqr [W/m2]" << endl;
//         for (int idx = 0; idx<qr.size(); idx++){

//             qrfile << gridPoint[idx] << "\t" << qr[idx] << endl;
//         }

//         qrfile.close();

//     }
//     if(absorptionName == "Grey" || absorptionName == "WSGG" || absorptionName == "WSGGJohansson" || absorptionName == "WSGGBordbar" ){
//         if (dqrfile.is_open()){
//             dqrfile << "x [m]\tdqr [W/m3]\tX-CO2\tX-H2O\ttemperature [K]\t kappa [1/m]" << endl;

//             for (int idx = 0; idx<dqr.size(); idx++){
//                 dqrfile << gridMidPoint[idx] << "\t" << dqr[idx] <<  "\t" << XCO2[idx] << "\t" << XH2O[idx] << "\t" << T[idx+1] << "\t" << kappaOutput[idx] << endl;
//             }

//             dqrfile.close();
//         }
//     }
//     else{
//     if (dqrfile.is_open()){
//             dqrfile << "x [m]\tdqr [W/m3]\tX-CO2\tX-H2O\ttemperature [K]" << endl;

//             for (int idx = 0; idx<dqr.size(); idx++){
//                 dqrfile << gridMidPoint[idx] << "\t" << dqr[idx] <<  "\t" << XCO2[idx] << "\t" << XH2O[idx] << "\t" << T[idx+1] << endl;
//             }

//             dqrfile.close();
//         }

//     }

// }


vector<double> RadiationModel::GetKappa(){
    return absorption->getKappa(); 
}

vector<vector<double>> RadiationModel::GetKappa(int i){
    return absorption->getKappa(i); 
}

