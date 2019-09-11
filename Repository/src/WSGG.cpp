#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>      // std::istringstream
#include <cmath> 
#include <algorithm>

#include "WSGG.h"
#include "RadiativeFunctions.h"

using namespace std;

WSGG::WSGG(){

    b = vector<vector<double>> (25, vector<double>(12, 0));
    // read the coefficients for wsgg
    ReadData();

    kappaOutput = vector <double> (0, 0.); 
}


/**
 * Read the coefficients for wsgg superposition model of Cassol (2013)
 */

void WSGG::ReadData(){

    ifstream infile("./data/bands.dat"); //maybe add a file to specify the file where are the coeffs
    
    if (infile) {
        string line;
        int nline = 0;
        while (getline(infile, line)) {
            istringstream iss(line);

            // hard coding 12 is number of coeffs per band
            
            for(int i =0; i<12; i++){
                iss >> b[nline][i];
            }
            nline++;
        }
    }
}

double WSGG::a(){

}

int WSGG::NbBands(){
    return NBands;
}

double WSGG::a(int iBand, double T){

    return (b[iBand][2]*pow(T,0)*1e-0 + b[iBand][4]*pow(T,1) + b[iBand][6]*pow(T,2) + b[iBand][8]*pow(T,3) + 
        b[iBand][10]*pow(T,4))*(b[iBand][3]*pow(T,0)*1e-0 + b[iBand][5]*pow(T,1) + b[iBand][7]*pow(T,2) + 
        b[iBand][9]*pow(T,3) + b[iBand][11]*pow(T,4));
}

double WSGG::kappa(int iBand, double pH2O, double pCO2){
    
    if (pCO2 < 1E-6 && pH2O < 1E-6){
        return 0.;
    }
    else if(pCO2<1E-6){
        return (b[iBand][0]*pH2O);
    }
    else if (pH2O < 1E-6){
        return (b[iBand][1]*pCO2);
    }
    
    else{
        return (b[iBand][0]*pH2O+b[iBand][1]*pCO2);
    }
}

vector<vector<double>> WSGG::solve(vector<double> mu, int NPoints, vector <double> T, vector <double> P, 
    vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign){

    vector<vector<vector<double>>> I(mu.size(), vector<vector<double>> (NbBands(), vector<double> (NPoints,0)));
    vector<vector<double>> ITot(mu.size(), vector<double> (NPoints,0));

    if(sign=="+"){

        kappaOutput = vector <double> (XCO2.size(), 0.); 

    } 

    for(int muj=0; muj<mu.size(); muj++){

        for(int i=0; i<NPoints; i++){

            for(int iBand =0; iBand<NbBands(); iBand++){

                if(i==0){

                    I[muj][iBand][i] = a(iBand,T[i])*BlackBodyIntensity(T[i]);
                }
                else{
                    double pH2O = XH2O[i-1]*P[i-1];
                    double pCO2 = XCO2[i-1]*P[i-1];

                    I[muj][iBand][i] = a(iBand,T[i])*BlackBodyIntensity(T[i]) + 
                    (I[muj][iBand][i-1] - a(iBand,T[i])*BlackBodyIntensity(T[i]))*
                    exp(-kappa(iBand,pH2O,pCO2)*DeltaX[i-1]/mu[muj]); 

                    if(sign =="+"){
                        kappaOutput[i-1] += a(iBand,T[i]) * kappa(iBand,pH2O,pCO2);
                    }

                }

                ITot[muj][i] = (iBand==0) ? I[muj][iBand][i] : ITot[muj][i] + I[muj][iBand][i]; 

            }
        }

        if(sign == "-"){

            reverse(ITot[muj].begin(), ITot[muj].end());
        }
    }

    return ITot; 

}

vector<double> WSGG::getKappa(){

    return kappaOutput;

}

vector < vector<double> > WSGG::getKappa(int){

    // return kappaOutput;

}
