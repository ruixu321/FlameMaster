#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>      // std::istringstream
#include <cmath> 
#include <tuple>
#include <algorithm>

#include "SNB.h"
#include "RadiativeFunctions.h"

using namespace std;

#define H2OSize1 449
#define H2OSize2 44002
#define Ncols 49
#define CO2Size1 323
#define CO2Size2 31654



SNB::SNB(){

    // b = vector<vector<double>> (NBands, vector<double>(12, 0));
    
    //dataCO2 = vector<vector<double>> (323*2, vector<double>(49, 1));

    ReadData();
}



inline void SNB::dataInterpolation(double T, int & IT, double & RT){

    //int IT;
    //double RT;
    if(T>300){
        if(T<5000){
            RT = (T-300.)/100.;
            IT = int(RT+1.0E-6);
            RT = RT-IT;
            IT = IT +1;
        }
        else{
            RT = 1.;
            IT = 46;
        }
    }
    else{

        RT = 0;
        IT = 0;
    }

    //return make_tuple(RT, IT);
}

void SNB::readH2O(){
    std::ifstream infileH2O("data/SNBH2O");
    if (infileH2O) {
        string line;
        int nline = 0;
        while (getline(infileH2O, line)) {
            istringstream iss(line);
            for(int i =0; i<Ncols; i++){
                iss >> dataH2O[nline*Ncols+i];
            }
            nline++;
        }
    }
}

void SNB::readCO2(){
    std::ifstream infileH2O("data/SNBCO2");
    if (infileH2O) {
        string line;
        int nline = 0;
        while (getline(infileH2O, line)) {
            istringstream iss(line);
            for(int i =0; i<Ncols; i++){
                iss >> dataCO2[nline*Ncols+i];
            }
            nline++;
        }
    }
}

void SNB::ReadData(){
    readCO2();
    readH2O();

}
double SNB::a(){

}

inline int SNB::NbBands(){
    return NBands;
}

double SNB::a(int iBand, double T){

    return 0;
}

double SNB::kappa(int iBand, double pH2O, double pCO2){

    return 0;
}

double SNB::transmissivity(const vector <double>& T, const vector <double>& P, const vector <double>& XH2O, const vector <double>& XCO2, const vector <double>& delta, int iBand, double mu, int idxStart, int idxEnd){

    double tau = 1.;

    int iBandCO2 = iBand - 9;

    double XN2 = 0;

    int IT;
    double RT;

    double GAMH2O = 0; 
    double YCH2O = 0;
    double XBH2O = 0;
    double YBH2O = 0;
    double YKH2O = 0;
    double ZCH2O = 0;
    double ZBH2O = 0;
    double ZKH2O = 0;
    double XKH2O = 0;
    double XDH2O = 0;

    double GAMCO2 = 0; 
    double YCCO2 = 0;
    double YKCO2 = 0;
    double XBCO2 = 0;
    double YBCO2 = 0;
    double ZCCO2 = 0;
    double ZBCO2 = 0;
    double ZKCO2 = 0;
    double XKCO2 = 0;
    double XDCO2 = 0;

    double SCCO2 = 0.;
    double SBCO2 = 0.;
    double SKCO2 = 0.;

    double SCH2O = 0.;
    double SKH2O = 0.;
    double SBH2O = 0.;

    int indexH2O1 = 0;
    int indexH2O2 = 0;
    int indexCO21 = 0;
    int indexCO22 = 0;

    double Tnum;
    double Tdenum;
    double T296 = Tnum/Tdenum;


    if(iBand < 449){

        for (int idxPoint = idxStart; idxPoint < idxEnd; idxPoint++){ 

            XN2 = 1. - XCO2[idxPoint] - XH2O[idxPoint];

            // find corresponding data in the database based on the index idx

            dataInterpolation(T[idxPoint+1], IT, RT); // don t forget about BC 

            //H2O
            //

            indexH2O1 = iBand*Ncols + IT;
            indexH2O2 = (iBand+H2OSize1)*Ncols + IT;
            indexCO21 = iBandCO2*Ncols + IT;
            indexCO22 = (iBandCO2+CO2Size1)*Ncols + IT;

            XKH2O=dataH2O[indexH2O1]+RT*(dataH2O[indexH2O1+1]-dataH2O[indexH2O1]);
            XDH2O=dataH2O[indexH2O2]+RT*(dataH2O[indexH2O2+1]-dataH2O[indexH2O2]);

            GAMH2O = P[idxPoint]*0.066*(7.*sqrt(296./T[idxPoint+1])*XH2O[idxPoint]+1.2*(XH2O[idxPoint]+XN2)+1.6*XCO2[idxPoint])*sqrt(296./T[idxPoint+1]);

            YCH2O =XH2O[idxPoint]*P[idxPoint]*delta[idxPoint]*100/mu;
            XBH2O=2*GAMH2O*XDH2O;
            YKH2O=YCH2O*XKH2O;
            YBH2O=YKH2O*XBH2O;
            SCH2O=SCH2O+YCH2O;

            if (SCH2O > 1.E-12){
                SKH2O=SKH2O+YKH2O;
                SBH2O=SBH2O+YBH2O;
                ZCH2O=SCH2O;
                ZKH2O=SKH2O/ZCH2O;
                ZBH2O=SBH2O/SKH2O;
                XKH2O=1.+2.*ZKH2O*ZCH2O/ZBH2O;
                tau = exp(-ZBH2O*(sqrt(XKH2O)-1.));
            }


            //CO2
            if ((iBand  >= 9) && (iBand < 323 + 9)){


                XKCO2=dataCO2[indexCO21]+RT*(dataCO2[indexCO21+1]-dataCO2[indexCO21]);
                XDCO2=dataCO2[indexCO22]+RT*(dataCO2[indexCO22+1]-dataCO2[indexCO22]);

                // calculate the transmissivity
                 
                Tnum = 296.*296.*296.*296.*296.*296.*296.*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1];
                Tdenum = 296.*296.*296.*296.*296.*296.*296.*296.*296.*296.*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1]*T[idxPoint+1];

                GAMCO2=P[idxPoint]*(0.07*XCO2[idxPoint]+0.058*XN2+0.10*XH2O[idxPoint])*pow(296./T[idxPoint+1],0.7);
                YCCO2 =XCO2[idxPoint]*P[idxPoint]*delta[idxPoint]*100/mu;
                XBCO2=2*GAMCO2*XDCO2;
                YKCO2=YCCO2*XKCO2;
                YBCO2=YKCO2*XBCO2;
                SCCO2=SCCO2+YCCO2;

                if (SCCO2 > 1.E-12){
                    SKCO2=SKCO2+YKCO2;
                    SBCO2=SBCO2+YBCO2;
                    ZCCO2=SCCO2;
                    ZKCO2=SKCO2/ZCCO2;
                    ZBCO2=SBCO2/SKCO2;
                    XKCO2=1+2*ZKCO2*ZCCO2/ZBCO2;
                    tau=tau*exp(-ZBCO2*(sqrt(XKCO2)-1.));
                }

            } // end CO2

        } // end H2O
    }

    return tau;

}

inline double SNB::nu(int iBand){

    return (50 + iBand*25); 
}



vector<vector<double>> SNB::solve(vector<double> mu, int NPoints, vector <double> T, vector <double> P, vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign){


    //vector <double> taui(NPoints, 1.);
    //vector <double> tauip1(Npoints, 1.);

    double taui = 1.;
    double tauip1 = 1.;

    //vector<vector<vector<double>>> I(mu.size(), vector<vector<double>> (NbBands(), vector<double> (NPoints,0)));
    //
    double I = 0.;
    vector<vector<double>> ITot(mu.size(), vector<double> (NPoints,0)); 


    for(int muj=0; muj<mu.size(); ++muj){

        for(int iBand =0; iBand<NbBands(); ++iBand){

		cout << muj << " " << nu(iBand) << endl;

		for(int i=0; i<NPoints; i++){

			if (i == 0){
				I = BlackbodyIntensityWaveNumber(nu(iBand), T[0]);
			}

			else{
				I = BlackbodyIntensityWaveNumber(nu(iBand), T[0])*transmissivity(T, P, XH2O, XCO2, DeltaX, iBand, mu[muj], 0, i);

				taui =  transmissivity(T, P, XH2O, XCO2, DeltaX, iBand, mu[muj], 0, i);

				for(int j=0; j<i;++j){

					if(j+1 == i){

						I +=  BlackbodyIntensityWaveNumber(nu(iBand), T[j+1])*(1-taui); 

					}
					else{

						tauip1 = transmissivity(T, P, XH2O, XCO2, DeltaX, iBand, mu[muj], j+1, i);
						I +=  BlackbodyIntensityWaveNumber(nu(iBand), T[j+1])*(tauip1-taui); 
					}

					taui = tauip1;

				}

				ITot[muj][i] = (iBand==0) ? I : ITot[muj][i] + I; 

			}

		}
	}

	if(sign == "-"){

		reverse(ITot[muj].begin(), ITot[muj].end());

	}
    }


    return ITot; 

}
