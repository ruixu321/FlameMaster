#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>      // std::istringstream
#include <cmath> 
#include <algorithm>

#include "WSGGBordbar.h"
#include "RadiativeFunctions.h"

#define ncoeffs 24

using namespace std;

WSGGBordbar::WSGGBordbar(){

    //c = vector<vector<vector<double>>> (4, vector<vector<double>>(5, vector<double>(5, 0.)));
    //d = vector<vector<double>> (4, vector<double>(5,0.));
    // read the coefficients for wsgg
    //ReadData();


    c={7.4129560000e-001,-5.2444410000e-001,5.8228600000e-001,-2.0969940000e-001,2.4203120000e-002, 
	    -9.4126520000e-001,2.7995770000e-001,-7.6723190000e-001,3.2040270000e-001,-3.9101740000e-002,
	    8.5318660000e-001,8.2307540000e-002,5.2894300000e-001,-2.4684630000e-001,3.1093960000e-002,
	    -3.3428060000e-001,1.4749870000e-001,-4.1606890000e-001,1.6976270000e-001,-2.0406600000e-002,
	    4.3143620000e-002,-6.8862170000e-002,1.1097730000e-001,-4.2086080000e-002,4.9188170000e-003,
	    1.5520730000e-001,-4.8621170000e-001,3.6680880000e-001,-1.0555080000e-001,1.0585680000e-002,
	    6.7556480000e-001,1.4092710000e+000,-1.3834490000e+000,4.5752100000e-001,-5.0197600000e-002,
	    -1.1253940000e+000,-5.9131990000e-001,9.0854410000e-001,-3.3342010000e-001,3.8423610000e-002,
	    6.0405430000e-001,-5.5338540000e-002,-1.7330140000e-001,7.9160830000e-002,-9.8933570000e-003,
	    -1.1054530000e-001,4.6466340000e-002,-1.6129820000e-003,-3.5398350000e-003,6.1212770000e-004,
	    2.5502420000e-001,3.8054030000e-001,-4.2497090000e-001,1.4294460000e-001,-1.5740750000e-002,
	    -6.0654280000e-001,3.4940240000e-001,1.8535090000e-001,-1.0136940000e-001,1.3024410000e-002,
	    8.1238550000e-001,-1.1020090000e+000,4.0461780000e-001,-8.1182230000e-002,6.2981010000e-003,
	    -4.5322900000e-001,6.7844750000e-001,-3.4326030000e-001,8.8308830000e-002,-8.4152210000e-003,
	    8.6930930000e-002,-1.3069960000e-001,7.4144640000e-002,-2.0292940000e-002,2.0109690000e-003,
	    -3.4519940000e-002,2.6567260000e-001,-1.2253650000e-001,3.0015080000e-002,-2.8205250000e-003,
	    4.1120460000e-001,-5.7283500000e-001,2.9244900000e-001,-7.9807660000e-002,7.9966030000e-003, 
	    -5.0559950000e-001,4.5795590000e-001,-2.6164360000e-001,7.6484130000e-002,-7.9083560000e-003,
	    2.3175090000e-001,-1.6567590000e-001,1.0526080000e-001,-3.2193470000e-002,3.3869650000e-003, 
	    -3.7549080000e-002,2.2951930000e-002,-1.6004720000e-002,5.0463180000e-003,-5.3643260000e-004}; 

    d={3.4042880000e-002,6.5230480000e-002,-4.6368520000e-002,1.3868350000e-002,-1.4449930000e-003,
	    3.5094570000e-001,7.4651380000e-001,-5.2930900000e-001,1.5944230000e-001,-1.6632610000e-002,
	    4.5707400000e+000,2.1680670000e+000,-1.4989010000e+000,4.9171650000e-001,-5.4299900000e-002,
	    1.0981690000e+002,-5.0923590000e+001,2.3432360000e+001,-5.1638920000e+000,4.3938890000e-001};

    kappaOutput = vector <double> (0, 0.); 

    kappaOutput1 = vector<vector<double>> ( 2 , vector <double> (0, 0.) ); 
    kappaOutput2 = vector<vector<double>> ( 2 , vector <double> (0, 0.) );
    kappaOutput3 = vector<vector<double>> ( 2 , vector <double> (0, 0.) );
    kappaOutput4 = vector<vector<double>> ( 2 , vector <double> (0, 0.) );

}


/**
 * Read the coefficients for wsgg superposition model of Cassol (2013)
 */

void WSGGBordbar::ReadData(){

	//cout << "Start reading coefficients" << endl; 
	//ifstream infile("./data/Bordbar.dat"); //maybe add a file to specify the file where are the coeffs

	//if (infile) {
	//	string line;

	//	int nline = 0;

	//	for (int i =0; i<4; ++i){
	//		for (int j=0; j<6; ++j){

	//			getline(infile, line);
	//			istringstream iss(line);

	//			if(j!=5){
	//				//cout << "line " << nline << endl;
	//				for(int k=0; k<5; ++k){
	//					iss >> c[i][j][k];
	//				}
	//			}
	//			else{	
	//				for(int k=0; k<5; ++k){
	//					iss >> d[i][k];
	//				}
	//			}
	//			nline++;
	//		}
	//	}
	//}

	//cout << "End reading coefficients" << endl; 
}

double WSGGBordbar::a(){

	return 0.;

}

int WSGGBordbar::NbBands(){
	return NBands;
}

double WSGGBordbar::a(int iBand, double T){

	return 0.;
}

double WSGGBordbar::a(int i, double T, double MR){
	return b(i,0,MR) + b(i,1,MR)*pow(T,1) + b(i,2,MR)*pow(T,2) + b(i,3,MR)*pow(T,3) + b(i,4,MR)*pow(T,4);
}

double WSGGBordbar::b(int i,  int j, double MR){

	//cout << i << " " << j << " " << c[i][j][0] << " " <<  c[i][j][1] << " " <<  c[i][j][2] << " " <<  c[i][j][3] << " " <<  c[i][j][4] << endl;
	return c[i*25+j*5+0] + c[i*25+j*5+1]*MR + c[i*25+j*5+2]*pow(MR,2) + c[i*25+j*5+3]*pow(MR,3) + c[i*25+j*5+4]*pow(MR,4);

}


double WSGGBordbar::kappa(int i, double MolarFracH2O, double MolarFracCO2){

	double MR = MolarFracH2O / MolarFracCO2;

	if (MolarFracCO2 < 1E-6 && MolarFracH2O < 1E-6){
        return 0.0;
    }
    else{
		return d[i*5+0] + d[i*5+1]*MR + d[i*5+2]*pow(MR,2) + d[i*5+3]*pow(MR,3) + d[i*5+4]*pow(MR,4);
    }

}

vector<vector<double>> WSGGBordbar::solve(vector<double> mu, int NPoints, vector <double> T, vector <double> P, 
	vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign){

	vector<vector<vector<double>>> I(mu.size(), vector<vector<double>> (NbBands(), vector<double> (NPoints,0)));
	vector<vector<double>> ITot(mu.size(), vector<double> (NPoints,0)); 
	double MR;

    if(sign=="+"){

        kappaOutput = vector <double> (XCO2.size(), 0.);

        kappaOutput1 = vector<vector<double>> (2, vector <double> (XCO2.size(), 0.) );
    	kappaOutput2 = vector<vector<double>> (2, vector <double> (XCO2.size(), 0.) );
		kappaOutput3 = vector<vector<double>> (2, vector <double> (XCO2.size(), 0.) );
		kappaOutput4 = vector<vector<double>> (2, vector <double> (XCO2.size(), 0.) );	


    } 


	for(int muj=0; muj<mu.size(); muj++){

		for(int i=0; i<NPoints; i++){

			// compute a coefficients
			vector <double> acoeffs (5,0.); 
			vector <double> kappa_coeffs (5,0.); 
			double sum_of_coeffs = 0;

			for(int iBand =1; iBand<NbBands(); iBand++){
				if (i == 0){
					MR = 2.0;
				}
				else{
					MR = XH2O[i]/XCO2[i];
					if (MR > 4.0){ MR=4.0; }
					if (MR < 0.01){ MR=0.01; }
					kappa_coeffs[iBand] = kappa(iBand-1, XH2O[i-1], XCO2[i-1]);
				}
				acoeffs[iBand] = a(iBand-1,T[i]/1200., MR);	
				sum_of_coeffs += acoeffs[iBand];
			}

			acoeffs[0] = 1 - sum_of_coeffs;
			kappa_coeffs[0] = 0; //clear gas

			for(int iBand =0; iBand<NbBands(); iBand++){

				if(i==0){

					I[muj][iBand][i] = acoeffs[iBand]*BlackBodyIntensity(T[i]);
				}
				else{
					I[muj][iBand][i] = acoeffs[iBand]*BlackBodyIntensity(T[i]) + 
					(I[muj][iBand][i-1] -acoeffs[iBand]*BlackBodyIntensity(T[i]))*
					exp(-kappa_coeffs[iBand] * (XH2O[i-1]+XCO2[i-1]) * (P[i-1]) * DeltaX[i-1]/mu[muj]); 

	                if(sign =="+" && iBand != 0){
	                	
	                    kappaOutput[i-1] += kappa_coeffs[iBand] * (XH2O[i-1]+XCO2[i-1]) * (P[i-1]) ;

	                    switch (iBand) {
							case 1: 
								kappaOutput1[0][i-1] = kappa_coeffs[iBand] * (XH2O[i-1]+XCO2[i-1]) * (P[i-1]) ;
								kappaOutput1[1][i-1] = acoeffs[iBand] * kappaOutput1[0][i-1] ;
								break;
							case 2: 
								kappaOutput2[0][i-1] = kappa_coeffs[iBand] * (XH2O[i-1]+XCO2[i-1]) * (P[i-1]) ;
								kappaOutput2[1][i-1] = acoeffs[iBand] * kappaOutput2[0][i-1] ;
								break;
							case 3: 
								kappaOutput3[0][i-1] = kappa_coeffs[iBand] * (XH2O[i-1]+XCO2[i-1]) * (P[i-1]) ;
								kappaOutput3[1][i-1] = acoeffs[iBand] * kappaOutput3[0][i-1] ;
								break;
							case 4: 
								kappaOutput4[0][i-1] = kappa_coeffs[iBand] * (XH2O[i-1]+XCO2[i-1]) * (P[i-1]) ;
								kappaOutput4[1][i-1] = acoeffs[iBand] * kappaOutput4[0][i-1] ;
								break;
							default:
								cout << "Error in calculating kappaOutput*:   " << "iBand = " << iBand << endl;	
								break;
						}
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


vector<double> WSGGBordbar::getKappa(){

    return kappaOutput;

}

vector < vector<double> > WSGGBordbar::getKappa(int i){

	switch (i) {
		case 1: 
			return kappaOutput1;
		case 2: 
			return kappaOutput2;
		case 3: 
			return kappaOutput3;
		case 4: 
			return kappaOutput4;	
		default:
			cout << "Error in WSGGBordbar::getKappa* " << endl;
			break;
	}
}
