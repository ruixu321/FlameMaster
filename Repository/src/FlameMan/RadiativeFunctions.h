#ifndef RADIATIVEFUNCTIONS_H
#define RADIATIVEFUNCTIONS_H

#include <string>
#include <iostream>
#include <vector>

using namespace std;


#define C1 3.7418E-16
#define C2 1.4388

inline double BlackbodyIntensityWaveNumber(int nu,double T){

    return C1*(nu*100)*(nu*100)*(nu*100)/(exp(C2*nu/T) -1)/M_PI*100*25;

}

inline double BlackBodyIntensity(double T){

    const double sigma = 5.670E-8;
    return sigma*pow(T,4)/M_PI;
}

#endif
