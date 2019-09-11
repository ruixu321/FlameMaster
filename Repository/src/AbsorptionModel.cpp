#include <string>
#include <iostream>

#include "AbsorptionModel.h"
#include "WSGG.h"
#include "SNB.h"
#include "Grey.h"
#include "WSGGJohansson.h"
#include "WSGGBordbar.h"


using namespace std;

AbsorptionModel *AbsorptionModel::make_absorptionModel(string name){

    if(name == "WSGG"){
        // cout << " Creation of the WSGG model " << endl;
        return new WSGG;
    }
    else if(name == "SNB"){
        // cout << " Creation of the SNB model " << endl;
        return new SNB;
    }
    else if(name == "Grey"){
        // cout << " Creation of the Grey model " << endl;
        return new Grey;
    }
    else if(name == "WSGGJohansson"){
        // cout << " Creation of the WSGGJohansson model " << endl;
        return new WSGGJohansson;
    }
    else if(name == "WSGGBordbar"){
        // cout << " Creation of the WSGGBordbar model " << endl;
        return new WSGGBordbar;
    }
    else{
        cout << "error: no " << name << " Model" << endl;
    }
}

