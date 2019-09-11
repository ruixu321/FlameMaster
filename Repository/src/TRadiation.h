#ifndef TRADIATION_H__
#define TRADIATION_H__

#include "Config.h"
#include "Constants.h"
#include "Input.h"
#include "RadiationModel.h"

class TRadiation {
public:
	TRadiation( TInputDataPtr input ) : fExpFuncCoeffCO2(-8.888e-4),fExpFuncCoeffH2O(-1.546e-3)
					,fAlphaCoeffCO2(46.241e-5),fAlphaCoeffH2O(22.6e-5){ InitRadiation( input ); };
	~TRadiation( void );

	void 		ComputeRadiationOnePoint ( Double *radiation, Double temp, Double *Y, Double *molarMass, Double	density );
	void 		ComputeRadiationOnePoint ( Double *radiation, Double temp, Double *Y, Double *molarMass, Double	density, 
                                                                           Double HRR, Double RadiativeFrac , Double *kappa);
	void 		ComputeRadiationOnePoint ( Double *radiation, double dqr);
	void 		ComputeRadiationRadiation( Double *temp, vector<double> XH2O, vector<double> XCO2, int nSpecies, 
                                                         vector<double> physicalGrid, double Delta, double chiRef, vector<double> HRR, string RadiationName);
    void        ComputeRadiationUnsteady ( vector<double>temp, vector<double> XH2O, vector<double> XCO2, int nSpecies, 
                                                         vector<double> physicalGrid, double Delta, double chiRef, vector<double> HRR, string RadiationName);

    vector <double> dqrF; 
    vector <double> kappa; 

    vector<vector<double>> kappa1; 
    vector<vector<double>> kappa2; 
    vector<vector<double>> kappa3; 
    vector<vector<double>> kappa4; 

protected:
	void 	InitRadiation( TInputDataPtr input );

	int			fH2OIndex;
	int			fCO2Index;
	int			fCH4Index;
	int			fCOIndex;
	const Double	fExpFuncCoeffCO2; // = -8.888e-4
	const Double	fExpFuncCoeffH2O; // = -1.546e-3
	const Double	fAlphaCoeffCO2; // = 46.241e-5
	const Double	fAlphaCoeffH2O; // = 22.6e-5

    // radiation model
    //
    //RadiationModel* radModel;
};
typedef TRadiation *TRadiationPtr;


class T0DRadiation : public TRadiation {
    public:
        T0DRadiation( TInputDataPtr input ) : TRadiation( input ) { InitT0DRadiation( input ); };
        ~T0DRadiation( void ) {};

        Double		GetRadiation( void ) { return fRadiation; };
        Double		GetKappa( int i ) { return fKappa; };
        void		SetRadiation( Double temp, Double *Y, Double *molarMass, Double density ) 
        { ComputeRadiationOnePoint( &fRadiation, temp, Y, molarMass, density ); };

        
    private:
        void 	InitT0DRadiation( TInputDataPtr /*input*/ ) {};
        Double	fRadiation;
        Double	fKappa;
};
typedef T0DRadiation *T0DRadiationPtr;



class T1DRadiation : public TRadiation {
    public:
        T1DRadiation( TInputDataPtr input ) : TRadiation( input ) { InitT1DRadiation( input );};
        ~T1DRadiation( void );

        VectorPtr	GetRadiation( void ) { return fRadiation; };
        VectorPtr	GetKappa() { return fKappa; };
        // VectorPtr   GetKappa( int i ) { return fKappa; };
        template<class Flame>
        void		FillJacRadiation( Double coeff, Flame* flame, NodeInfoPtr nodeInfo );
    private:
        void 	    InitT1DRadiation( TInputDataPtr input );
        VectorPtr	fRadiation;		//	length is nOfGridPoints
        VectorPtr	fKappa;		    //	length is nOfGridPoints
};
typedef T1DRadiation *T1DRadiationPtr;


    template<class Flame>
void T1DRadiation::FillJacRadiation( Double coeff, Flame* flame, NodeInfoPtr nodeInfo )
{
    static int init = 0;

    if ( !init ) {
        fprintf( stderr, "###warning: function T1DRadiation::FillJacRadiation not available for RADCAL radiation model\n" );
        fprintf( stderr, "###warning: might influence convergence, not results\n" );
        init = 1;
    }
}

#endif // TRADIATION_H__
