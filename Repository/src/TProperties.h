#ifndef TPROPERTIES_H__
#define TPROPERTIES_H__

#include "Config.h"
#include "Constants.h"
#include "Input.h"
#include "TRadiation.h"
#include "FlameNode.h"
#include "RadiationModel.h"

class TProperties {	// contains properties of the gas mixture
// TFlame contains a function for the complete calculation of all species and mixture properties
public:
	template<class Species>
	TProperties( TInputDataPtr input, Species* species ) 
			{ InitProperties( input, species ); }
	~TProperties( void );

	TRadiationPtr 	GetRadiation( void ) { return fRadiation; };
	void 			ComputeMixtureMolarMass( Double &mixMolarMass, Double *Y, Double *molarMass, int nSpecies );
	template<class Species>
	void			ComputeDCpDT( Double &dCpdT, Double *Y, Double temp, Species* species );
/*	// 0D stuff*/
/*	void			CompMixtureProps( Double &mixHeatCapacity, Double &density, Double *heatCapacity, Double temp, Double mixMolarMass, Double pressure );*/
    string          fRadiationName; //Rui

protected:
	template<class Species>
	void 			InitProperties( TInputDataPtr input, Species* species );
/*	template<class Species>
	void			FillProperties( Species* species );*/
	
	TRadiationPtr	fRadiation;			
/*	Double			**sqrtMjOverMi;			//	contains sqrt( Mj/Mi ); last element is sqrtMjOverMi[nOfSpecies][nOfSpecies]*/
/*	Double			**sqrtMiOverMjPlusOne;	//	contains 0.3535534 / sqrt( Mi/Mj + 1 ); last element is sqrtMiOverMjPlusOne[nOfSpecies][nOfSpecies]*/
};
typedef TProperties *TPropertiesPtr;



class T0DProperties : public TProperties {	// contains properties of the gas mixture
public:
	template<class Species>
	T0DProperties( TInputDataPtr input, Species* species )
		: TProperties( input, species )
			{ InitT0DProperties( input ); }
	~T0DProperties( void );
	string          fRadiationName;
	Double			GetDensity( void ) { return fDensity; };
	Double			&GetDensityRef( void ) { return fDensity; };
	Double			GetMixHeatCapacity( void ) { return fMixHeatCapacity; };
	Double			GetMixConductivity( void ) { return fMixConductivity; };
	Double	 		GetMixViscosity( void ) { return fMixViscosity; };
	Double			GetMixMolarMass( void ) { return fMixMolarMass; };
	Double			&GetMixMolarMassRef( void ) { return fMixMolarMass; };
	Double			GetPressure( void ) { return fPressure; };
	Double			&GetPressureRef( void ) { return fPressure; };

	string          GetRadiationName(void){return fRadiationName;};  

	T0DRadiationPtr GetRadiation( void ) { return fRadiation; };

	template<class Species>
	void			CompMixtureProps( Double *heatCapacity, Double *conductivity, Double *mu, Double *Y, Double temp
			, Double &pressure, Double &density, EqOfState what, int nSpeciesIn, Species* species );

private:
	void 			InitT0DProperties( TInputDataPtr input );

	Double			fPressure;
	Double			fDensity;
	Double			fMixHeatCapacity;
	Double			fMixViscosity;		// [kg / (m*s)]
	Double			fMixConductivity;
	Double			fMixMolarMass;
	T0DRadiationPtr	fRadiation;
};
typedef T0DProperties *T0DPropertiesPtr;

class T1DProperties : public TProperties {	// contains properties of the gas mixture
public:
	template<class Species>
	T1DProperties( TInputDataPtr input, Species* species ) 
			: TProperties( input, species )
			{ InitT1DProperties( input ); }
	~T1DProperties( void );

	VectorPtr 		GetViscosity( void ) { return fViscosity; };
	VectorPtr 		GetDensity( void ) { return fDensity; };
	VectorPtr 		GetConductivity( void ) { return fConductivity; };
	VectorPtr 		GetHeatCapacity( void ) { return fHeatCapacity; };
	VectorPtr 		GetMolarMass( void ) { return fMolarMass; };
	VectorPtr 		GetDCpdT( void ) { return fdCpdT; };
	T1DRadiationPtr GetRadiation( void ) { return fRadiation; };
	VectorPtr	 	GetDCpDT( void ) { return fdCpdT; };
	template<class Species>
	void 			CompMixtureProps( TFlameNodePtr flameNode, Double *Y, Double temp
										, Double pressure, Species* species );
	void			PrintProperties( int k );
	void			PrintProperties( TNewtonPtr bt );
private:
	void 			InitT1DProperties( TInputDataPtr input );

	VectorPtr		fViscosity;		// [kg / (m*s)]		//	vector of length maxGridPoints+2, where left boundary is fViscosity[-1] 
	VectorPtr		fDensity;		// [kg / m^3]		//	vector of length maxGridPoints+2, ...
	VectorPtr		fConductivity;		// [W / (m*K)]	//	vector of length maxGridPoints+2, ...
	VectorPtr		fHeatCapacity;		// [J / (kg*K)]	//	vector of length maxGridPoints+2, ...
	VectorPtr		fMolarMass;		// [kg / kmol]		//	vector of length maxGridPoints+2, ...
	T1DRadiationPtr	fRadiation;		
	VectorPtr		fdCpdT;			// 				//	vector of length maxGridPoints, ...
};
typedef T1DProperties *T1DPropertiesPtr;

#include "TProperties.hpp"

#endif // TPROPERTIES_H__
