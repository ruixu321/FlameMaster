#ifndef NODEINFO_H__
#define NODEINFO_H__

#include "Config.h"  // typedef Double and Flag

class TFlameNode {
public:
	TFlameNode( void ) { InitTFlameNode(); };

// reaction
	Double		*tBodyConc;
	Double		*rateCoeff;
	Double		*tempReaction;
	Double		*pressureReaction;
	Double		*currRateCoeff;
	Flag		*kNewerThanW;
	Double		*reactionRate;
	Double		*YReaction;
	Double		*currReacRate;
	Double		**dMdY;
	
// species
	Double			*viscosity;
	Double			*heatCapacity;
	Double			*conductivity;
	Double			*enthalpy;
	Double			*diffusivityPrev;
	Double			*diffusivity;
	Double			*diffusivityNext;
	Double			**diffTherm;
	Double			*productionRate;
	Double			*savedY;
	Double			*savedDeltaiY;
	Double			*sumDiff;
	Double			*deltaI;
	Double			*tempProp;
	Double			*pressureProp;
	Double			***binDij;
	Double			**GijOverWj;
	Double			**OmegaDOverDCoeff;
	Double			**DThermConst;

// properties
	Double	*mixViscosity;
	Double	*mixDensity;
	Double	*mixConductivity;
	Double	*mixHeatCapacity;
	Double	*mixMolarMass;
	
// flame
	Double				**Y;
	Double				*temp;
	Double				*radiation;
	Double				*kappa; //Alexis 28/06/2018
	Double				*diffCorr;
	Double				viscosityInf;	// viscosity at right boundary
	Double				rhoInf;			// density at right boundary

// soot
	Double				**Pij;
	Double				*sumPi;
	Double				*pahMoments;
	Double				*moments;
	Double				*pahReactionRate;
//	Double				*sootReactionRate;
	Double				*diffSoot;
	Double				*sootSource;
	Double				**dMdx;
	
private:
	void InitTFlameNode( void ) {;};
};
typedef TFlameNode *TFlameNodePtr;

#endif // NODEINFO_H__
