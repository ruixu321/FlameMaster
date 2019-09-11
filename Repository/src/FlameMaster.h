#ifndef __FLAMEMASTER_H__
#define __FLAMEMASTER_H__

/* This header is convenience header to gather all configurations
 *
 */

/* to add a new flametype one has to add the class with all initializers
 * for constants and pure virtual functions, a constructor, a destructor and an
 * initializer the flametype has to be added to the enumeration list 'FlameType' 
 * (in Input.h) and an identifier for the input specification has to be introduced 
 * in file 'ReadAddData.flex' in the startcondition 'scReadFlameType' flame specific
 * initialization of the class FirstInput has to be done in function
 * 'FirstInput::InitFirstInput' and for the class TInputData in
 * 'TInputData::InitInputData' the print functions of FirstInput and TInputData
 * should be updated
 */

#include "SteadyStates.h"
#include "TTransportModel.h"

#include "T0DIsoChor.h"
#include "T0DPSR.h"  //Krithika
#include "TCountDiffFlameCont.h"
#include "TCountDiffFlameMix.h"
#include "TCountDiffFlamePhys.h"
#include "TCountDiffFlameSim.h"
#include "TCountDiffPhysEigen.h"
#include "TCountPremFlamePhys.h"
#include "TCountPremFlameSim.h"
#include "TFlameSheet.h"
#include "TTransFlamelet.h"
#include "TUnstrPremFlamePhys.h"
#ifdef COMPILE_FORTRAN_SRC
#include "TEquilibrium.h"
#include "TTrans1DIsoChor.h"
#endif  //  COMPILE_FORTRAN_SRC

#endif  // __FLAMEMASTER_H__
