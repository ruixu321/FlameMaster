/*
	TofZ.h: header file for TofZ.c
	
	Revision history:
		7/22/91		first version by jg
		7/28/91		implementation of the derivative function (jg)
					wrong slope at Z = § corrected (jg)
*/

#ifndef __TofZ__
#define __TofZ__

typedef struct TofZParams {
	Double	Tu,				/* unburnt temperature in [K]		*/
			Tmax,			/* maximum temperature in [K]		*/
			Zst,			/* stoichiometric mixture fraction	*/
			alpha,			/* left crossover value of Z		*/
			beta;			/* right crossover value of Z		*/
} TofZParams, *TofZParamsPtr;


#ifdef __cplusplus
extern "C" {
#endif

Double temperature( Double z );
Double dTdZ( Double z );
void set_coefficients( TofZParamsPtr ptr );
void get_coefficients( TofZParamsPtr ptr );

#ifdef __cplusplus
}
#endif


#endif	/* __TofZ__ */
