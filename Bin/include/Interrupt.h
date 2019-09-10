#ifndef __Interrupt__
#define __Interrupt__

enum Interrupt	{ kNoSignal, kExit, kNewGrid };

#if defined (applec) || defined (powerc)
#	pragma once
#endif

#ifdef __cplusplus
extern "C" {
#endif
	extern int gExit;
	void InstallInterruptHP( void );
#ifdef __cplusplus
}
#endif

#endif	/* __Interrupt__ */

