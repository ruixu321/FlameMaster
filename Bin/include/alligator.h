/*
	alligator.h: header file for the "simple alligator" allocation package
	
	History:
		24.01.94	general revision (jg)
		27.12.92	1st version (jg)
*/

#ifndef __alligator__
#define __alligator__

#ifdef __cplusplus
extern "C" {
#endif
extern void *allocate( unsigned int size, const char *file, int line );
extern void fatalerrormsg(const char *s, const char *file, int line );
extern void alligator_delete( void *ptr );
#ifdef __cplusplus
}
#endif

/*	NEW(t)  allocates memory for 1 data structure of type t.
	Example:	float *fPtr = NEW(float);
*/
#define	NEW(t)		((t*)allocate(sizeof(t),__FILE__,__LINE__))

/*	NEW2(n,t)  allocates memory for an array of n elements of
	a data structure of type t.
	Example:	float *fPtr = NEW2(10,float);
*/
#define	NEW2(n,t)	((t*)allocate((unsigned)((n)*sizeof(t)),__FILE__,__LINE__))

/*	DELETE(ptr)  deallocates a data object previously allocated
	by the NEW or NEW2 macro.  ptr may be the NULL pointer.
*/
#define DELETE(ptr)	alligator_delete( ptr )

#define FATAL(s)	fatalerrormsg((s),__FILE__,__LINE__)

#endif	/* __alligator__ */

