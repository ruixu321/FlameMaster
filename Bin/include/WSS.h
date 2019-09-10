
#ifndef __WSS__
#define __WSS__

#include "alligator.h"
#include "jgTypes.h"

typedef uint address;

typedef struct WSS {
	int			size;			/*	size of the stack (in bytes)				*/
	address		*base;			/*	base address of the stack					*/
	address		*top;			/*	top address of the stack (== base + size)	*/
	address		*cur;			/*	current stack pointer (top > cur >= base)	*/
} WSS, *WSSPtr;

WSSPtr NewWSS( int size );
void DisposeWSS( WSSPtr wss );
void *_WSSPush( WSSPtr wss, int size, const char *file, int line );
void WSSPop( WSSPtr wss );
int WSSRemove( WSSPtr wss, void *to_be_freed );
int WSSFreeMem( WSSPtr wss );
void WSSTrace( WSSPtr wss, int print );

#define WSSPush( wss, n )	_WSSPush(wss,n,__FILE__,__LINE__)

#endif	/*	__WSS__	*/

