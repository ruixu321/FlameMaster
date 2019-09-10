/*
	ListTool.h
*/

#ifndef __LTHeader__
#define __LTHeader__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#if defined (applec) || defined (powerc)
#include <Types.h>
#include <CursorCtl.h>
#include <Memory.h>
/*pascal unsigned long TickCount( void ) = 0xA975;*/
#endif

#define kMaxHeader			8192
#define streq(s1,s2)		(strcmp((s1),(s2))==0)
#ifndef FALSE
#define FALSE	0
#endif
#ifndef TRUE
#define TRUE	1
#endif

#define kMaxLenOfName		528
#define kMaxLenOfContext	528
#define kMaxLenOfString		800
#define kMaxLenOfUnit		244
#define kSafety				320
#define kTokenSeparator		'\n'
#define kTokenSepStr		"\n"
#define kMaxParameters		5000
#define kMaxTokens			10240
#define kMaxLabelChars		100240
#define kPool				40096	/* number of Doubles which are allocated each time */
#define kAtTrailer			(-2)	/* signals end of body when trailer is present */

typedef double			**doubleHdl, *doublePtr;

enum variableKind { kTextType, kPhysicalQuantity };

enum tokenKind {			/* list of allowed tokens */
	kSymbol,
	kNumber,
	kString,
	kAssignment,
	kLeftParen,
	kRightParen,
	kUnit,
	kBegin,
	kEnd,
	kHeader,
	kBody,
	kTrailer,
	kEndOfList
};

typedef struct _physicalQuantity {
	double value;
	char  unit[kMaxLenOfUnit];
} physicalQuantity;

typedef struct _parameter {
	char identifier[kMaxLenOfName];
	char context[kMaxLenOfContext];
	enum variableKind tag;
	union {
		char string[kMaxLenOfString];
		physicalQuantity quantity;
	} what;
} parameter;

typedef struct StartProfile {
	char		*labels;
	doublePtr	data;
	int			gridPoints;
	int			variables;
} StartProfile, *StartProfilePtr;


/*/typedef unsigned long	ulong;*/


/* Prototypes */

#ifdef __cplusplus
extern "C" {
#endif
/* ListTool.c:*/
StartProfilePtr ReadStartProfiles( StartProfilePtr sp, FILE *infp );
int GetStartProfileGridPoints( FILE *infp );
void CleanReadStartProfile(void);
parameter *GetParameter( const char *identifier );
void printHeader(char *fname);

/* LTHeader.c:*/
void parseHeader(char *stream, FILE *fp);

/* LTBody.c:*/
void parseBody(FILE *fp);
void saveData(FILE *fp);

/* LTUtilities.h:*/
double GetVariable( const char *varName );
char *lowercase(char *str);
void FatalError(const char *str);
void Warning(char *str);
void Message(char *str);
void printHelp(char *tname);
int moveto(char *str, FILE *fp);
int isReserved(char *str);
void strip(char *str);

/* LTSymbols.h*/
void checkSymbols(char *fn);
#ifdef __cplusplus
}
#endif

#endif /*__LTHeader__*/

#ifndef noErr
#define noErr 0
#endif
