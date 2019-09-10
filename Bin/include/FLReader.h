/*
	FLReader.h
	Header file for Flamelet Reader package

	 Josef Gttgens, Peter Terhoeven, Ian Herwono, 1993
	
	Version 1.b2		07/21/93
*/

#include "Config.h"

#if defined COMP_MSVC
//#include "../../WindowsPort/VisualStudio.h"
#include "../../WindowsPort/wingetopt.h"
#endif

#ifndef __FLReader__
#define __FLReader__

#ifdef applec
#pragma once
#endif

enum { 				/* used for type determination */
	kFLRBoolean, 
	kFLRInteger, 
	kFLRFloat, 
	kFLRString, 
	kFLRArray
};

enum {		/* used to specify the section range when writing a flamelet file */
	kFLRAll, 			/*  kFLRAll 	: all sections 		*/
	kFLRHeader, 		/*  kFLRHeader 	: Header section 	*/
	kFLRBody, 			/*  kFLRBody 	: Body section 		*/
	kFLRTrailer			/*  kFLRTrailer	: Trailer section  	*/
};


typedef struct FLRSymbol {
	int 	type;
	union 	{  int 		bool;
			   int 		i;
			   Double 	f;
			   char 	*s;
			   VectorPtr v;
			  } val;
	char	*id;
	char	*unit;
	char	*context;
	int		detached;
} FLRSymbol,*FLRSymbolPtr;


/*	Global lists													*/
extern ListPtr
	gFLRSymbolList,					/*  global FLRSymbol list 		*/
	gFLRNumberList, 				/*  global Number list 			*/
	gFLRIntList, 					/*  global Integer list 		*/
	gFLRFloatList, 					/*  global Float list 			*/
	gFLRStringList, 				/*  global String list 			*/
	gFLRArrayList;					/*  global Array (Vector) list 	*/

#ifdef __cplusplus
extern "C" {
#endif

void InitFLReader(void);
	/*	initializes the FLReader package								*/
int ReadFlamelet(FILE *fp);
	/*	reads the flamelet file "fname", all global lists are valid
		after reading the flamelet file									*/
void WriteFlameletFile( FILE *fp,int section );
	/*	writes output file in flamelet format							*/
void CleanupFLReader(void);
	/*	cleans up the storage associated with the FLReader Package		*/
void DetachItem(void *ptr, void *);
	/*	set detached = TRUE, so this FLRSymbol's value will 
		not be deleted when calling CleanupFLReader						*/
void PrintSymbolList( FILE *fp);
	/*	prints out the global FLRSymbol list							*/
ListPtr GetPatternedItems(char *pattern);
	/*	collects the item(s) which match the pattern and 
		returns a list of matched items	 								*/

int ChangeString(const char *all,char *newstring);
	/*	changes string value to newstring (all: reg. expr. pat.),
		returns TRUE if found											*/
int ChangeNumber(const char *all,const Double val, char *unit);
	/*	changes number value to val (all: reg. expr. pat.),
		returns TRUE if found											*/
int ChangeArray(const char *all,VectorPtr val,char *unit);
	/*	changes array to val (all: reg. expr. pat.),
		returns TRUE if found											*/
int DeleteSymbol(const char *all);
	/*	deletes a FLRSymbol item  (all: reg. expr. pat.),
		returns TRUE if found											*/

void PrintSymbol( void *ptr, void *aux );
	/*	prints out the FLRSymbol to fp (ForAll)							*/
int GetString(const char *all,char *val);	
	/*  gets string value in val (all: reg. epxr. pat.), 
		returns TRUE if found											*/
int GetNumber( const char *all, Double *val, char *unit );
	/*  gets number value in val (all: reg. epxr. pat.), 
		returns TRUE if found											*/
int GetInteger( const char *all, int *val, char *unit );
	/*  gets integer value in val (all: reg. epxr. pat.), 
		returns TRUE if found											*/
int GetFloat(const char *all,Double *val,char *unit);
	/*  gets float value for in val (all: reg. epxr. pat.), 
		returns TRUE if found											*/
int GetArray(const char *all,VectorPtr *val,char *unit);
	/*  gets Array in val (all: reg. epxr. pat.), 
		returns TRUE if found											*/
int FLRCompare( const char *ptr1, const char *ptr2 );
	/*	case insensitive string comparison								*/
	
#ifdef __cplusplus
}
#endif


#endif	/* __FLReader__ */
