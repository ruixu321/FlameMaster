#ifndef __MapMan__
#define __MapMan__

using namespace std; // GB

/*
	MapMan.h:	Header file for the MapMan package.

	MapMan is a package for manipulation of collections of data sets.
	Supported tasks are:
		o	shift profiles
		o	scale profiles
		....
*/


//		Constants
//--------------------------------------------------------------------------------
const int	kStringLen = 24;


//		Forward declarations
//--------------------------------------------------------------------------------
class MMDataBag;


//		Class MMDataSet
//--------------------------------------------------------------------------------

class MMDataSet
{
	friend ostream& operator<<( ostream&, MMDataSet& );		

public:
	enum XType { zeroValue, lastValue, zeroGradient, lastGradient };
	
					MMDataSet();
					MMDataSet( Double*, int, int = 1, const char* = 0, 
						MMDataBag* = 0 );
					MMDataSet( VectorPtr, const char* = 0, MMDataBag* = 0 );
					//~MMDataSet();

	void			FreeMemory( void );
	int				IsValid( void ) const;
	int				IsSameDataSet( const Double* v, int len, int off, 
						const char* name ) const;
	int				Length( void ) const;
	int				Offset( void ) const;
	const char*	Name( void ) const;	
	void			SetXType( XType );
	XType			GetXType( void );
	
	//	Mapping functions.
	virtual void	Map( Double* y_new, int len, int offset=1 );
	virtual void	Map( VectorPtr );
	virtual VectorPtr	Map( void );
	
	
private:
	Double*		fVec;
	int				fLen;
	int				fOffset;
	char			fName[kStringLen];
	XType			fExtrapolate;	// one field enough?
	MMDataBag*		fBag;			// we need a reference to its bag
	
	const Double*	GetData( void ) const;
};


inline int MMDataSet::IsValid( void ) const {
	return (fVec != NULL) && (fLen > 0);
}

inline const char* MMDataSet::Name( void ) const {
	return (const char*)fName;
}

inline int MMDataSet::Length( void ) const { return fLen; }
inline int MMDataSet::Offset( void ) const { return fOffset; }


//		Class MMDataBag
//--------------------------------------------------------------------------------

class MMDataBag
{
	friend ostream& operator<<( ostream&, MMDataBag& );		

public:
						MMDataBag( int maxItems );
	virtual				~MMDataBag();
	void				Initialize( void );
	MMDataSet&			operator[]( int index );
	MMDataSet&			operator[]( const char* name);	// not yet implemented
	int					Insert( Double* v, int len, int offset = 1, const char *name = 0 );
	int					Insert( VectorPtr v, const char *name = 0 );
						// returns 0 if successful, -1 otherwise

	//	Old and new independent variables.
	virtual void		SetNewInpedVar( Double *, int, int off = 1, 
							const char *name = 0 );
	virtual void		SetNewInpedVar( VectorPtr, const char *name = 0 );
	virtual void		SetNewInpedVar( MMDataSet& );
	virtual void		SetOldInpedVar( Double *, int, int off = 1, 
							const char *name = 0 );
	virtual void		SetOldInpedVar( VectorPtr, const char *name = 0 );
	virtual void		SetOldInpedVar( MMDataSet& );
	const MMDataSet* 	GetNewInpedVar( void ) const;
	const MMDataSet* 	GetOldInpedVar( void ) const;

	virtual void		FreeMemory( void );
	int					NumElems( void ) { return fCurrent; }
	int					Length( void ) const;

private:
	int					fCurrent;
	int					fElems;
	int					fArrayLength;	// all arrays must have the same length
	MMDataSet*			fPool;
	MMDataSet			fNewIndepVar;
	MMDataSet			fOldIndepVar;
};

inline int MMDataBag::Length( void ) const
	{ return fPool[0].Length(); }
inline const MMDataSet* MMDataBag::GetNewInpedVar( void ) const
	{ return (const MMDataSet*)&fNewIndepVar; }
inline const MMDataSet* MMDataBag::GetOldInpedVar( void ) const
	{ return (const MMDataSet*)&fOldIndepVar; }


//		Class NodeMover
//--------------------------------------------------------------------------------

class NodeMover
{
	friend ostream& operator<<( ostream&, NodeMover& );		

public:
	enum fromType { fromLeftSide, fromRightSide };

						NodeMover();
						NodeMover( Double* oldVar, int len );
						NodeMover( VectorPtr oldVar );
						NodeMover( Double* oldVar, Double* newVar, int len, 
							int n, fromType from );
						NodeMover( VectorPtr oldVar, VectorPtr newVar, 
							int n, fromType from );
												
	virtual int		MoveIt( void );
	virtual int		MoveIt( Double* newVar, int len, int n, fromType from );
	virtual int		MoveIt( VectorPtr newVar, int n, fromType from );
	
private:
	int					fLen;			// old and new have same length
	int					fNumToMove;		// number of nodes to move
	fromType			fFrom;			// move from left of right edge
	Double*			fOld;			// old array
	Double*			fNew;			// new array
	
	void				SetOldVar( Double* oldVar, int len );
	void				SetNewVar( Double* newVar, int len );
	int					BadMembers( void );
	virtual int		MoveNodesFromLeft( void );
	virtual int		MoveNodesFromRight( void );
};

#endif	/* __MapMan__ */
