/*
	List.h: Header file for List functions
*/

#ifndef __List__
#define __List__

#ifndef __alligator__
#include "alligator.h"
#endif

typedef struct Link {
	struct Link	*fNext;			/* pointer to next list item			*/
	struct Link	*fPrev;			/* pointer to previous item				*/
	void			*fItem;			/* pointer to data object				*/
} Link, *LinkPtr;

typedef struct List {
	LinkPtr	fHead;				/* pointer to the "head" of the list	*/
	LinkPtr	fTail;
	int		fNumItems;			/* current number of items in list		*/
	char	fLabel[40];			/* optional label						*/
} List, *ListPtr;

enum { kListHead, kListTail };


/* List.c */
ListPtr NewList(const char *);
LinkPtr NewLink( LinkPtr prev, LinkPtr next, void *item );
void AddItem(ListPtr, void *, int where);
void RemoveItem(ListPtr, void *);
void RemoveList(ListPtr);
void DeleteList(ListPtr);
void ForAll(ListPtr, void (*DoIt )(void *, void *), void *, int start);
void *FindItem(ListPtr, void *, int (*cmp )(void *, void *));
int ItemsInList(ListPtr);
ListPtr CloneList(ListPtr, const char *new_label);
void DisposeList(ListPtr);		/* currently equivalent to RemoveList	*/

#endif	/* __List__ */
