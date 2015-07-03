#ifndef _PAIR_H_DEFINED_
#define _PAIR_H_DEFINED_
struct Pair
{
	int a, b;
};
typedef struct Pair Pair;

void make_Pair (Pair * t, int x, int y);

int PairComparator (const void *_elem1,
		    const void *_elem2);

#endif
