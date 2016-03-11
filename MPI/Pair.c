#include "Pair.h"
void
make_Pair (Pair * t, int x, int y)
{
	t->a = x;
	t->b = y;
}

int
PairComparator (const void *_elem1, const void *_elem2)
{
	Pair *elem1 = (Pair *) _elem1;
	Pair *elem2 = (Pair *) _elem2;
	if (elem1->a == elem2->a)
	  {
		  if (elem1->b < elem2->b)
			  return 0;
		  else
			  return 1;
	  }

	if (elem1->a < elem2->a)
		return 0;
	else
		return 1;
}
