#include "PtrVector.h"
#include <stdlib.h>

void
makePtrVector (PtrVector * vec)
{
	vec->arr = NULL;
	vec->capacity = 0;
	vec->size = 0;
}


void
PtrVectorDoubleCapacity (PtrVector * vec)
{
	vec->arr =
		(void **) realloc (vec->arr,
				   (vec->capacity) *
				   2 * sizeof (void *));
	vec->capacity *= 2;
}

void
PtrVectorPB (PtrVector * vec, void *el)
{
	if (vec->capacity == 0)
	  {
		  vec->arr =
			  (void **) realloc (vec->arr,
					     1 * sizeof (void *));
		  vec->capacity = 1;
	  }
	if (vec->size + 1 > vec->capacity)
	  {
		  PtrVectorDoubleCapacity (vec);
	  }
	vec->arr[vec->size] = el;
	vec->size += 1;
}

void
PtrVectorDeallocate (PtrVector * vec)
{
	free (vec->arr);
	vec->capacity = 0;
	vec->size = 0;
}

typedef struct MTVector MTVector;
