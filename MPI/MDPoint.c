#include "AxesTranspose.h"
#include "MDPoint.h"
#include <stdlib.h>
#include <string.h>


void
make_MDPoint (MDPoint * X, int dimensions)
{
	X->coordinates = calloc (dimensions, sizeof (coord_t));
	X->own_data_id = 0;
};

void
MDPointRemove (MDPoint * X)
{
	free (X->coordinates);
	X->own_data_id = 0;
}

int
MDPointPack (char *ptr, MDPoint * X, int dims)
{
	int i;
	memcpy (ptr, X->coordinates, sizeof (coord_t) * dims);
	ptr += dims * sizeof (coord_t);
	memcpy (ptr, &X->own_data_id, sizeof (tag_t));
	return sizeof (coord_t) * dims + sizeof (tag_t);
}

int
MDPointUnpack (char *ptr, MDPoint * Res, int dims)
{
	make_MDPoint (Res, dims);
	memcpy (Res->coordinates, ptr, sizeof (coord_t) * dims);
	ptr += dims * sizeof (coord_t);
	memcpy (&((Res)->own_data_id), ptr, sizeof (tag_t));
	return dims * sizeof (coord_t) + sizeof (tag_t);
}

int MDPointComparator(void *_a, void *_b, int dims) {
	MDPoint *a = (MDPoint*)_a;
	MDPoint *b = (MDPoint*)_b;
	int res = memcmp(a->coordinates,b->coordinates,sizeof(coord_t)*dims);
	if(res < 0)
		return 1;
	else
		return 0;
}



