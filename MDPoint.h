#include "AxesTranspose.h"
#ifndef MDPOINT_H_INCLUDED
#define MDPOINT_H_INCLUDED
#define MDPOINT_NETTO_SIZE (sizeof(coord_t*) + sizeof(tag_t))
#define MDPOINT_FULL_SIZE(D) (sizeof(coord_t)*D + sizeof(tag_t))

struct MDPoint
{
	coord_t *coordinates;
	tag_t own_data_id;
};

typedef struct MDPoint MDPoint;
void make_MDPoint (MDPoint * X, int dimensions);
void MDPointRemove (MDPoint * X);
int MDPointPack (char *ptr, MDPoint * X,
		 int dims);
int MDPointUnpack (char *ptr, MDPoint * Res,
		   int dims);

#endif
