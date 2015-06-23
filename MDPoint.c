#include "AxesTranspose.h"
#include "MDPoint.h"
#include <stdlib.h>

void make_MDPoint(MDPoint *X, int dimensions) {
	X->coordinates = calloc(dimensions, sizeof(coord_t));
	X->own_data_id = 0;
};

void MDPointRemove(MDPoint *X) {
	free(X->coordinates);
	X->own_data_id = 0;
}


