#include "AxesTranspose.h"
#ifndef MDPOINT_H_INCLUDED
#define MDPOINT_H_INCLUDED
struct MDPoint {
        coord_t* coordinates;
        id_t own_data_id;
};

typedef struct MDPoint MDPoint;
void make_MDPoint(MDPoint *X, int dimensions);
void MDPointRemove(MDPoint *X);

#endif
