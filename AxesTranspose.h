#ifndef HLAXESTRANSPOSE_H_INCLUDED
#define HLAXESTRANSPOSE_H_INCLUDED

typedef unsigned int coord_t;
#define MPI_COORD_T MPI_UNSIGNED_LONG_LONG

//unsigned int hilpos_t
/*
typedef unsigned int hilpos_t;
#define MPI_HILPOS_T MPI_UNSIGNED
#define HILPOS_T_MAX_VALUE ((hilpos_t)1<<31)
#define HILPOS_EPS 1
//change to 2.0
#define HILPOS_BS_2 2
//change to 1.0
#define HILPOS_BS_1 1 
//change to HILPOS_EPS
#define HILPOS_MARGIN 0*/

//double hilpos_t
typedef double hilpos_t;
#define MPI_HILPOS_T MPI_DOUBLE
#define HILPOS_T_MAX_VALUE ((hilpos_t)1e200)
#define HILPOS_EPS ((hilpos_t)(1e-15))
#define HILPOS_BS_1 ((hilpos_t)0)
#define HILPOS_BS_2 ((hilpos_t)2)
#define HILPOS_MARGIN ((hilpos_t)(1e-1))


typedef unsigned int tag_t;
#define MPI_TAG_T MPI_UNSIGNED
typedef unsigned int amount_t;

hilpos_t GetHCoordinate( coord_t *, coord_t *, coord_t*, int, int);

#endif
