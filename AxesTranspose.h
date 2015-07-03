#include<limits.h>
#include<stdlib.h>
#include<string.h>
#ifndef HLAXESTRANSPOSE_H_INCLUDED
#define HLAXESTRANSPOSE_H_INCLUDED

typedef unsigned int coord_t;
#define COORD_T_MAX UINT_MAX
#define MPI_COORD_T MPI_UNSIGNED

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

static inline hilpos_
GetHCoordinate (coord_t * Z, coord_t * X, coord_t * Y, int b, int n)
{				// position,bits,dimensions, Z,Y is used
	memcpy (X, Z, sizeof (coord_t) * n);
	coord_t M = 1 << (b - 1), P, Q, t;
	int i, j;
	for (Q = M; Q > 1; Q >>= 1)
	  {
		  P = Q - 1;
		  for (i = 0; i < n; i++)
			  if (X[i] & Q)
				  X[0] ^= P;
			  else
			    {
				    t = (X[0] ^ X[i]) & P;
				    X[0] ^= t;
				    X[i] ^= t;
			    }
	  }
	for (i = 1; i < n; i++)
		X[i] ^= X[i - 1];
	t = 0;
	for (Q = M; Q > 1; Q >>= 1)
		if (X[n - 1] & Q)
			t ^= Q - 1;
	for (i = 0; i < n; i++)
		X[i] ^= t;
	int nob = sizeof (coord_t) * 8;
	int wsk = 0;
	int wskbits = 0;
	Y[0] = 0;
	for (i = b - 1; i >= 0; i--)
	  {
		  for (j = 0; j < n; j++)
		    {
			    Y[wsk] *= 2;
			    Y[wsk] += ((X[j] >> i) & 1);
			    wskbits++;
			    if (wskbits == nob)
			      {
				      wskbits = 0;
				      wsk++;
				      Y[wsk] = 0;
			      }
		    }
	  }

	hilpos_t res = 0;
	hilpos_t akt_val = 1;
	for (i = n - 1; i >= 0; i--)
	  {
		  res += akt_val * Y[i];
		  akt_val *= (hilpos_t) (1 << b);
	  }

	return res;
}

#endif
