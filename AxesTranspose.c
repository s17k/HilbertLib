#include<string.h>
#include<stdlib.h>
#include<stdio.h>
typedef unsigned int coord_t;
#define MPI_COORD_T MPI_UNSIGNED 
#define COORD_T_MAX_VALUE ((coord_t)1<<31)

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


typedef unsigned int id_t;
#define MPI_ID_T MPI_UNSIGNED
typedef unsigned int amount_t;

// coordinates(X) ---> HilbertLibPosition(X) (changes X from coordinates in n dimensions to HilbertLibPosition)
hilpos_t GetHCoordinate( coord_t* X, coord_t *Y, int b, int n) { // position,bits,dimensions, X is destroyed, Y is used
	coord_t M = 1 << (b-1), P, Q, t;
	int i,j;
	for(Q=M;Q>1;Q>>=1) {
		P = Q -1;
		for(i=0;i<n;i++)
			if(X[i]&Q) 
				X[0]^=P;
			else 
				{t = (X[0]^X[i]) & P;X[0] ^= t; X[i] ^= t; }
	}
	for( i=1;i<n;i++) 
		X[i] ^= X[i-1];
	t = 0;
	for(Q=M;Q>1;Q>>=1)
		if(X[n-1] & Q)
			t ^= Q-1;
	for(i=0;i<n;i++) X[i]^=t;
	//coord_t *Y = calloc(n,sizeof(coord_t));
	int nob = sizeof(coord_t)*8;
	int wsk = 0;
	int wskbits =0;
	Y[0] = 0;
	for(i=b-1;i>=0;i--) {
		for(j=0;j<n;j++) {
			Y[wsk] *= 2;
			Y[wsk] += ((X[j]>>i)&1);
			wskbits++;
			if(wskbits == nob) {
				wskbits = 0;
				wsk++;
				Y[wsk] = 0;
			}
		}
	}
	
	hilpos_t res = 0;
	hilpos_t akt_val = 1;
	for(i=n-1;i>=0;i--) {
		res += akt_val * Y[i];
		akt_val *= (hilpos_t)(1<<b);
	}

	return res;
	//free(Y);
}

/*typedef HilbertLibPosition coord_t*;
int HilbertLibPositionComparator(const HilbertLibPosition *A, const HilbertLibPosition *B, int size) { // maybe inlined
	int r = memcmp(*A,*B,size);
	if(r <= 0) return 1;
	else return 0;
}
HilbertLibPosition* makeHilbertLibPosition(int dimensions) {
	HilbertLibPosition* res = calloc(dimensions,sizeof(coord_t);
	return res;
}
void HilbertLibPositionFree(HilbertLibPosition *X) {
	free(X);
}*/
