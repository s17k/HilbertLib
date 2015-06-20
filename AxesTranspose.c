#include<string.h>
#include<stdlib.h>
typedef unsigned int coord_t;
#define MPI_COORD_T MPI_UNSIGNED 
#define COORD_T_MAX_VALUE ((coord_t)1<<31)
typedef unsigned int hilpos_t;
#define MPI_HILPOS_T MPI_UNSIGNED
#define HILPOS_T_MAX_VALUE ((hilpos_t)1<<31)
#define HILPOS_EPS 1
typedef unsigned int id_t;
#define MPI_ID_T MPI_UNSIGNED
typedef unsigned int amount_t;

// coordinates(X) ---> HilbertLibPosition(X) (changes X from coordinates in n dimensions to HilbertLibPosition)
void AxestoTranspose( coord_t* X, int b, int n) { // position,bits,dimensions
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
	coord_t *Y = calloc(n,sizeof(coord_t));
	int nob = sizeof(coord_t)*8;
	int wsk = 0;
	int wskbits =0;
	for(i=b-1;i>=0;i--) {
		for(j=0;j<n;j++) {
			Y[wsk] *= 2;
			Y[wsk] += ((X[j]>>i)&1);
			wskbits++;
			if(wskbits == nob) {
				wskbits = 0;
				wsk++;
			}
		}
	}
	memcpy(X,Y,sizeof(Y));
	free(Y);
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
