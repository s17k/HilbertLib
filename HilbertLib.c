#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include "AxesTranspose.c"


struct MDPoint {
	coord_t* coordinates;
	void* own_data_ptr;
};
typedef struct MDPoint MDPoint;

void make_MDPoint(MDPoint *X, int dimensions) {
	X->coordinates = calloc(dimensions, sizeof(coord_t));
	X->own_data_ptr = NULL;
};

MDPoint* HomePtr;
coord_t *HilbertPos;
int HilbertLibCurveSortComparator(const void *_elem1, const void *_elem2) {
	MDPoint* *__elem1 = (MDPoint**)_elem1;
	MDPoint* *__elem2 = (MDPoint**)_elem2;
	MDPoint* elem1 = *__elem1;
	MDPoint* elem2 = *__elem2;

	int pos = elem1-HomePtr;
	int pos2 = elem2-HomePtr;
	if(HilbertPos[pos] < HilbertPos[pos2])
		return 1;
	else
		return 0;
}

// X - data, Datasize - |data|, b - bits, n - |dimensions|
void HilbertLibNodeCurveSort(MDPoint *X, MDPoint *SortedData, coord_t *HCoordinates, int Datasize, int b, int n) { 	
	MDPoint* res = calloc(Datasize,sizeof(MDPoint)); // miejsce do operowania AxestoTranspose
	int i=0;
	for(i=0;i<Datasize;i++) { // saving i-th Hilbert Coordinates in res[i].coords[0]
		make_MDPoint(&res[i], n);
		memcpy(res[i].coordinates,X[i].coordinates,sizeof(coord_t)*n);
		AxestoTranspose(res[i].coordinates,b,n);
	}
	MDPoint* *ptrs = calloc(Datasize,sizeof(MDPoint*));
	for(i=0;i<Datasize;i++) 
		ptrs[i] = &X[i];
	HomePtr = X;
	coord_t* first_elem = calloc(Datasize,sizeof(coord_t));
	for(i=0;i<Datasize;i++) 
		first_elem[i] = res[i].coordinates[0];
	HilbertPos = first_elem;
	qsort(ptrs,Datasize,sizeof(MDPoint*),HilbertLibCurveSortComparator);
	HCoordinates = first_elem;
	//free(first_elem); // needed for qsort only (but returned)
	SortedData = calloc(Datasize,sizeof(MDPoint));
	for(i=0;i<Datasize;i++) {
		SortedData[i] = *(ptrs[i]);
	}
	free(X);
	free(ptrs);
	free(res);
}

//coord_t* HilbertLibMakeBins(MDPoint *X, MDPoint) {}

int HilbertLibRootBinSearch(// how many nodes have hcoordinates <= Right
coord_t *HCoordinates, 
int Datasize, 
coord_t Right) 
{ 	
	int bsleft = 0, bsright = Datasize, bsmiddle;  // binary search left,right,middle
	// properly it is binsearching the first node which have hcoordinate > Right
	while(bsleft < bsright) {	
		bsmiddle = (bsleft + bsright)/2;
		if(HCoordinates[bsmiddle] <= Right) {
			bsleft = bsmiddle+1;
		} else {
			bsright = bsmiddle;
		}
	}
	return bsleft;
}

// counting all nodes(hc) which satisfy : Left < hc <= Right
int HilbertLibRootHowMany(coord_t *HCoordinates, int Datasize, coord_t Left, coord_t Right) { 
	return 
		HilbertLibRootBinSearch(HCoordinates,Datasize,Right) - 
		HilbertLibRootBinSearch(HCoordinates,Datasize,Left);
}

