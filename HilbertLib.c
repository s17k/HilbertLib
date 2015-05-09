#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include "AxesTranspose.c"


struct MDPoint {
	coord_t* coords;
	void* own_data_ptr;
};
typedef struct MDPoint MDPoint;

void* home_ptr;
HilbertLibPosition* HilbertPos;
int comp(const void *elem1, const void *elem2) {
	int pos = elem1-home_ptr;
	int pos2 = elem2-home_ptr;
	return HilbertLibPositionComparator(HilbertLibPosition[pos], HilbertLibPosition[pos2]);
}

// X - data, ds - |data|, b - bits, n - |dimensions|
HilbertLibPosition* HilbertLibCurveSort(MDPoint *X, MDPoint* *SortedData, int ds, int b, int n) { 	
	HilbertLibPosition* res = calloc(ds,sizeof(HilbertLibPosition));
	int i=0;
	for(i=0;i<ds;i++) {
		res[i] = makeHilbertLibPosition(n);
		memcpy(res[i],X[i].coords,sizeof(coord_t)*n);
		AxestoTranspose(res[i],b,n);
	}
	void* ptrs = calloc(ds,sizeof(void*));
	for(i=0;i<ds;i++) 
		ptrs[i] = &X[i];
	home_ptr = (void*) X;
	HilbertPos = res;
	qsort(ptrs,ds,sizeof(*ptrs),comp);
	*SortedData = calloc(ds,sizeof(MDPoint));
	for(i=0;i<ds;i++) {
		(*SortedData)[i] = *(ptrs[i]);
	}
	free(X);
	free(ptrs);
	return res;
}

coord* HilbertLibMakeBins(MDPoint *X, MDPoint) {}

int HilbertLibBinSearch(HilbertLibPosition *X, int ds, HilbertLibPosition right) { // last <= right
	int left = 0, right = ds - 1, middle;
	while(HilbertLibPositionComparator()) {
		middle = left;
	}
}

int HilbertLibHowMany(HilbertLibPosition *X, int ds, HilbertLibPosition left, HilbertLibPosition right) {
	return HilbertLibBinSearch(X,ds,right) - HilbertLibBinSearch(X,ds,left);
}
