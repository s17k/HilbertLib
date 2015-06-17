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
void HilbertLibNodeCurveSort(
MDPoint *X, // particles represented as MDPoints (input)
MDPoint* *SortedData, //Sorted MDPoints to return (output)
coord_t* *HCoordinates, // Hilbert coordinates to return (output)
int Datasize, // |data| - number of particles in this node (input)
int b,// precision (input)
int n // dimensions (input)
) 
{ 	
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
	*HCoordinates = calloc(Datasize,sizeof(coord_t));
	//free(first_elem); // needed for qsort only (but returned)
	(*SortedData) = calloc(Datasize,sizeof(MDPoint));
	for(i=0;i<Datasize;i++) {
		(*SortedData)[i] = *(ptrs[i]);
		(*HCoordinates)[i] = first_elem[ptrs[i]-(&X[0])];

	}
	free(X);
	free(ptrs);
	free(res);
}

// how many particles have hcoordinates <= Right
int HilbertLibNodeBinSearch(
coord_t *HCoordinates, 
int Datasize, 
coord_t Right) 
{ 	
	// binary search left,right,middle
	int bsleft = 0, bsright = Datasize, bsmiddle;  	
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
int HilbertLibNodeHowMany(coord_t *HCoordinates, int Datasize, coord_t Left, coord_t Right) { 
	return 
		HilbertLibNodeBinSearch(HCoordinates,Datasize,Right) - 
		HilbertLibNodeBinSearch(HCoordinates,Datasize,Left);
}

struct BinsBox {
	coord_t* data;
	size_t size;
	size_t capacity;
};

typedef struct BinsBox BinsBox;

void BinsBoxDoubleCapacity(BinsBox *b) {
	b->data = (coord_t*)realloc((void*)b->data, b->capacity*2);
	b->capacity*=2;
}

void BinsBoxPushBack(BinsBox *b, coord_t x) {
	b->size++;
	if(b->size > b->capacity) 
		BinsBoxDoubleCapacity(b);
	b->data[b->size-1] = x;
}

coord_t BinsBoxGet(BinsBox *b, size_t index) {
	return b->data[index];
}

void BinsBoxDeallocate(BinsBox *b) {
	free(b->data);
	b->size = 0;
	b->capacity = 0;
}

typedef struct BinsBox BinsBox;

// divide points into bins[root] 

int HilbertLibGetNOfParticles(int ProcessCount, int PointCount, int RootRank) {
	int *recvbuf = calloc(ProcessCount,sizeof(int));
	int sendbuf = PointCount;
	MPI_Gather(&sendbuf,1,MPI_INT,recvbuf,1,MPI_INT,RootRank,MPI_COMM_WORLD);
	int suma = 0;
	int i;
	for(i=0;i<ProcessCount;i++) 
		suma += recvbuf[i];
	return suma;
}

coord_t HilbertLibCalculateNextBoundary(
coord_t a, coord_t b, 
coord_t *MyHCoordinates,
int MyPointsCount,
int particlesRate,
int nodesCount,
int RootRank
) {
	coord_t bsleft = a, bsright = b, bsmiddle;
	int* hm = calloc(nodesCount,sizeof(int));
	int i = 0;
	int suma;
	while(bsleft < bsright) {
		bsmiddle = (bsleft + bsright)/2;			
		coord_t * sendBuff = calloc(2,sizeof(coord_t));
		sendBuff[0] = a;
		sendBuff[1] = bsmiddle;
		coord_t * recvBuff = calloc(2,sizeof(coord_t));
		MPI_Scatter(sendBuff,2,MPI_UNSIGNED,recvBuff,2,MPI_UNSIGNED,RootRank,MPI_COMM_WORLD);
		int singlebuff = HilbertLibNodeHowMany(MyHCoordinates,MyPointsCount,sendBuff[0],sendBuff[1]);
		MPI_Gather(&singlebuff,1,MPI_INT,hm,1,MPI_INT,RootRank,MPI_COMM_WORLD);
		suma = 0;
		for(i=0;i<nodesCount;i++) {
			suma += hm[i];		
		}
		if(suma > particlesRate) {
			bsright = bsmiddle-1;
		} else {
			bsleft = bsmiddle;
		}
	}
	return bsmiddle;
}

void HilbertLibRootMakeBins(
int RootRank, // rank of root
size_t NodesCount, // number of nodes
coord_t *MyHCoordinates, // Hilbert Coordinates of root points
size_t MyPointsCount, // Number of points assigned to root
int b, // |precision bits|
int dimensions // |dimensions|
)
{
	int allPointsCount = HilbertLibGetNOfParticles(NodesCount,MyPointsCount,RootRank);
	int particlesRate = allPointsCount/NodesCount;
	int i = 0;
	// in binsBoundaries[i] there will be last Hilbert Coordinate for process i
	coord_t* binsBoundaries = calloc(NodesCount,sizeof(coord_t));
	coord_t lastused = 0;
	for(i=0;i<NodesCount;i++) {
		binsBoundaries[i] =
			HilbertLibCalculateNextBoundary(
				lastused,(1<<b)-1,
				MyHCoordinates,
				MyPointsCount,
				particlesRate,
				NodesCount,
				RootRank
			);
	}
}

void HilbertLibMPISendNOfParticles(int what, int RootRank) {
	int sendBuff = what;
	MPI_Gather(&sendBuff,1,MPI_INT,NULL,NULL,RootRank,MPI_COMM_WORLD);
}

//divide points into bins[node]
HilbertLibNodeMakeBins(coord_t *MyHCoordinates, size_t MyParticlesCount, int RootRank) {
	//Send number of particles i have
	HilbertLibMPISendNoOfParticles(MyParticlesCount);	
	void * recvBuff = (void*) calloc(2,sizeof(coord_t));
	while(true) {
		MPI_Scatter(NULL,NULL,NULL,recvBuff,2,MPI_INT,RootRank,MPI_COMM_WORLD);
		if(recvBuff[0] < 0) {
			break;
		}
		int sendbuf = HilbertLibNodeHowMany(MyHCoordinates,MyParticlesCount,recvBuff[0],recvBuff[1]);
		MPI_Gather(&sendbuf,1,MPI_INT,NULL,NULL,RootRank,MPI_COMM_WORLD);		
	}
}


