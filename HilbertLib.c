#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>
#include "AxesTranspose.c"
#include "BinsBox.c"


struct MDPoint {
	coord_t* coordinates;
	void* own_data_ptr;
};


typedef struct MDPoint MDPoint;

void make_MDPoint(MDPoint *X, int dimensions) {
	X->coordinates = calloc(dimensions, sizeof(coord_t));
	X->own_data_ptr = NULL;
};

void MDPointRemove(MDPoint *X) {
	free(X->coordinates);
	free(X->own_data_ptr);
}

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
		return 0;
	else
		return 1;
}

// X - data, Datasize - |data|, b - bits, n - |dimensions| [Node]
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
		printf("Znajduje H dla %d,",i);
		AxestoTranspose(res[i].coordinates,b,n);
		printf(" Wyszlo mi %d\n",res[i].coordinates[0]);
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
	//free(X);
	free(ptrs);
	free(res);
}

// how many particles have hcoordinates <= Right [Node]
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

// counting all nodes(hc) which satisfy : Left < hc <= Right [Node]
int HilbertLibNodeHowMany(coord_t *HCoordinates, int Datasize, coord_t Left, coord_t Right) { 
	return 
		HilbertLibNodeBinSearch(HCoordinates,Datasize,Right) - 
		HilbertLibNodeBinSearch(HCoordinates,Datasize,Left);
}


// divide points into bins[Root] 

int HilbertLibGetNOfParticles(int ProcessCount, int PointCount, int RootRank) {
	int *recvbuf = calloc(ProcessCount,sizeof(int));
	int sendbuf = PointCount;
	MPI_Gather(&sendbuf,1,MPI_INT,recvbuf,1,MPI_INT,RootRank,MPI_COMM_WORLD);
	int suma = 0;
	int i;
	for(i=0;i<ProcessCount;i++) 
		suma += recvbuf[i];
	return suma;
	free(recvbuf);
}

// Calculates the next boundary [Root]
coord_t HilbertLibCalculateNextBoundary( 
coord_t a, coord_t b, 
coord_t *MyHCoordinates,
int MyPointsCount,
int particlesRate,
int nodesCount,
int RootRank,
int* how_many_used // how many used is not filled
) {
	coord_t bsleft = a, bsright = /*b*/1000000000, bsmiddle;
	int* hm = calloc(nodesCount,sizeof(int));
	int i = 0;
	int suma;
	while(bsleft < bsright) {
		printf("bsleft : %d, bsright %d\n",bsleft,bsright);
		bsmiddle = (bsleft + bsright + 1)/2;			
		coord_t * sendBuff = calloc(2,sizeof(coord_t));
		sendBuff[0] = a;
		sendBuff[1] = bsmiddle;
		printf("Wysyłam parę (%d,%d)\n",a,bsmiddle);
		coord_t * recvBuff = calloc(2,sizeof(coord_t));
		MPI_Scatter(sendBuff,2,MPI_COORD_T,recvBuff,2,MPI_COORD_T,RootRank,MPI_COMM_WORLD);
		int singlebuff = HilbertLibNodeHowMany(
			MyHCoordinates,
			MyPointsCount,
			sendBuff[0],
			sendBuff[1]
		);
		free(sendBuff);
		free(recvBuff);

		MPI_Gather(&singlebuff,1,MPI_INT,hm,1,MPI_INT,RootRank,MPI_COMM_WORLD);
		suma = 0;
		for(i=0;i<nodesCount;i++) {
			suma += hm[i];		
		}
		printf("na tym przedziale jest %d\n",suma);
		if(suma > particlesRate) {
			bsright = bsmiddle-1;
		} else {
			bsleft = bsmiddle;
		}
	}
	free(hm);
	*how_many_used = 0;
	return bsmiddle;
}

// Make Bins [Root]
coord_t* HilbertLibRootMakeBins(
int RootRank, // rank of root
size_t NodesCount, // number of nodes
coord_t *MyHCoordinates, // Hilbert Coordinates of root points
size_t MyPointsCount, // Number of points assigned to root
int b, // |precision bits|
int dimensions // |dimensions|
)
{
	int allPointsCount = 
		HilbertLibGetNOfParticles(
			NodesCount,
			MyPointsCount,
			RootRank
		);
	printf("AllPointsCount : %d\n",allPointsCount);
	int particlesRate = allPointsCount/NodesCount;
	int i = 0;
	// in binsBoundaries[i] there will be last Hilbert Coordinate for process i
	coord_t* binsBoundaries = calloc(NodesCount,sizeof(coord_t));
	coord_t lastused = 0;
	for(i=0;i<NodesCount;i++) {
		int how_many_used = 0;
		printf("Licze Boundary dla %d ...   \n ", i);
		binsBoundaries[i] =
			HilbertLibCalculateNextBoundary(
				lastused,(1<<b)-1,
				MyHCoordinates,
				MyPointsCount,
				particlesRate,
				NodesCount,
				RootRank,
				&how_many_used
			);
		printf("I wyszlo mi %d ... \n ", binsBoundaries[i]);
		allPointsCount -= how_many_used;
		// -i-1 commented, becauase how_many_used is not updated
		particlesRate = allPointsCount/(NodesCount/*-i-1*/); 
	}
	// last scatter to indicate the end of queries
	int* sendbuf = calloc(2,sizeof(int));
	sendbuf[0] = sendbuf[1] = -1;
	int* recvbuf = calloc(2,sizeof(int));
	MPI_Scatter(
		sendbuf,
		2,
		MPI_INT,
		recvbuf,
		2,
		MPI_INT,
		RootRank,
		MPI_COMM_WORLD
	);
	free(sendbuf);
	free(recvbuf);
	return binsBoundaries;
}

//Sends # of particles to Root [Node]
void HilbertLibSendNOfParticles(int what, int RootRank) {
	int sendBuff = what;
	MPI_Gather(&sendBuff,1,MPI_INT,NULL,0,NULL,RootRank,MPI_COMM_WORLD);
}

//divide points into bins[Node]
void HilbertLibNodeMakeBins(coord_t *MyHCoordinates, size_t MyParticlesCount, int RootRank) {
	//Send number of particles i have
	HilbertLibSendNOfParticles(MyParticlesCount, RootRank);	
	coord_t *recvBuff = calloc(2,sizeof(coord_t));
	while(true) {
		MPI_Scatter(
			NULL,
			0,
			NULL,
			recvBuff,
			2,
			MPI_COORD_T,
			RootRank,
			MPI_COMM_WORLD
		);
		if(recvBuff[0] < 0) {
			break;
		}
		printf("Dostałem %d %d\n",recvBuff[0],recvBuff[1]);
		int sendbuf = HilbertLibNodeHowMany(
			MyHCoordinates,
			MyParticlesCount,
			recvBuff[0],
			recvBuff[1]
		);
		printf("wiec odpowiadam %d\n",sendbuf);
		MPI_Gather(
			&sendbuf,
			1,
			MPI_INT,
			NULL,
			0,
			NULL,
			RootRank,
			MPI_COMM_WORLD
		);		
	}
}

// [Root]
void HilbertLibSendBoundariesToAll(int* Boundaries, int Size, int RootRank) {
	coord_t* recvBuff = calloc(Size,sizeof(coord_t));
	MPI_Scatter(
		Boundaries,
		Size,
		MPI_COORD_T,
		recvBuff,
		Size,
		MPI_COORD_T,
		RootRank,
		MPI_COMM_WORLD
	);
	free(recvBuff);
}

//[Node]
int * HilbertLibRecvBoundariesFromRoot(int Size, int RootRank) {
	coord_t* recvBuff = calloc(Size,sizeof(coord_t));
	MPI_Scatter(
		NULL,
		0,
		0,
		recvBuff,
		Size,
		MPI_COORD_T,
		RootRank,
		MPI_COMM_WORLD
	);
	return recvBuff;
}

MDPoint getMDPointFromRawBuffer(coord_t* *Place, int Dimensions) {
	MDPoint nowy;
	make_MDPoint(&nowy,Dimensions);
	int i;
	for(i=0;i<Dimensions;i++) {
		nowy.coordinates[i] = *(*Place);
		(*Place)+=1;
	}
	return nowy;
}

// Relocate the points, according to Hilbert Curve
#define SENDING_TAG1  111
void HilbertLibRelocate(
	MDPoint *Data,
	coord_t *HCoordinates,
	coord_t *Boundaries,
	int MyPointsCount,
	int ProcessCount, 
	int Dimensions,
	MDPoint* *NewData,
	int *NewDataCount
) {
	//first send the amounts of particles to send later
	int* sendAmounts = calloc(ProcessCount,sizeof(int));
	int i = 0;
	int wsk = 0;
	int wskBuf = 0;
	int j;
	coord_t* sendBuf = calloc(MyPointsCount*Dimensions,sizeof(coord_t));
	for(i=0;i<ProcessCount;i++) {
		for(j=0;j<Dimensions;j++) {
			sendBuf[wskBuf] = Data[i].coordinates[j];
			wskBuf++;
		}
		while(HCoordinates[wsk] <= Boundaries[i]) {
			wsk++;
			sendAmounts[i]++;
		}
	}
	coord_t* recvAmounts = calloc(ProcessCount,sizeof(int));
	MPI_Alltoall(
		sendAmounts,
		1,
		MPI_INT,
		recvAmounts,
		1,
		MPI_INT,
		MPI_COMM_WORLD
	);
	// All processes now know who will give them what amount of particles
	int pref = 0;
	// iSending
	MPI_Request *requestList, requestNull;
	requestList = calloc(ProcessCount,sizeof(MPI_Request));
	for(i=0;i<ProcessCount;i++) {
		if(sendAmounts[i] == 0)
			continue;
		MPI_Isend(
			sendBuf+pref,
			sendAmounts[i],
			MPI_COORD_T,
			i,
			SENDING_TAG1,
			MPI_COMM_WORLD,
			&requestNull
		);

		pref += sendAmounts[i]* Dimensions;
	}
	// iRecving
	int myNewPointsSize = 0;
	for(i=0;i<ProcessCount;i++) {
		myNewPointsSize += recvAmounts[i];
	}
	*NewDataCount = myNewPointsSize;
	(*NewData) = calloc(myNewPointsSize,sizeof(MDPoint));
	coord_t* recvNewDataBuf = calloc(myNewPointsSize*Dimensions,sizeof(coord_t));
	pref = 0;
	for(i=0;i<ProcessCount;i++) {
		if(recvAmounts[i] == 0)
			continue;
		MPI_Irecv(
			recvNewDataBuf+pref,
			recvAmounts[i],
			MPI_COORD_T,
			i,
			SENDING_TAG1,
			MPI_COMM_WORLD,
			&requestList[i]
		);
			
		pref += recvAmounts[i]* Dimensions;
	}
	// Waiting
	//MPI_Status *statuses = calloc(ProcessCount,sizeof(MPI_Status));
	MPI_Waitall(ProcessCount,requestList,MPI_STATUSES_IGNORE);
	coord_t* Data_ptr = recvNewDataBuf;
	for(i=0;i<ProcessCount;i++) {
		for(j=0;j<recvAmounts[i];j++) {
			(*NewData)[i] = 
				getMDPointFromRawBuffer(&Data_ptr,Dimensions);
		}
	}

}

int main (int argc, char *argv[]) {
	// Initialization
	MPI_Init(&argc, &argv);
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("Witaj, jestem procesem # %d\n",rank);
	#define ROOT 0
	#define DIMENSIONS 3
	#define BITS_PRECISION 5
	
	// Random Input Generation
	srand(time(NULL)+rank);
	int MyPointsCount = rand()%2+1;
	int i,j;
	MDPoint *MyPoints = calloc(MyPointsCount,sizeof(MDPoint));
	for(i=0;i<MyPointsCount;i++) {
		make_MDPoint(&MyPoints[i],DIMENSIONS);
		for(j=0;j<DIMENSIONS;j++) {
			MyPoints[i].coordinates[j] = rand()%(1<<BITS_PRECISION);
		}
	}
	
	//Printing genereated points
	/*for(i=0;i<MyPointsCount;i++) {
		printf("Punkt #%d : ",i);
		for(j=0;j<DIMENSIONS;j++) {
			printf("%u ", MyPoints[i].coordinates[j]);
		}
		printf("\n");

	}*/

	MDPoint *SortedData;
	coord_t *HCoordinates;
	HilbertLibNodeCurveSort(
		MyPoints,
		&SortedData,
		&HCoordinates,
		MyPointsCount,
		BITS_PRECISION,
		DIMENSIONS
	);
	//Removed, because SortedData is using this data now
	//for(i=0;i<MyPointsCount;i++) 
	//	MDPointRemove(&MyPoints[i]);
	free(MyPoints);
	MyPoints = NULL;
	//Printing Sorted Points, and their HCoordinates
	for(i=0;i<MyPointsCount;i++) {
		printf("Punkt #%d : ",i);
		for(j=0;j<DIMENSIONS;j++) {
			printf("%d ", SortedData[i].coordinates[j]);
		}
		printf("  | H : %d",HCoordinates[i]);
		printf("\n");

	}
	
	int * boundaries;
	if(rank == ROOT) {
		boundaries = HilbertLibRootMakeBins(
			ROOT,
			size,
			HCoordinates,
			MyPointsCount,
			BITS_PRECISION,
			DIMENSIONS
		);
		HilbertLibSendBoundariesToAll(boundaries,size,ROOT);
	} else {
		HilbertLibNodeMakeBins(
			HCoordinates,
			MyPointsCount,
			ROOT
		);
		boundaries = HilbertLibRecvBoundariesFromRoot(size,ROOT);
	}
	MDPoint* NewData;
	int NewDataCount;
	HilbertLibRelocate(
		SortedData,
		HCoordinates,
		boundaries,
		MyPointsCount,
		size,
		DIMENSIONS,
		&NewData,
		&NewDataCount
	);

	printf("Żegnaj, jestem procesem #%d\n",rank);
	MPI_Finalize();
	return 0;
}


