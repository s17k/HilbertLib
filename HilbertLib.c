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
		//printf("Znajduje H dla %d,",i);
		AxestoTranspose(res[i].coordinates,b,n);
		//printf(" Wyszlo mi %d\n",res[i].coordinates[0]);
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
	free(first_elem);
	
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



void HilbertLibNodeGetMINMAX(coord_t* HCoordinates, int Datasize, coord_t *MIN, coord_t *MAX) {
	int i;
	(*MIN) = 0;
	(*MAX) = COORD_T_MAX_VALUE;
	for(i = 0; i < Datasize;i++) {
		if(HCoordinates[i] < *MIN)
			*MIN = HCoordinates[i];
		if(HCoordinates[i] > *MAX)
			*MAX = HCoordinates[i];
	}
}

// divide points into bins[Root] 
int HilbertLibGetNOfParticles(
	int ProcessCount, 
	int PointCount, 
	int RootRank, 
	coord_t* MIN, 
	coord_t* MAX,
	coord_t *HCoordinates
) {
	int *recvbuf = calloc(ProcessCount,sizeof(int));
	int sendbuf = PointCount;
	MPI_Gather(&sendbuf,1,MPI_INT,recvbuf,1,MPI_INT,RootRank,MPI_COMM_WORLD);
	int suma = 0;
	int i;
	for(i=0;i<ProcessCount;i++) 
		suma += recvbuf[i];
	return suma;
	free(recvbuf);

	HilbertLibNodeGetMINMAX(HCoordinates,PointCount,MIN,MAX);
	coord_t sendBuf = calloc(2,sizeof(coord_t));
	sendBuf[0] = MIN;
	sendBuf[1] = MAX;
	coord_t *recvBuf = calloc(2*ProcessCount,sizeof(coord_t));
	MPI_Gather(sendBuf,2,MPI_COORD_T,recvBuf,2,MPI_COORD_T,RootRank,MPI_COMM_WORLD);
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
	coord_t * sendBuff = calloc(2,sizeof(coord_t));
	while(bsleft < bsright) {
		//printf("bsleft : %d, bsright %d\n",bsleft,bsright);
		bsmiddle = (bsleft + bsright + 1)/2;			
		sendBuff[0] = a;
		sendBuff[1] = bsmiddle;
		//printf("Wysyłam parę (%d,%d)\n",a,bsmiddle);
		MPI_Bcast(sendBuff,2,MPI_COORD_T,RootRank,MPI_COMM_WORLD);
		int singlebuff = HilbertLibNodeHowMany(
			MyHCoordinates,
			MyPointsCount,
			sendBuff[0],
			sendBuff[1]
		);

		MPI_Gather(&singlebuff,1,MPI_INT,hm,1,MPI_INT,RootRank,MPI_COMM_WORLD);
		suma = 0;
		for(i=0;i<nodesCount;i++) {
			suma += hm[i];		
		}
		//printf("na tym przedziale jest %d\n",suma);
		if(suma > particlesRate) {
			bsright = bsmiddle-1;
		} else {
			bsleft = bsmiddle;
		}
	}	
	// Getting the information about how many points are in (a,bsleft>
	sendBuff[0] = a;
	sendBuff[1] = bsleft;
	MPI_Bcast(sendBuff,2,MPI_COORD_T,RootRank,MPI_COMM_WORLD);

	int singlebuff = HilbertLibNodeHowMany(
		MyHCoordinates,
		MyPointsCount,
		sendBuff[0],
		sendBuff[1]
	);
	MPI_Gather(&singlebuff,1,MPI_INT,hm,1,MPI_INT,RootRank,MPI_COMM_WORLD);
	suma = 0;
	for(i=0;i<nodesCount;i++) {
		suma += hm[i];
	}
	(*how_many_used) = suma;
	free(sendBuff);
	free(hm);
	return bsleft;
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
	coord_t MaxCoord, MinCoord;
	int allPointsCount = 
		HilbertLibGetNOfParticles(
			NodesCount,
			MyPointsCount,
			RootRank,
			&MaxCoord,
			&MinCoord,
			MyHCoordinates
		);
	printf("AllPointsCount : %d\n",allPointsCount);
	int particlesRate = allPointsCount/NodesCount;
	int i = 0;
	// in binsBoundaries[i] there will be last Hilbert Coordinate for process i
	coord_t* binsBoundaries = calloc(NodesCount,sizeof(coord_t));
	coord_t lastUsed = 0;
	for(i=0;i<NodesCount;i++) {
		int how_many_used = 0;
		printf("Licze Boundary dla %d ...   \n ", i);
		binsBoundaries[i] =
			HilbertLibCalculateNextBoundary(
				lastUsed,(1<<b)-1,
				MyHCoordinates,
				MyPointsCount,
				particlesRate,
				NodesCount,
				RootRank,
				&how_many_used
			);
		printf("I wyszlo mi %d ... \n ", binsBoundaries[i]);
		allPointsCount -= how_many_used;
		lastUsed = binsBoundaries[i];
		if(i!=(NodesCount-1))
			particlesRate = allPointsCount/(NodesCount-i-1); 
	}
	// last scatter to indicate the end of queries
	coord_t* sendbuf = calloc(2,sizeof(int));
	sendbuf[0] = sendbuf[1] = -1;
	MPI_Bcast(
		sendbuf,
		2,
		MPI_COORD_T,
		RootRank,
		MPI_COMM_WORLD
	);
	free(sendbuf);
	return binsBoundaries;
}

//Sends # of particles to Root [Node]
void HilbertLibSendNOfParticles(int what, int RootRank, coord_t* HCoordinates) {
	int sendBuff = what;
	MPI_Gather(&sendBuff,1,MPI_INT,NULL,0,NULL,RootRank,MPI_COMM_WORLD);
	coord_t* sendBuff2 = calloc(2,sizeof(coord_t));
	coord_t MIN,MAX;

	sendBuff[0] = 
	MPI_Gather(sendBuff2,2,MPI_COORD_T,NULL,0,NULL,RootRank,MPI_COMM_WORLD);
	free(sendBuff2);
}

//divide points into bins[Node]
void HilbertLibNodeMakeBins(coord_t *MyHCoordinates, size_t MyParticlesCount, int RootRank) {
	//Send number of particles i have
	HilbertLibSendNOfParticles(MyParticlesCount, RootRank);	
	coord_t *recvBuff = calloc(2,sizeof(coord_t));
	//printf("Zaczynam odpowiadac\n");
	while(true) {
		MPI_Bcast(
			recvBuff,
			2,
			MPI_COORD_T,
			RootRank,
			MPI_COMM_WORLD
		);
		if(recvBuff[0] < 0) {
			break;
		}
		//printf("Dostałem %d %d\n",recvBuff[0],recvBuff[1]);
		int sendbuf = HilbertLibNodeHowMany(
			MyHCoordinates,
			MyParticlesCount,
			recvBuff[0],
			recvBuff[1]
		);
		//printf("wiec odpowiadam %d\n",sendbuf);
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
	free(recvBuff);
}

// [Root]
void HilbertLibSendBoundariesToAll(coord_t* Boundaries, int Size, int RootRank) {
	printf("Rozsylam\n");
	MPI_Bcast(
		Boundaries,
		Size,
		MPI_COORD_T,
		RootRank,
		MPI_COMM_WORLD
	);
	printf("Rozeslalem\n");
}


//[Node]
coord_t * HilbertLibRecvBoundariesFromRoot(int Size, int RootRank) {
	printf("Wchodze\n");
	coord_t* recvBuff = calloc(Size,sizeof(coord_t));
	MPI_Bcast(
		recvBuff,
		Size,
		MPI_COORD_T,
		RootRank,
		MPI_COMM_WORLD
	);
	printf("Wychodze\n");
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
	free(sendAmounts);
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
	free(sendBuf);
	free(recvAmounts);
	free(requestList);
	free(recvNewDataBuf);


}

void HilbertLibPartition(
	MDPoint *MyPoints, 
	int MyPointsCount, 
	int RootRank, 
	int Dimensions, 
	int BitsPrecision,
	int rank,
	int size,
	MDPoint* *NewDataPtr,
	int*  NewDataSize
) {
	MDPoint *SortedData;
	coord_t *HCoordinates;
	HilbertLibNodeCurveSort(
		MyPoints,
		&SortedData,
		&HCoordinates,
		MyPointsCount,
		BitsPrecision,
		Dimensions
	);
	free(MyPoints);
	MyPoints = NULL;
	//Printing Sorted Points, and their HCoordinates
	int i,j;
	for(i=0;i<MyPointsCount;i++) {
		printf("Punkt #%d : ",i);
		for(j=0;j<Dimensions;j++) {
			printf("%d ", SortedData[i].coordinates[j]);
		}
		printf("  | H : %d",HCoordinates[i]);
		printf("\n");

	}
	
	coord_t * boundaries = NULL;
	if(rank == RootRank) {
		boundaries = HilbertLibRootMakeBins(
			RootRank,
			size,
			HCoordinates,
			MyPointsCount,
			BitsPrecision,
			Dimensions
		);
		HilbertLibSendBoundariesToAll(boundaries,size,RootRank);
	} else {
		HilbertLibNodeMakeBins(
			HCoordinates,
			MyPointsCount,
			RootRank
		);
		boundaries = HilbertLibRecvBoundariesFromRoot(size,RootRank);
	}/*
	for(i=0;i<size;i++) {
		printf("boundaries[%d] = %d\n",i,boundaries[i]);
	}
	MDPoint* NewData;
	int NewDataCount;
	HilbertLibRelocate(
		SortedData,
		HCoordinates,
		boundaries,
		MyPointsCount,
		size,
		Dimensions,
		&NewData,
		&NewDataCount
	);
	for(i=0;i<MyPointsCount;i++) {
		MDPointRemove(&SortedData[i]);
	}
	(*NewDataPtr) = NewData;
	(*NewDataSize) = NewDataCount;*/
	free(SortedData);
	free(HCoordinates);
	free(boundaries);
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
	srand(rank+size);
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
	MDPoint *NewData = NULL;
	int NewDataCount = 0;

	HilbertLibPartition( // MyPoints is freed
		MyPoints,
		MyPointsCount,
		ROOT,
		DIMENSIONS,
		BITS_PRECISION,
		rank,
		size,
		&NewData,
		&NewDataCount
	);

	for(i=0;i<NewDataCount;i++) {
		MDPointRemove(&NewData[i]);
	}
	free(NewData);
	printf("Żegnaj, jestem procesem #%d\n",rank);
	MPI_Finalize();
	return 0;
}


