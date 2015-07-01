#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>
#include "AxesTranspose.h"
#include "HilbertLib.h"
#include <assert.h>
#include "MDPoint.h"
#include "MyTree.h"
#include "Pair.h"

MDPoint* HomePtr;
hilpos_t *HilbertPos;
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
hilpos_t* *HCoordinates, // Hilbert coordinates to return (output)
int Datasize, // |data| - number of particles in this node (input)
int b,// precision (input)
int n // dimensions (input)
) 
{
	if(Datasize == 0) {
		(*SortedData) = NULL;
		(*HCoordinates) = NULL;
		return;
	}
	int i=0;
	coord_t* tmp = calloc(n,sizeof(coord_t));
	coord_t* tmp2 = calloc(n,sizeof(coord_t));
	hilpos_t* first_elem = calloc(Datasize,sizeof(hilpos_t));
	for(i=0;i<Datasize;i++) { // saving i-th Hilbert Coordinates in res[i].coords[0]
		first_elem[i] = 
			GetHCoordinate(X[i].coordinates,tmp2,tmp,b,n);
	}
	free(tmp);
	free(tmp2);
	MDPoint* *ptrs = calloc(Datasize,sizeof(MDPoint*));
	for(i=0;i<Datasize;i++) 
		ptrs[i] = &X[i];
	HomePtr = X;
	HilbertPos = first_elem;
	qsort(ptrs,Datasize,sizeof(MDPoint*),HilbertLibCurveSortComparator);
	*HCoordinates = calloc(Datasize,sizeof(hilpos_t));
	//free(first_elem); // needed for qsort only (but returned)
	(*SortedData) = calloc(Datasize,sizeof(MDPoint));
	for(i=0;i<Datasize;i++) {
		(*SortedData)[i] = *(ptrs[i]);
		(*HCoordinates)[i] = first_elem[ptrs[i]-(&X[0])];
	}
	//free(X);
	free(ptrs);
	free(first_elem);
	
}

// how many particles have hcoordinates <= Right [Node]
int HilbertLibNodeBinSearch(
hilpos_t *HCoordinates, 
int Datasize, 
hilpos_t Right) 
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
int HilbertLibNodeHowMany(hilpos_t *HCoordinates, int Datasize, hilpos_t Left, hilpos_t Right) { 
	if(Datasize == 0)
		return 0;
	return 
		HilbertLibNodeBinSearch(HCoordinates,Datasize,Right) - 
		HilbertLibNodeBinSearch(HCoordinates,Datasize,Left);
}



void HilbertLibNodeGetMINMAX(hilpos_t* HCoordinates, int Datasize, hilpos_t *MIN, hilpos_t *MAX) {
	int i;
	(*MIN) = HILPOS_T_MAX_VALUE;
	(*MAX) = 0;
	if(Datasize == 0) return;
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
	hilpos_t* MIN, 
	hilpos_t* MAX,
	hilpos_t *HCoordinates
) {
	int *recvbuf = calloc(ProcessCount,sizeof(int));
	int sendbuf = PointCount;
	MPI_Gather(&sendbuf,1,MPI_INT,recvbuf,1,MPI_INT,RootRank,MPI_COMM_WORLD);
	int suma = 0;
	int i;
	for(i=0;i<ProcessCount;i++) 
		suma += recvbuf[i];
	
	free(recvbuf);

	HilbertLibNodeGetMINMAX(HCoordinates,PointCount,MIN,MAX);
	hilpos_t* sendBuf = calloc(2,sizeof(hilpos_t));
	sendBuf[0] = *MIN;
	sendBuf[1] = *MAX;
	hilpos_t *recvBuf = calloc(2*ProcessCount,sizeof(hilpos_t));
	MPI_Gather(sendBuf,2,MPI_HILPOS_T,recvBuf,2,MPI_HILPOS_T,RootRank,MPI_COMM_WORLD);
	
	//gather Min and Max from other processes
	hilpos_t MinRes = HILPOS_T_MAX_VALUE, MaxRes = 0;
	for(i=0;i<ProcessCount*2;i++) {
		if(recvBuf[i] < MinRes)
			MinRes = recvBuf[i];
		if(recvBuf[i+1] > MaxRes)
			MaxRes = recvBuf[i+1];
	}
	*MIN = MinRes;
	*MAX = MaxRes;
        //printf("MIN %f MAX %f\n",*MIN,*MAX);
	free(sendBuf);
	free(recvBuf);
        return suma;
}

// Calculates the next boundary [Root]
hilpos_t HilbertLibCalculateNextBoundary( 
hilpos_t a, hilpos_t b, 
hilpos_t *MyHCoordinates,
int MyPointsCount,
int particlesRate,
int nodesCount,
int RootRank,
int* how_many_used // how many used is not filled
) {
	hilpos_t bsleft = a, bsright = b, bsmiddle;
	int* hm = calloc(nodesCount,sizeof(int));
	int i = 0;
	int suma;
	hilpos_t * sendBuff = calloc(2,sizeof(hilpos_t));
	while(bsright/bsleft > 1 + HILPOS_EPS) { // CHANGE FOR UNSIGNED INTS
                
		bsmiddle = (bsleft + bsright + HILPOS_BS_1)/HILPOS_BS_2;			
                //printf("bsleft : %f, bsright %f, bsmiddle %f\n",bsleft,bsright,bsmiddle);
		sendBuff[0] = a;
		sendBuff[1] = bsmiddle;
		//printf("Wysyłam parę (%d,%d)\n",a,bsmiddle);
		MPI_Bcast(sendBuff,2,MPI_HILPOS_T,RootRank,MPI_COMM_WORLD);
		int singlebuff = HilbertLibNodeHowMany(
			MyHCoordinates,
			MyPointsCount,
			sendBuff[0],
			sendBuff[1]
		);

		MPI_Reduce(&singlebuff,&suma,1,MPI_INT,MPI_SUM,RootRank,MPI_COMM_WORLD);
		/*MPI_Gather(&singlebuff,1,MPI_INT,hm,1,MPI_INT,RootRank,MPI_COMM_WORLD);
		suma = 0;
		for(i=0;i<nodesCount;i++) {
			suma += hm[i];		
		}*/
		//printf("na tym przedziale jest %d\n",suma);
		if(suma > particlesRate) {
			bsright = bsmiddle-HILPOS_EPS;
		} else {
			bsleft = bsmiddle;
		}
	}	
	// Getting the information about how many points are in (a,bsleft>
	sendBuff[0] = a;
	sendBuff[1] = bsleft;
	MPI_Bcast(sendBuff,2,MPI_HILPOS_T,RootRank,MPI_COMM_WORLD);

	int singlebuff = HilbertLibNodeHowMany(
		MyHCoordinates,
		MyPointsCount,
		sendBuff[0],
		sendBuff[1]
	);

	MPI_Reduce(&singlebuff,&suma,1,MPI_INT,MPI_SUM,RootRank,MPI_COMM_WORLD);
	/*MPI_Gather(&singlebuff,1,MPI_INT,hm,1,MPI_INT,RootRank,MPI_COMM_WORLD);
	suma = 0;
	for(i=0;i<nodesCount;i++) {
		suma += hm[i];
	}*/
	(*how_many_used) = suma;
	free(sendBuff);
	free(hm);
	return bsleft;
}

// Make Bins [Root]
hilpos_t* HilbertLibRootMakeBins(
int RootRank, // rank of root
size_t NodesCount, // number of nodes
hilpos_t *MyHCoordinates, // Hilbert Coordinates of root points
size_t MyPointsCount, // Number of points assigned to root
int b, // |precision bits|
int dimensions // |dimensions|
)
{
	hilpos_t MaxCoord, MinCoord;
	int allPointsCount = 
		HilbertLibGetNOfParticles(
			NodesCount,
			MyPointsCount,
			RootRank,
			&MinCoord,
			&MaxCoord,
			MyHCoordinates
		);
	//printf("AllPointsCount : %d\n",allPointsCount);
	int particlesRate = allPointsCount/NodesCount;
	int i = 0;
	// in binsBoundaries[i] there will be last Hilbert Coordinate for process i
	hilpos_t* binsBoundaries = calloc(NodesCount,sizeof(hilpos_t));
	hilpos_t lastUsed = MinCoord-HILPOS_EPS;
	for(i=0;i<NodesCount;i++) {
		int how_many_used = 0;
		//printf("Licze Boundary dla %d ...   \n ", i);
		binsBoundaries[i] =
			HilbertLibCalculateNextBoundary(
				lastUsed,MaxCoord,
				MyHCoordinates,
				MyPointsCount,
				particlesRate,
				NodesCount,
				RootRank,
				&how_many_used
			);
		//printf("I wyszlo mi %d ... \n ", binsBoundaries[i]);
		allPointsCount -= how_many_used;
		lastUsed = binsBoundaries[i];
		if(i!=(NodesCount-1))
			particlesRate = allPointsCount/(NodesCount-i-1); 
	}
	// last scatter to indicate the end of queries
	hilpos_t* sendbuf = calloc(2,sizeof(hilpos_t));
	sendbuf[0] = 1;
        sendbuf[1] = 0;
	MPI_Bcast(
		sendbuf,
		2,
		MPI_HILPOS_T,
		RootRank,
		MPI_COMM_WORLD
	);
	free(sendbuf);
	return binsBoundaries;
}

//Sends # of particles to Root [Node]
void HilbertLibSendNOfParticles(int DataSize, int RootRank, hilpos_t* HCoordinates) {
	int sendBuff = DataSize;
	MPI_Gather(&sendBuff,1,MPI_INT,NULL,0,NULL,RootRank,MPI_COMM_WORLD);
	hilpos_t* sendBuff2 = calloc(2,sizeof(hilpos_t));
	hilpos_t MIN,MAX;
	HilbertLibNodeGetMINMAX(HCoordinates, DataSize, &MIN, &MAX);
	sendBuff2[0] = MIN;
	sendBuff2[1] = MAX;
	MPI_Gather(sendBuff2,2,MPI_HILPOS_T,NULL,0,NULL,RootRank,MPI_COMM_WORLD);
	free(sendBuff2);
}

//divide points into bins[Node]
void HilbertLibNodeMakeBins(hilpos_t *MyHCoordinates, size_t MyParticlesCount, int RootRank) {
	//Send number of particles i have
	HilbertLibSendNOfParticles(MyParticlesCount, RootRank, MyHCoordinates);	
	hilpos_t *recvBuff = calloc(2,sizeof(hilpos_t));
	//printf("Zaczynam odpowiadac\n");
	while(true) {
		MPI_Bcast(
			recvBuff,
			2,
			MPI_HILPOS_T,
			RootRank,
			MPI_COMM_WORLD
		);
		if(recvBuff[1] < recvBuff[0]) {
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
		MPI_Reduce(
			&sendbuf,
			NULL,
			1,
			MPI_INT,
			MPI_SUM,
			RootRank,
			MPI_COMM_WORLD
		);
		//MPI_Gather(&sendbuf,1,MPI_INT,NULL,0,0,RootRank,MPI_COMM_WORLD);
	}
	free(recvBuff);
}

// [Root]
void HilbertLibSendBoundariesToAll(hilpos_t* Boundaries, int Size, int RootRank) {
	//printf("Rozsylam\n");
	MPI_Bcast(
		Boundaries,
		Size,
		MPI_HILPOS_T,
		RootRank,
		MPI_COMM_WORLD
	);
	//printf("Rozeslalem\n");
}


//[Node]
hilpos_t * HilbertLibRecvBoundariesFromRoot(int Size, int RootRank) {
	//printf("Wchodze\n");
	hilpos_t* recvBuff = calloc(Size,sizeof(hilpos_t));
	MPI_Bcast(
		recvBuff,
		Size,
		MPI_HILPOS_T,
		RootRank,
		MPI_COMM_WORLD
	);
	//printf("Wychodze\n");
	return recvBuff;
}

MDPoint getMDPointFromRawBuffer(coord_t* Place, int Dimensions) { // can be easily fasten
	MDPoint nowy;
	make_MDPoint(&nowy,Dimensions);
	int i;
	for(i=0;i<Dimensions;i++) {
		nowy.coordinates[i] = *(Place+i);
	}
	return nowy;
}

// Relocate the points, according to Hilbert Curve
#define SENDING_TAG1  111
#define SENDING_TAG2  222
void HilbertLibRelocate(
	MDPoint *Data,
	hilpos_t *HCoordinates,
	hilpos_t *Boundaries,
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
	tag_t* tagSendBuf = calloc(MyPointsCount,sizeof(tag_t));
	coord_t* sendBuf = calloc(MyPointsCount*Dimensions,sizeof(coord_t));
	for(i=0;i<MyPointsCount;i++) {
		for(j=0;j<Dimensions;j++) {
			sendBuf[wskBuf] = Data[i].coordinates[j];
			wskBuf++;
		}
		tagSendBuf[i] = Data[i].own_data_id;
                while(Boundaries[wsk]/HCoordinates[i] < 1 - HILPOS_EPS*2) { // Change with Unsigned
                    wsk++;
                    assert(wsk < ProcessCount);
                }
                sendAmounts[wsk]++;
                    
	}
	int* recvAmounts = calloc(ProcessCount,sizeof(int));
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
	
	for(i=0;i<ProcessCount;i++) {
		if(sendAmounts[i] == 0)
			continue;
		MPI_Isend(
			sendBuf+pref,
			sendAmounts[i]*Dimensions,
			MPI_COORD_T,
			i,
			SENDING_TAG1,
			MPI_COMM_WORLD,
			&requestNull
		);
		MPI_Isend(
			tagSendBuf + (pref/Dimensions),
			sendAmounts[i],
			MPI_TAG_T,
			i,
			SENDING_TAG2,
			MPI_COMM_WORLD,
			&requestNull
		);

		pref += sendAmounts[i]* Dimensions;
	}
	free(sendAmounts);
	// iRecving
	int myNewPointsSize = 0;
        int zeros = 0;
	for(i=0;i<ProcessCount;i++) {
		myNewPointsSize += recvAmounts[i];
                if(recvAmounts[i] == 0)
                    zeros += 1;
	}
	int actual_size = (ProcessCount - zeros)*2;
        if(actual_size > 0)
            requestList = calloc(actual_size,sizeof(MPI_Request));
        else   
            requestList = NULL;
        
	*NewDataCount = myNewPointsSize;
	(*NewData) = calloc(myNewPointsSize,sizeof(MDPoint));
	coord_t* recvNewDataBuf = calloc(myNewPointsSize*Dimensions,sizeof(coord_t));
	tag_t* recvTagTBuf = calloc(myNewPointsSize,sizeof(tag_t));
	pref = 0;
        unsigned int reqptr = 0;
	for(i=0;i<ProcessCount;i++) {
		if(recvAmounts[i] == 0)
			continue;
		MPI_Irecv(
			recvNewDataBuf+pref,
			recvAmounts[i]*Dimensions,
			MPI_COORD_T,
			i,
			SENDING_TAG1,
			MPI_COMM_WORLD,
			&requestList[reqptr*2]
		);
		MPI_Irecv(
			recvTagTBuf + (pref/Dimensions),
			recvAmounts[i],
			MPI_TAG_T,
			i,
			SENDING_TAG2,
			MPI_COMM_WORLD,
			&requestList[reqptr*2+1]
		);
                reqptr++;
			
		pref += recvAmounts[i]* Dimensions;
	}
	// Waiting
        if(ProcessCount-zeros != 0)
            MPI_Waitall(actual_size,requestList,MPI_STATUSES_IGNORE);
        int li =0;
	for(i=0;i<ProcessCount;i++) {
		for(j=0;j<recvAmounts[i];j++) {
			make_MDPoint(&((*NewData)[li]),Dimensions);
				//getMDPointFromRawBuffer(recvNewDataBuf+li*Dimensions,Dimensions);
                        memcpy(((*NewData)[li]).coordinates,recvNewDataBuf+(li*Dimensions),sizeof(coord_t)*Dimensions);
			(*NewData)[li].own_data_id = recvTagTBuf[li];
                        li+=1;          
		}
	}
	free(sendBuf);
	free(tagSendBuf);
	free(recvAmounts);
	free(requestList);
	free(recvNewDataBuf);
	free(recvTagTBuf);
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
	hilpos_t *HCoordinates;
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
	/*for(i=0;i<MyPointsCount;i++) {
		//printf("Punkt #%d : ",i);
		for(j=0;j<Dimensions;j++) {
			printf("%d ", SortedData[i].coordinates[j]);
		}
		//printf("  | H : %d",HCoordinates[i]);
		printf("\n");

	}*/
	
	hilpos_t * boundaries = NULL;
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
	}
	/*if(rank == 0)
	for(i=0;i<size;i++) {
		printf("boundaries[%d] = %e\n",i,boundaries[i]);
	}*/
	
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
	(*NewDataSize) = NewDataCount;
	free(SortedData);
	free(HCoordinates);
	free(boundaries);
}

MTNode* HilbertLibPrepareNodeForQueries (
	MDPoint *Data,
	int DataSize,
	int Dimensions
) {
	MTNode* root = calloc(1,sizeof(MTNode));
	makeMTNode(root,0,0);
	MDPoint* *temp = calloc(DataSize,sizeof(MDPoint*));
	int i;
	for(i=0;i<DataSize;i++) {
		temp[i] = &Data[i];
	}
	MTmake(
		root,
		temp,
		DataSize,
		Dimensions,
		0
	);
	free(temp);
	return root;
}

void exchangeNumberOfQueries (
	int* *RecvSend,
	int NumberOfProcesses
) {
	*RecvSend = calloc(NumberOfProcesses,sizeof(int));
	MPI_Allgather(
		&NumberOfProcesses,
		1,
		MPI_INT,
		*(RecvSend),
		1,
		MPI_INT,
		MPI_COMM_WORLD
	);
}

char* sendQuery(
	coord_t *LD,
	coord_t *RD,
	int NumberOfProcesses,
	int Dimensions,
	MDPoint* *Res, 
	int* ResSize,
	int MyRank,
	int Counter, // unique (for query) number in <0,n-1>, where n is # of querie
	MPI_Request* Request,
	char* *BigBuff, // should be NULL on first query and freed after
	int numberOfMyQueries
) {
	int number_of_bytes = sizeof(coord_t) * 2 * Dimensions + sizeof(int);
	if(*BigBuff == NULL) {
		*BigBuff = malloc(number_of_bytes * numberOfMyQueries * sizeof(coord_t));
	}
	char* buff = (*BigBuff) + Counter * number_of_bytes * sizeof(coord_t);

	memcpy(buff,LD,sizeof(coord_t)*Dimensions);
	memcpy(buff+sizeof(coord_t)*Dimensions,RD,sizeof(coord_t)*Dimensions);
	memcpy(buff+sizeof(coord_t)*Dimensions*2,&Counter,sizeof(int));
	
	MPI_Ibcast(
		buff,
		number_of_bytes,
		MPI_BYTE,
		MyRank,
		MPI_COMM_WORLD,
		Request
	);

}

void answerQueries(
	int NumberOfProcesses,
	int Dimensions,
	MDPoint *Data,
	int DataSize,
	MTNode *Root,
	int* RecvCount,
	int MyRank,
	MDPoint*** *SelfQueriesResult,
	int* *SelfQueries
) {
	int i,j,k,l;
	int number_of_bytes = sizeof(coord_t) * 2 * Dimensions + sizeof(int);
	char* buff = malloc(number_of_bytes);
	int numberOfQueries = 0;
	for(i=0;i<NumberOfProcesses;i++) 
			numberOfQueries += RecvCount[i];
	char* *results = calloc(numberOfQueries,sizeof(char*));
	int* boolSend = calloc(DataSize,sizeof(int));
	for(i=0;i<DataSize;i++)
		boolSend[i] = -1;
	char* *realBuffs = calloc(NumberOfProcesses,sizeof(char*));
	int* *infoBuffs = calloc(NumberOfProcesses,sizeof(int*));
	for(i=0;i<NumberOfProcesses;i++) {
		if(RecvCount[i] == 0)
				continue;
		int allResults = 0; // counted twice 
		int resultsCount = 0; // counted once
		MDPoint** *Res = calloc(RecvCount[i],sizeof(MDPoint**));
		int *resSize = calloc(RecvCount[i],sizeof(int));
		int *queryRank = calloc(RecvCount[i],sizeof(int));
		for(j=0;j<RecvCount[i];j++) {
			MPI_Request req; 
			MPI_Ibcast(
				buff,
				1,
				MPI_INT,
				i,
				MPI_COMM_WORLD,
				&req
			);
			MPI_Wait(&req,MPI_STATUS_IGNORE);
			coord_t *LD = (coord_t*)buff;
			coord_t *RD = (coord_t*)(buff + Dimensions*sizeof(coord_t));
			queryRank[j] = *((int*)(buff + sizeof(coord_t) * Dimensions * 2));
			MTQuery(
				Root,
				LD,
				RD,
				&Res[j],
				&resSize[j],
				Dimensions
			);
			allResults += resSize[j];
			for(k=0;k<resSize[j];k++) {
				if(boolSend[Res[j][k]-Data] != i) {
					boolSend[Res[j][k]-Data] = i;
					resultsCount++;
				}
			}
		}
		Pair *toSort = calloc(allResults,sizeof(Pair));
		int* whoseIsThatPoint = calloc(allResults,sizeof(int)); // can be omitted
		int wsk = 0;
		for(j=0;j<RecvCount[i];j++) {
			for(k=0;k<resSize[j];k++) {
				make_Pair(&toSort[wsk], Res[j][k]-Data, queryRank[j]);
				wsk+=1;
			}
		}
		qsort(toSort, allResults, sizeof(Pair), PairComparator);
		wsk = 0;
		for(j=0;j<allResults;j++) {
			whoseIsThatPoint[wsk] = toSort[j].b;
			if(j == 0 || toSort[j].a != toSort[j-1].a)
				whoseIsThatPoint[wsk] *= (-1);
			wsk+=1;
		}
		free(Res);
		int *info_buff = calloc(2,sizeof(int));
		info_buff[0] = allResults;
		info_buff[1] = resultsCount;
		MPI_Isend(
			info_buff,
			2,
			MPI_INT,
			i,
			765,
			MPI_COMM_WORLD,
			NULL
		);
		infoBuffs[i] = info_buff;
		int real_buff_size = 
			info_buff[0]*sizeof(int) + info_buff[1]*MDPOINT_FULL_SIZE(Dimensions);
		char* real_buff = malloc(real_buff_size);
		char* real_buff_ptr = real_buff;
		realBuffs[i] = real_buff;
		for(j=0;j<allResults;j++) {
			if(j == allResults-1 || toSort[j].a != toSort[j+1].a) {
				int offset = 
					MDPointPack(
						real_buff_ptr,
						&Data[toSort[j].a],
						Dimensions
					);
				real_buff_ptr += offset;

			}
		}
		memcpy(real_buff_ptr,whoseIsThatPoint,sizeof(int)*allResults); 
		MPI_Isend(
			real_buff,
			real_buff_size,
			MPI_BYTE,
			i,
			766,
			MPI_COMM_WORLD,
			NULL
		);
		free(whoseIsThatPoint);		
		free(queryRank);
		for(j=0;j<RecvCount[i];j++) {
			free(Res[j]);
		}
		free(resSize);
		if(i == MyRank) {
			
		}
	}
	for(i=0;i<NumberOfProcesses;i++) {
		if(RecvCount[i] == 0)
			continue;
		free(realBuffs[i]);
		free(infoBuffs[i]);
		free(results[i]);
	}

	free(buff);
}

int abs(int x) {
	if(x < 0)
		return -x;
	else
		return x;
}

void recvQueries( // Add MPI_TEST_SOME to fasten
	MDPoint* *NewNeighbours,
	int *NewNeighboursSize,
	MDPoint** *Results,
	int* *ResultSize,
	int Dimensions,
	int ProcessCount,
	int QueriesCount
) {
	int i,j;
	PtrVector *VecResults = calloc(QueriesCount,sizeof(PtrVector));
	for(i=0;i<QueriesCount;i++) {
		makePtrVector(&VecResults[i]);
	}
	int* *cntBuffers = calloc(ProcessCount,sizeof(int*));
	MPI_Request req;
	for(i=0;i<ProcessCount;i++) {
		cntBuffers[i] = calloc(2,sizeof(int));
		int* cntBuf = cntBuffers[i];
		MPI_Irecv(
			cntBuf,
			2,
			MPI_INT,
			i,
			765,
			MPI_COMM_WORLD,
			&req
		);
		MPI_Wait(&req,MPI_STATUS_IGNORE);
		*NewNeighboursSize += cntBuf[1];
	}
	(*NewNeighbours) = calloc(*NewNeighboursSize,sizeof(MDPoint));
	int newNeighboursPtr = 0;
	int lastNewNeighboursPtr = 0;
	for(i=0;i<ProcessCount;i++) {
		int* cntBuf = cntBuffers[i];
		int big_buff_size = 
			cntBuf[0] * sizeof(int) + 
			cntBuf[1] * MDPOINT_FULL_SIZE(Dimensions);
		char* big_buff = malloc(big_buff_size);
		MPI_Irecv(
			big_buff,
			big_buff_size,
			MPI_BYTE,
			i,
			766,
			MPI_COMM_WORLD,
			&req
		);
		MPI_Wait(&req,MPI_STATUS_IGNORE);
		char * big_buff_ptr = big_buff;
		for(i=0;i<cntBuf[1];i++) {
			int offset = 
				MDPointUnpack(
					big_buff_ptr,
					&(((*NewNeighbours)[newNeighboursPtr])),
					Dimensions
				);
			newNeighboursPtr += 1;
			big_buff_ptr += offset;
		}
		
		lastNewNeighboursPtr = newNeighboursPtr;
		int pointPtr = -1;
		for(j=0;j<cntBuf[0];j++) {
			int actVal = *((int*)big_buff_ptr);
			if(actVal < 0)
				pointPtr += 1;
			actVal = abs(actVal);
			PtrVectorPB(&VecResults[actVal],&((*NewNeighbours)[lastNewNeighboursPtr + pointPtr]));
			big_buff_ptr += sizeof(int);
		}

		free(big_buff);
		free(cntBuf);
		// next points from input 
	}
	free(cntBuffers);
	free(VecResults);
	for(i=0;i<QueriesCount;i++) {
		(*Results)[i] = calloc(VecResults[i].size,sizeof(MDPoint*));
		memcpy((*Results)[i], VecResults[i].arr, VecResults[i].size * sizeof(MDPoint*));
		PtrVectorDeallocate(&VecResults[i]);
	}
}






