#include "MDPoint.h"
#include "MyTree.h"
#ifndef HILBERTLIBDEFINED
#define HILBERTLIBDEFINED
void HilbertLibPartition(
	MDPoint *MyPoints,
	int MyPointsCount,
	int RootRank,
	int Dimensions,
	int BitsPrecision,
	int rank,
	int size,
	MDPoint * *NewDataPtr,
	int * NewDataSize
);

MTNode* HilbertLibPrepareNodeForQueries (
	MDPoint *Data,
	int DataSize,
	int Dimensions
);

void exchangeNumberOfQueries (
	int* *RecvSend,
	int NumberOfProcesses
);

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
);

void answerQueries(
	int NumberOfProcesses,
	int Dimensions,
	MDPoint *Data,
	int DataSize,
	MTNode *Root,
	int* RecvCount,
	int MyRank,
	MDPoint*** *SelfQueriesResult,
	int* *SelfQueriesResultCount,
	int* *SelfQueriesRank,
	int* SelfQueriesCount
);

void recvQueries( 
	MDPoint* *NewNeighbours,
	int *NewNeighboursSize,
	MDPoint** *Results,
	int Dimensions,
	int ProcessCount,
	int QueriesCount,
	int selfQueriesCount,
	MDPoint*** SelfQueriesResult,
	int* SelfQueriesResultCount,
	int* SelfQueriesRank,
	int  SelfQueriesCount,
	int  MyRank
);
#endif
