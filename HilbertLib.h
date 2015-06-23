#include "MDPoint.h"
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
#endif
