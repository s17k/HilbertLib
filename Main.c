#include<mpi.h>
#include<stdlib.h>
#include<stdio.h>
#include "HilbertLib.h"
int main (int argc, char *argv[]) {
	// Initialization
	MPI_Init(&argc, &argv);
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	#define ROOT 0
	#define DIMENSIONS 3
	#define BITS_PRECISION 30
	
	MDPoint * MyPoints;
	int MyPointsCount;
	int i,j;
	// Random Input Generation
	if(1 | (rank != 0)) {
		MyPointsCount = 10000;
		srand(rank+size);
		MyPoints = calloc(MyPointsCount,sizeof(MDPoint));
		for(i=0;i<MyPointsCount;i++) {
			make_MDPoint(&MyPoints[i],DIMENSIONS);
			for(j=0;j<DIMENSIONS;j++) {
				MyPoints[i].coordinates[j] = rand()%(1ll<<BITS_PRECISION-1);
			}
			MyPoints[i].own_data_id = MyPoints[i].coordinates[0]^17;
		}
	} else {
		MyPointsCount = 11659;
		MyPoints = calloc(MyPointsCount,sizeof(MDPoint));
		for(i=0;i<MyPointsCount;i++) {
			make_MDPoint(&MyPoints[i],DIMENSIONS);
			for(j=0;j<DIMENSIONS;j++) {
				scanf("%d ", &MyPoints[i].coordinates[j]);
			}
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

	int *newPointsCount = calloc(size,sizeof(int));
	MPI_AlltoAll(&NewDataCount, 1, MPI_INT, newPointsCount, 1, MPI_INT, MPI_COMM_WORLD);

	//printf("%d\n",NewDataCount);
	MPI_File *fh;
	MPI_File_open(MPI_COMM_WORLD, "mainResult.dat", MPI_MODE_RDWR, MPI_INFO_NULL, fh);
	MPI_File_set_size(fh,0);

	int myoffset = 0;
	for(i=0;i<rank;i++)
		myoffset += newPointsCount[i];
	MPI_File_seek(offset,1,MPI_SEEK_SET);
	for(i=0;i<NewDataCount;i++) {
		for(j=0;j<DIMENSIONS;j++) {
			//printf("%d ", NewData[i].coordinates[j]);
			MPI_File_write(
				fh,

				
		}
		//printf("data_id = %d",NewData[i].own_data_id);
		printf("%d", rank);
		printf("\n");*/
		MDPointRemove(&NewData[i]);
	}
	free(NewData);
	MPI_File_close(fh);
	MPI_Finalize();
	return 0;
}
