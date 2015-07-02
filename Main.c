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
#define calloc(a,b) (a==0 ? NULL : calloc(a,b))
#include<mpi.h>
#include<stdlib.h>
#include<stdio.h>
#include "HilbertLib.h"

#define ROOT 0
#define DIMENSIONS 3
#define BITS_PRECISION 20


coord_t rand_coord_t() {
	return rand()%(1ll<<BITS_PRECISION-1);
}

int main (int argc, char *argv[]) {
	// Initialization
	
	MPI_Init(&argc, &argv);


	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	FILE * File;
	char *arr = calloc(30,sizeof(char));;
	sprintf(arr,"MainOutput%d",rank);
	File = fopen(arr,"w");
	free(arr);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MDPoint * MyPoints;
	int MyPointsCount;
	int i,j;
	// Random Input Generation
	if(1 | (rank != 0)) {
		MyPointsCount = 5;
		srand(time(0) + rank+size);
		MyPoints = calloc(MyPointsCount,sizeof(MDPoint));
		for(i=0;i<MyPointsCount;i++) {
			make_MDPoint(&MyPoints[i],DIMENSIONS);
			for(j=0;j<DIMENSIONS;j++) {
				MyPoints[i].coordinates[j] = 
					rand_coord_t();
			}
			MyPoints[i].own_data_id = 
				MyPoints[i].coordinates[0]^17;
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
	for(i=0;i<MyPointsCount;i++) {
		fprintf(File,"Punkt #%d : ",i);
		for(j=0;j<DIMENSIONS;j++) {
			fprintf(File,"%u ", MyPoints[i].coordinates[j]);
		}
		fprintf(File,"\n");

	}
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
	fprintf(File,"NewDataCount[%d] = %d\n",rank,NewDataCount);
	for(i=0;i<NewDataCount;i++) {

		fprintf(File,"rank:#%d point:#%d @@@   ",rank,i);
		for(j=0;j<DIMENSIONS;j++) {
			fprintf(File,"%d ", NewData[i].coordinates[j]);
		}
		fprintf(File,"\n");
	}
	
	MTNode * MyTreeRoot = HilbertLibPrepareNodeForQueries (
		NewData,
		NewDataCount,
		DIMENSIONS
	);
	int *number_of_queries = calloc(size,sizeof(int));
	MDPoint *QueryPoints;
	int QueryPointsSize;
	MPI_Request req1;
	char * BigBuff = NULL;
	MDPoint*** SelfQueriesResult;
	int* SelfQueriesResultCount;
	int* SelfQueriesRank;
	int SelfQueriesCount;
	

	if(rank == 0) {
		exchangeNumberOfQueries(
			&number_of_queries,
			1
		);
		coord_t *LD = calloc(DIMENSIONS,sizeof(coord_t));
		coord_t *RD = calloc(DIMENSIONS,sizeof(coord_t));
		for(i=0;i<DIMENSIONS;i++) {
			LD[i] = rand_coord_t();
			RD[i] = rand_coord_t();
			if(LD[i] > RD[i]) {
				LD[i] ^= RD[i];
				RD[i] ^= LD[i];
				LD[i] ^= RD[i];
			}
			assert(LD[i] <= RD[i]);
		}
		sendQuery(
			LD,
			RD,
			size,
			DIMENSIONS,
			&QueryPoints,
			&QueryPointsSize,
			rank,
			0,
			&req1,
			&BigBuff,
			1
		);
		MPI_Request_free(&req1);
	} else {
		exchangeNumberOfQueries(
			&number_of_queries,
			size
		);
	}
	answerQueries(
		size,
		DIMENSIONS,
		NewData,
		NewDataCount,
		MyTreeRoot,
		number_of_queries,
		rank,
		&SelfQueriesResult,
		&SelfQueriesResultCount,
		&SelfQueriesRank,
		&SelfQueriesCount
	);
	MDPoint* NewNeighbours;
	int NewNeighboursSize;
	MDPoint* *Results;
	int* ResultsSize;

	if(rank == 0) {
		recvQueries(
			&NewNeighbours,
			&NewNeighboursSize,
			&Results,
			DIMENSIONS,
			size,
			1,
		
		)
	}
	for(i=0;i<NewDataCount;i++) {
		MDPointRemove(&NewData[i]);
	}
	MTDelete(MyTreeRoot);
	free(NewData);
	free(BigBuff);
	fclose(File);
	MPI_Finalize();
	
	return 0;
	/*int *newPointsCount = calloc(size,sizeof(int));
	MPI_AlltoAll(&NewDataCount, 1, MPI_INT, newPointsCount, 1, MPI_INT, MPI_COMM_WORLD);

	//fprintf("%d\n",NewDataCount);
	MPI_File *fh;
	MPI_File_open(MPI_COMM_WORLD, "mainResult.dat", MPI_MODE_RDWR, MPI_INFO_NULL, fh);
	MPI_File_set_size(fh,0);

	int myoffset = 0;
	for(i=0;i<rank;i++)
		myoffset += newPointsCount[i];
	MPI_File_seek(offset,1,MPI_SEEK_SET);
	char *line = calloc(128,sizeof(char));
	
	for(i=0;i<NewDataCount;i++) {
		
		for(j=0;j<DIMENSIONS;j++) {
			sfprintf(line, "%d ", NewData[i].coordinates[j]);
			MPI_File_write(
				fh,

				
		}
		//fprintf("data_id = %d",NewData[i].own_data_id);
		sfprintf(line,"%d", rank);
		sfprintf(line,"\n");
		fprintf("%s",line);
	}*/

	//free(line);*/
	//MPI_File_close(fh);
}
