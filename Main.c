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
	#define ROOT 0
	#define DIMENSIONS 3
	#define BITS_PRECISION 20
	
	MDPoint * MyPoints;
	int MyPointsCount;
	int i,j;
	// Random Input Generation
	if(1 | (rank != 0)) {
		MyPointsCount = 30;
		srand(time(0) + rank+size);
		MyPoints = calloc(MyPointsCount,sizeof(MDPoint));
		for(i=0;i<MyPointsCount;i++) {
			make_MDPoint(&MyPoints[i],DIMENSIONS);
			for(j=0;j<DIMENSIONS;j++) {
				MyPoints[i].coordinates[j] = 
					rand()%(1ll<<BITS_PRECISION-1);
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
		MDPointRemove(&NewData[i]);	
	}
	free(NewData);
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
