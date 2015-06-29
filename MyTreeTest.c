#include "MyTree.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//#include <mpi.h> //temporarly
#define DIMENSIONS 3
#define DATASIZE 6


int main (int argc, char* argv[]) {
	//MPI_Init(&argc,&argv);
	srand(time(NULL));
	MTNode Root;
	makeMTNode(&Root,0,0);
	MDPoint *Data = calloc(DATASIZE,sizeof(MDPoint));
	MDPoint* *DataPtr = calloc(DATASIZE,sizeof(MDPoint*));
	int i=0;
	for(i=0;i<DATASIZE;i++) {
		make_MDPoint(&Data[i],DIMENSIONS);
		int j = 0;
		for(j=0;j<DIMENSIONS;j++) {
			Data[i].coordinates[j] = rand()%8;
			printf("%d ", Data[i].coordinates[j]);
		}
		DataPtr[i] = &Data[i];
		printf("\n");
	}
	MTmake(&Root,DataPtr,DATASIZE,DIMENSIONS,0);
	coord_t *LD = calloc(DIMENSIONS,sizeof(coord_t));
	coord_t *RD = calloc(DIMENSIONS,sizeof(coord_t));
	printf("OKREŚL OBSZAR POSZUKIWAN (l1,r1), (l2,r2) ..., (ln,rn)\n");
	for(i=0;i<DIMENSIONS;i++) {
		printf("Wymiar %d : ", i);
		scanf("%d %d", &LD[i], &RD[i]);
	}
	/*LD[0] = 0, LD[1] = 0, LD[2] = 0;
	RD[0] = 7, RD[1] = 7, RD[2] = 7;*/
	MDPoint* *Res = NULL;
	int ResSize =0;
	MTQuery(&Root,LD,RD,&Res,&ResSize, DIMENSIONS);
	printf("Znalezione punkty : \n");
	for(i=0;i<ResSize;i++) {
		printf("%d: ", i);
		int j;
		for(j=0;j<DIMENSIONS;j++) {
			printf("%d ",Res[i]->coordinates[j]); 
		}
		printf("\n");
	}
	free(Res);
	free(LD);
	free(RD);
	for(i=0;i<DATASIZE;i++)
		MDPointRemove(&Data[i]);
	free(Data);
	MTDelete(&Root);
	//MPI_Finalize();
}
