#include "MyTree.h"
#include <assert.h>
#include <stdio.h>

void makeMTNode(MTNode *foo, int dimdiv, coord_t val) {
	foo->dim = dimdiv;
	foo->val = val;
	foo->left = NULL, foo->right = NULL;
}

void coordinatesMINMAX(MDPoint* *Data, int DataSize, int Dim, coord_t *MIN, coord_t* MAX) {
	int i;
	*MIN = Data[0]->coordinates[Dim];
	*MAX = Data[0]->coordinates[Dim];
	for(i=1;i<DataSize;i++) {
		if(Data[i]->coordinates[Dim] > *MAX) {
			*MAX = Data[i]->coordinates[Dim];
		}
		else if(Data[i]->coordinates[Dim] < *MIN) {
			*MIN = Data[i]->coordinates[Dim];	
		}
	}
}

void MTmake(
	MTNode *Node, 
	MDPoint* *Data, 
	int DataSize, 
	int Dimensions, 
	int ActDim
) {
	int i;
	printf("Data\n");
	for(i=0;i<DataSize;i++) {
		printf("%d:  ",i);
		int j;
		for(j=0;j<Dimensions;j++) {
			printf("%lld ", Data[i]->coordinates[j]);
		}
		printf("\n");
	}	
	if(DataSize == 0)
		return;
	if(DataSize == 1) {
		Node->right = NULL;
		Node->left = calloc(DataSize,sizeof(MDPoint*));
		int i;
		for(i=0;i<DataSize;i++) 
			((MDPoint**)(Node->left))[i] = Data[i];
		Node->val = DataSize;
		return;
	}
	MDPoint* *leftData = NULL;
	MDPoint* *rightData = NULL; 
	coord_t pivot;
	int countSmaller;
	coord_t MIN,MAX;
	int started = -1;
	while(1) {
		if(started == ActDim) { // Everything is the same
			Node->right = NULL;
			Node->left = calloc(DataSize,sizeof(MDPoint*));
			int i;
			for(i=0;i<DataSize;i++) 
				((MDPoint**)(Node->left))[i] = Data[i];
			Node->val = DataSize;
			return;
		}
		pivot = Data[rand()%DataSize]->coordinates[ActDim];
		countSmaller = 0;
		for(i=0;i<DataSize;i++) {
			if(Data[i]->coordinates[ActDim] <= pivot) {
				countSmaller++;
			}
		}
		if(countSmaller == 0 || DataSize-countSmaller == 0) {
			coordinatesMINMAX(Data,DataSize,ActDim,&MIN,&MAX);
			if(MIN != MAX) {
				pivot = (MIN + MAX)/2;
				break;
			}
			started = ActDim;
			ActDim+=1;
			ActDim%=Dimensions;
		} else {
			break;		
		}
	}
	printf("pivot : %d Node->dim = %d\n",pivot,ActDim);
	Node->dim = ActDim;
	Node->val = pivot;
	countSmaller = 0;
	for(i=0;i<DataSize;i++) {
		if(Data[i]->coordinates[ActDim] <= pivot) {
			countSmaller++;
		}
	}	
	assert(countSmaller != 0);
	assert((DataSize - countSmaller) != 0);
	leftData = calloc(countSmaller,sizeof(MDPoint*));
	rightData = calloc(DataSize-countSmaller,sizeof(MDPoint*));
	int leftDataPtr = 0, rightDataPtr = 0;
	for(i=0;i<DataSize;i++) {
		if(Data[i]->coordinates[ActDim] <= pivot) {
			leftData[leftDataPtr] = Data[i];
			leftDataPtr+=1;
		} else {
			rightData[rightDataPtr] = Data[i];
			rightDataPtr+=1;
		}
	}
	Node->left = calloc(1,sizeof(MTNode));
	Node->right = calloc(1,sizeof(MTNode));
	makeMTNode(Node->left,0,0);
	makeMTNode(Node->right,0,0);
	MTmake(Node->left,leftData,countSmaller,Dimensions,(ActDim+1)%Dimensions);
	free(leftData);
	MTmake(Node->right,rightData,DataSize-countSmaller,Dimensions,(ActDim+1)%Dimensions);
	free(rightData);
}

void MTDelete(MTNode *Node) {
	if(Node->left != NULL && Node->right != NULL) {
		MTDelete(Node->left);
		MTDelete(Node->right);
		free(Node->left);
		free(Node->right);
	} else {
		//free((MDPoint**)Node->left);
	}
}

void MTQueryLocal(
	MTNode *Node,
	coord_t* LD,
	coord_t* RD,
	PtrVector *vec,
	int Dimensions
) {
	if(Node->right == NULL) {
		printf("Koncze\n");
		if(Node->left == NULL) 
			return;
		int i;
		printf("Znalazlem ");
		for(i=0;i<Dimensions;i++) {
			coord_t x = (((MDPoint**)(Node->left))[0])->coordinates[i];
			if(x < LD[i] || x > RD[i])
				return;
			printf("%lld \n",x);
		}
		for(i=0;i<(Node->val);i++) {
			PtrVectorPB(vec,((MDPoint**)Node->left)[i]);
		}
		return;
	}
	if(LD[Node->dim] <= Node->val) {
		printf("Szukam w lewej (%dta wspolredna jest <= %d\n",Node->dim,Node->val);
		MTQueryLocal(Node->left,LD,RD,vec,Dimensions);
	}
	if(RD[Node->dim] > Node->val) {
		printf("Szukam w prawej (%dta wspolredna jest > %d\n",Node->dim,Node->val);
		MTQueryLocal(Node->right,LD,RD,vec,Dimensions);
	}
}

void MTQuery(
	MTNode *Node, 
	coord_t* LD, 
	coord_t* RD, 
	MDPoint** *Res, 
	int* ResSize,
	int Dimensions
) {
	PtrVector vec;
	makePtrVector(&vec);
	MTQueryLocal(
		Node,
		LD,
		RD,
		&vec,
		Dimensions
	);
	*Res = calloc(vec.size,sizeof(MDPoint*));
	*ResSize = vec.size;
	printf("vec.size = %d\n", vec.size) ;
	memcpy(*Res,vec.arr,(*ResSize)*sizeof(MDPoint*));
}
	

