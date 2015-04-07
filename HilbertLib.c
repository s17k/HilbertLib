#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

typedef unsigned int coord_t;

typedef unsigned int coord_t;
void AxestoTranspose( coord_t* X, int b, int n) { // position,bits,dimensions
	coord_t M = 1 << (b-1), P, Q, t;
	int i,j;
	for(Q=M;Q>1;Q>>=1) {
		P = Q -1;
		for(i=0;i<n;i++)
			if(X[i]&Q) 
				X[0]^=P;
			else 
				{t = (X[0]^X[i]) & P;X[0] ^= t; X[i] ^= t; }
	}
	for( i=1;i<n;i++) 
		X[i] ^= X[i-1];
	t = 0;
	for(Q=M;Q>1;Q>>=1)
		if(X[n-1] & Q)
			t ^= Q-1;
	for(i=0;i<n;i++) X[i]^=t;
	coord_t *Y = calloc(1,sizeof(coord_t));
	Y[0] = 0;
	for(i=b-1;i>=0;i--) {
		for(j=0;j<n;j++) {
			Y[0] *= 2;
			Y[0] += ((X[j]>>i)&1);
		}
	}
	X[0] = Y[0];
	free(Y);
}

typedef unsigned int uint;
typedef long long ll;


// It sends all points' coordinates from each process to root
const int HilbertSendDataRank = 101; // Used rank (maybe +1) 
void HilbertSendData(int root, MDPoint *particles, ll *wages, uint particles_num, uint dimensions) {
	double *buff;
	buff = calloc(particles_num*dimensions,sizeof(double));

	uint i=0;
	uint wsk = 0;
	uint j=0;
	for(i=0;i<particles_num;i++) {
		for(j=0;j<dimensions;j++) {
			buff[wsk] = particles[i]->c[j];
			wsk++;
		}
	}
	int k = malloc(sizeof(int));
	k = particles_num;
	if(wages != NULL)
		k *= -1;
	// particles number sending and waiting
	MPI_Request particlesNumMPIRequest;
	MPI_Isend(&k,1,MPI_INT,root,HilbertSendDataRank,MPI_COMM_WORLD,&particlesNumMPIRequest);
	MPI_Status particlesNumStatus;
	MPI_Wait(&particlesNumMPIRequest,&particlesNumStatus);
	
	//particles sending and waiting
	MPI_Request particlesRequest;
	MPI_Isend(buff, particles_num*dimensions, MPI_LONG_LONG, root, HilbertSendDataRank+1,MPI_COMM_WORLD, &particlesRequest);
	MPI_Status particlesStatus;
	MPI_Wait(&particleRequest, &particleStatus);
	
	free(buff);
	free(k);

	if(wages != NULL) {
		ll *wagesbuff = calloc(particles_num,sizeof(ll));
		memcpy(wagesbuff,wages,particles_num);
		MPI_Request wagesRequest;
		MPI_Isend(wagesbuff, particles_num, root, 103,MPI_COMM_WORLD, &wagesRequest);
		MPI_Status wagesStatus;
		MPI_Wait(&wagesRequest, &wagesStatus);
		free(wagesbuff);
	}
}

int abs(int x) {
	return x > 0 ? x : -x;
}

struct HilbertData {
	MDPoint *particles;
	ll *wages;
	uint particles_num;
};

typedef struct HilbertData HilbertData;

void HilbertRecvData(int guy, uint dimensions, HilbertData *data) {
	int receiveCount;
	MPI_Request receiveCountRequest;
	MPI_Irecv(&receiveCount,1, MPI_INT, guy,101, MPI_COMM_WORLD, &receiveCountRequest);
	MPI_Status st;
	MPI_Wait(&receiveCountRequest,&st);

	bool ifwages = false;
	if(receiveCount < 0)
		ifwages = true;
	double *buff = callloc(dimensions*receiveCount,sizeof(double));
	MPI_Request buffRequest;
	MPI_Irecv(buff,dimensions*receiveCount,MPI_DOUBLE,guy,102,MPI_COMM_WORLD,&buffRequest);
	data->wages = NULL;
	if(ifwages) {
		MPI_Request wagesBuffRequest;
		MPI_Status wagesStatus;
		ll *wagesBuff = calloc(receiveCount,sizeof(long long));
		MPI_Irecv(wagesBuff,receiveCount,MPI_LONG_LONG,guy,103,MPI_COMM_WORLD,&wagesBuffRequest);
		MPI_Wait(&wagesBuffRequest,&wagesStatus);
		data->wages = wagesBuff;
		free(wagesBuff);
	}
	MPI_Wait(&buffRequest,&st);
	int i=0;
	data->particles = calloc(receiveCount,sizeof(MDPoint));
	int ptr = 0;
	for(i=0;i<receiveCount;i++) {
		uint j;
		data->particles->c = calloc(dimensions,sizeof(double));
		for(j=0;j<dimensions;j++) {
			data->particles[i]->c[j] = buff[ptr];
			ptr++;
		}
	}
	data->particles_num = receiveCount;
	free(buff);
}

void Partition(HilbertData *data) {
	
}
