#ifndef PTRVECTORDEFINED
#define PTRVECTORDEFINED
struct PtrVector
{
	void **arr;
	int capacity;
	int size;
};

typedef struct PtrVector PtrVector;

void makePtrVector (PtrVector *);
void PtrVectorPB (PtrVector *, void *);
void PtrVectorDeallocate (PtrVector *);
#endif
