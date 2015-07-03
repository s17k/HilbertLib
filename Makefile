HilbertLibVT: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc mpicc -o HilbertLib.o HilbertLib.c 
OBJ=MDPoint.o MyTree.o PtrVector.o Pair.o HilbertLib.o Main.o

CFLAGS=-O3 -Winline
%.o: %.c
	mpicc -c -o $@ $< $(CFLAGS)
MainDBG: $(OBJ)	
	mpicc $(OBJ) -o Main.x $(CFLAGS)
MainMPICH: Main.c MDPoint.c AxesTranspose.c MyTree.c PtrVector.c Pair.c 
	mpicc.mpich2 -c MDPoint.c -o MDPoint.o -g
	mpicc.mpich2 -c AxesTranspose.c -o AxesTranspose.o -g
	mpicc.mpich2 -c HilbertLib.c -o HilbertLib.o -g 
	mpicc.mpich2 -c Main.c -o Main.o -g
	mpicc.mpich2 -c MyTree.c -o MyTree.o -g
	mpicc.mpich2 -c PtrVector.c -o PtrVector.o -g
	mpicc.mpich2 -c Pair.c -o Pair.o -g
	mpicc.mpich2 MDPoint.o AxesTranspose.o HilbertLib.o Main.o MyTree.o Pair.o PtrVector.o -o Main.x -g 


MyTree: MyTreeTest.c
	cc -c MDPoint.c -o MDPoint.o -g
	cc -c AxesTranspose.c -o AxesTranspose.o -g
	cc -c MyTree.c -o MyTree.o -g
	cc -c PtrVector.c -o PtrVector.o -g
	cc -c MyTreeTest.c -o MyTreeTest.o -g
	cc MDPoint.o AxesTranspose.o PtrVector.o MyTree.o MyTreeTest.o -o MyTreeTest.x -g
clean:
	rm *.o
	rm *.x
HilbertLibFAST: HilbertLib.c
	/home/administrator/vt/bin/vtcc -o HilbertLib.o HilbertLib.c -O2
Run: HilbertLib.o
	mpiexec -n 10 HilbertLib.o
	mpicc -c MDPoint.c -o MDPoint.o $(CFLAGS)
	mpicc -c AxesTranspose.c -o AxesTranspose.o $(CFLAGS)
	mpicc -c HilbertLib.c -o HilbertLib.o $(CFLAGS)
	mpicc -c Main.c -o Main.o $(CFLAGS)
	mpicc -c MyTree.c -o MyTree.o $(CFLAGS)
	mpicc -c PtrVector.c -o PtrVector.o $(CFLAGS)
	mpicc -c Pair.c -o Pair.o $(CFLAGS)

