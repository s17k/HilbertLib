HilbertLibVT: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc mpicc -o HilbertLib.o HilbertLib.c 
MainDBG: Main.c
	cc -c MDPoint.c -o MDPoint.o -g
	cc -c AxesTranspose.c -o AxesTranspose.o -g
	mpicc -c HilbertLib.c -o HilbertLib.o -g 
	mpicc -c Main.c -o Main.o -g
	cc -c MyTree.c -o MyTree.o -g
	cc -c PtrVector.c -o PtrVector.o -g
	cc -c Pair.c -o Pair.o -g
	mpicc MDPoint.o AxesTranspose.o HilbertLib.o Main.o MyTree.o Pair.o PtrVector.o -o Main.x -g 
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
