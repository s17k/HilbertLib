HilbertLibVT: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc mpicc -o HilbertLib.o HilbertLib.c 
MainDBG: Main.c
	mpicc -c MDPoint.c -o MDPoint.o -g
	mpicc -c AxesTranspose.c -o AxesTranspose.o -g
	mpicc -c HilbertLib.c -o HilbertLib.o -g 
	mpicc -c Main.c -o Main.o -g
	mpicc -c MyTree.c -o MyTree.o
	mpicc -c PtrVector.c -o PtrVector.o
	mpicc MDPoint.o AxesTranspose.o HilbertLib.o Main.o -o Main.x -g 
MyTree: MyTreeTest.c
	mpicc -c MDPoint.c -o MDPoint.o -g
	mpicc -c AxesTranspose.c -o AxesTranspose.o -g
	mpicc -c MyTree.c -o MyTree.o -g
	mpicc -c PtrVector.c -o PtrVector.o -g
	mpicc -c MyTreeTest.c -o MyTreeTest.o -g
	mpicc MDPoint.o AxesTranspose.o PtrVector.o MyTree.o MyTreeTest.o -o MyTreeTest.x -g
clean:
	rm *.o
	rm *.x
HilbertLibFAST: HilbertLib.c
	/home/administrator/vt/bin/vtcc -o HilbertLib.o HilbertLib.c -O2
Run: HilbertLib.o
	mpiexec -n 10 HilbertLib.o
