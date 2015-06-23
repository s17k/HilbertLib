HilbertLibVT: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc mpicc -o HilbertLib.o HilbertLib.c 
MainDBG: Main.c
	mpicc -c MDPoint.c -o MDPoint.o -g
	mpicc -c AxesTranspose.c -o AxesTranspose.o -g
	mpicc -c HilbertLib.c -o HilbertLib.o -g 
	mpicc -c Main.c -o Main.o -g
	mpicc MDPoint.o AxesTranspose.o HilbertLib.o Main.o -o Main.x -g 
clean:
	rm *.o
	rm *.x
HilbertLibFAST: HilbertLib.c
	/home/administrator/vt/bin/vtcc -o HilbertLib.o HilbertLib.c -O2
Run: HilbertLib.o
	mpiexec -n 10 HilbertLib.o
