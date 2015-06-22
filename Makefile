HilbertLib: HilbertLib.c
	mpicc -o HilbertLib.o HilbertLib.c -g 
Run: HilbertLib.o
	mpiexec -n 10 HilbertLib.o
