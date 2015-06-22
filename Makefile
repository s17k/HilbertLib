HilbertLibDBG: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc mpicc -o HilbertLib.o HilbertLib.c -O2
HilbertLibFAST: HilbertLib.c
	/home/administrator/vt/bin/vtcc -o HilbertLib.o HilbertLib.c -O2
Run: HilbertLib.o
	mpiexec -n 10 HilbertLib.o
