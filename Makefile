HilbertLibVT: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc mpicc -o HilbertLib.o HilbertLib.c 
CC=mpicc
OBJ=MDPoint.o MyTree.o PtrVector.o Pair.o HilbertLib.o Main.o

CFLAGS=-O3 -Winline
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
main: $(OBJ)	
	$(CC) $(OBJ) -o Main.x $(CFLAGS)

MyTreeOBJ=MDPoint.o MyTree.o PtrVector.o MyTreeTest.o
MyTree: $(MyTreeOBJ)
	$(CC) $(MyTreeOBJ) -o MyTreeTest.x $(CFLAGS)
clean:
	rm -f *.o
	rm -f *.x
