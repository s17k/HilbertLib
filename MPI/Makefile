HilbertLibVT: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc mpicc -o HilbertLib.o HilbertLib.c 
CC=mpicc
OBJ=MDPoint.o MyTree.o PtrVector.o Pair.o HilbertLib.o Main.o

CFLAGS=-O3 -Winline
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
main: $(OBJ)	
	$(CC) $(OBJ) -o Main.x $(CFLAGS)

MyTreeTestOBJ=MDPoint.o MyTree.o PtrVector.o tests/MyTreeTest.o
MyTree: $(MyTreeTestOBJ)
	$(CC) $(MyTreeTestOBJ) -o MyTreeTest.x $(CFLAGS)
MyTreeAutoTestOBJ=MDPoint.o MyTree.o PtrVector.o MyTreeAutoTest.o
MyTreeAuto: $(MyTreeAutoTestOBJ)
	$(CC) $(MyTreeAutoTestOBJ) -o MyTreeTest.x $(CFLAGS)

	
clean:
	rm -rf *.o
	rm -rf *.x
