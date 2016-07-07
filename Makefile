CC=g++
CFLAGS=-I. -Wall
DEPS = matrixOperations.h readInput.h
OBJ = elfft.o matrixOperations.o readInput.o 

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

elfft: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)