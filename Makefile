include ./CBLAS/Makefile.in

CC=g++
CFLAGS=-I. -Wall
DEPS = matrixOperations.h readInput.h globalVariables.h printFunctions.h solverFunctions.h
OBJ = elfft.o matrixOperations.o readInput.o globalVariables.o printFunctions.o solverFunctions.o

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

elfft: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
