CC=g++
CFLAGS=-I. -Wall
DEPS = matrixOperations.h readInput.h globalVariables.h printFunctions.h
OBJ = elfft.o matrixOperations.o readInput.o globalVariables.o printFunctions.o

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

elfft: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)