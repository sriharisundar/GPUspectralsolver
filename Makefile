CC=g++ -std=c++11
CFLAGS=-I. -Wall -lfftw3 -lm
DEPS = matrixOperations.h readInput.h globalVariables.h printFunctions.h solverFunctions.h matrix.h
OBJ = matrixOperations.o readInput.o globalVariables.o printFunctions.o solverFunctions.o matrix.o

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

elfft: elfft.o $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

testinverse: testInverse.o $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

all: elfft
	make elfft