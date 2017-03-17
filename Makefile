CC=g++ -std=c++11
CFLAGS=-Wall -lfftw3 -lm
HEADERS=-I. 
LIBS= -L.
DEPS = matrixOperations.h readInput.h globalVariables.h printFunctions.h solverFunctions.h
OBJ = matrixOperations.o readInput.o globalVariables.o printFunctions.o solverFunctions.o 

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(HEADERS) $(LIBS) 

elfft: elfft.o $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(HEADERS) $(LIBS)

all: elfft
	make elfft

clean: 
	rm *.o
	rm elfft
