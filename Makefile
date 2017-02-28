CC=g++ -std=c++11
CFLAGS=-Wall -lfftw3 -lm -lclFFT -lOpenCL
HEADERS=-I. -I $(HOME)/clFFT-2.12.2-Linux-x64/include
LIBS= -L. -L $(HOME)/clFFT-2.12.2-Linux-x64/lib64 
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
