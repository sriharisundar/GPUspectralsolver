CC=g++ -std=c++11
CFLAGS=-Wall -lfftw3 -lm -lclFFT -lOpenCL
HEADERS=-I. -I $(HOME)/clFFT-2.12.2-Linux-x64/include -I $(HOME)/fftw/include/
LIBS= -L. -L $(HOME)/clFFT-2.12.2-Linux-x64/lib64 -L $(HOME)/fftw/lib/
DEPS = matrixOperations.h readInput.h globalVariables.h printFunctions.h solverFunctions.h
OBJ = matrixOperations.o readInput.o globalVariables.o printFunctions.o solverFunctions.o 

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(HEADERS) $(LIBS) 

fftw: $(OBJ)
	$(CC) -DUSEFFTW=1 -c elfft.cc $(CFLAGS) $(HEADERS) $(LIBS) 
	$(CC) -o elfft_fftw elfft.o $^ $(CFLAGS) $(HEADERS) $(LIBS) 

clfft: $(OBJ)
	$(CC) -DUSECLFFT=1 -c elfft.cc $(CFLAGS) $(HEADERS) $(LIBS) 
	$(CC) -o elfft_clfft elfft.o $^ $(CFLAGS) $(HEADERS) $(LIBS) 

all: elfft
	make fftw

clean: 
	rm *.o
	rm elfft_fftw
	rm elfft_clfft
