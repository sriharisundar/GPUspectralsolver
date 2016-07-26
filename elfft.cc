#include "readInput.h"
#include "globalVariables.h"
#include <iostream>
#include <fftw3.h>
using namespace std;

int main(int argc, char *argv[])
{	
	if (argc!=2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	initglobal();
	readinput(argv[1]);

	return 0;
}