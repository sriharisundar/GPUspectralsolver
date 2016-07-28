#include "readInput.h"
#include "globalVariables.h"
#include "solverFunctions.h"
#include <iostream>
#include <fstream>
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

	cout<<outputFile<<endl;
	fstream fieldsOut;
	fieldsOut.open(outputFile,ios::out);

	fstream errorOut;
	errorOut.open("err.out",ios::out);

	//augmentLagrangian();
	return 0;
}