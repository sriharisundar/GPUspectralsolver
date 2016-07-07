#include "readInput.h"
#include <iostream>
using namespace std;

#define N1 32
#define N2 32
#define N3 32

int main(int argc, char *argv[])
{	
	if (argc!=2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	readinput(argv[1],N1,N2,N3);

	return 0;
}