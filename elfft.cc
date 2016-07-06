#include "readInput.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{	
	if (argc!=2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	readinput(argv[1]);

	return 0;
}