#include "readInput.h"
#include "matrixOperations.h"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

void readtexture(char filename[]){

	string line;

	fstream textureIn;
	textureIn.open(filename, ios::in);

	while(!textureIn.eof()){
	getline(textureIn,line);
	//cout<<line<<endl;
	}

	double a[3][3];
	double angles[][3]={{0,60,0},{45,0,0}};

	transformationMatrix(a,angles[0],2);

	cout<<a[2][2]<<endl;

}

void readprops(char filename[]){
	char line[100];

	fstream textureIn;
	textureIn.open(filename, ios::in);

	textureIn.getline(line,99);
	cout<< line <<endl;
	textureIn.getline(line,99);
	cout<< line <<endl;
	textureIn.getline(line,99);
	cout<< line <<endl;

}

void readinput(char filename[100],int N1, int N2, int N3){
	
	char textureFile[100],propsFile[100],dummy[100];

	fstream maininputIn;
	maininputIn.open(filename, ios::in);

	maininputIn.getline(dummy,100);
	maininputIn.getline(textureFile,100);
	maininputIn.getline(dummy,100);
	maininputIn.getline(propsFile,100);

	readtexture(textureFile);
	readprops(propsFile);
}
