#include "readInput.h"
#include "matrixOperations.h"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

void readinput(char filename[100]){
	
	char textureFile[100],propsFile[100],dummy[100];
	int matrix[6][2]={{1,1},{2,2},{3,3},{2,3},{1,3},{1,2}};

	cout<<matrix[1][1]<<endl;

	fstream maininputIn;
	maininputIn.open(filename, ios::in);

	maininputIn.getline(dummy,100);
	maininputIn.getline(textureFile,100);
	maininputIn.getline(dummy,100);
	maininputIn.getline(propsFile,100);

	readtexture(textureFile);
	readprops(propsFile);
}

void readtexture(char filename[]){

	string line;

	fstream textureIn;
	textureIn.open(filename, ios::in);

	while(!textureIn.eof()){
	getline(textureIn,line);
	//cout<<line<<endl;
	}

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