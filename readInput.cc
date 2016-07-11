#include "readInput.h"
#include "matrixOperations.h"
#include "globalVariables.h"
#include "printFunctions.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h> //for atoi()

void readtexture(char filename[]){

	int i,j,k,gid,pid;
	double a[3][3];
	double phi1,Phi,phi2;
	double cvoxel[3][3][3][3];

	std::string line;

	std::fstream textureIn;
	textureIn.open(filename, std::ios::in);

	while(!textureIn.eof()){
		textureIn>>phi1>>Phi>>phi2>>i>>j>>k>>gid>>pid;
		euler[k-1][j-1][i-1][0]=phi1;
		euler[k-1][j-1][i-1][1]=Phi;
		euler[k-1][j-1][i-1][2]=phi2;
		grainID[k-1][j-1][i-1]=gid;
		phaseID[k-1][j-1][i-1]=pid;

		transformationMatrix(a,euler[k-1][j-1][i-1],2);
		transformFourthOrderTensor(cmat33,cvoxel,a);
	}

}

void readprops(char filename[]){
	char line[100];
	int i,j,k,l;
	double youngs,poisson,mu,lambda;

	std::fstream textureIn;
	textureIn.open(filename, std::ios::in);

	textureIn.getline(line,99);
	if(atoi(line)==0){
		std::cout<<"Not isotropic"<<std::endl;

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
				textureIn>>cmat[i][j];

		voigt(cmat,cmat33,1);
	}

	else{
		std::cout<<"Isotropic"<<std::endl;
		textureIn>>youngs>>poisson;

		mu=youngs/(2.0*(1+poisson));
		lambda=2.0*mu*poisson/(1.0-2.0*poisson);

		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					for(l=0;l<3;l++)
						cmat33[i][j][k][l]=lambda*identityR2[i][j]*identityR2[k][l]+2.0*mu*identityR4[i][j][k][l];

	}
}

void readinput(char filename[100]){
	
	char textureFile[100],propsFile[100],dummy[100];
	int n1,n2,n3;

	std::fstream maininputIn;
	maininputIn.open(filename, std::ios::in);

	maininputIn.getline(dummy,100);
	maininputIn.getline(textureFile,100);

	maininputIn.getline(dummy,100);
	maininputIn.getline(propsFile,100);
	
	maininputIn.getline(dummy,100);
	maininputIn>>n1>>n2>>n3;

	readprops(propsFile);
	readtexture(textureFile);


//	print2darray(cmat);
//	print4darray(cmat33);
}
