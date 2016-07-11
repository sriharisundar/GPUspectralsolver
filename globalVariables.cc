#include "globalVariables.h"
#include "printFunctions.h"
#include <cmath>

int n1,n2,n3;

double cmat[6][6],cmat33[3][3][3][3];

double identityR2[3][3]={{1,0,0},{0,1,0},{0,0,1}};

double identityR4[3][3][3][3];

double basis[3][3][6];

double euler[N3][N2][N1][3];

int grainID[N3][N2][N1],phaseID[N3][N2][N1];

void initglobal(void){

	double rsq2=1.0/sqrt(2);		
	double rsq3=1.0/sqrt(3);		
	double rsq6=1.0/sqrt(6);		

	//Initialize 4th order kronecker delta
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++)
				for(int l=0;l<3;l++)
					identityR4[i][j][k][l]=(identityR2[i][k]*identityR2[j][l]+identityR2[i][l]*identityR2[j][k])/2.0;
			for(int k=0;k<3;k++)
				basis[i][j][k]=0.0;
		}

    basis[0][0][1]=-rsq6;
    basis[1][1][1]=-rsq6;
    basis[2][2][1]= 2.0*rsq6;

    basis[0][0][0]=-rsq2;
    basis[1][1][0]= rsq2;

    basis[1][2][2]=rsq2;
    basis[2][1][2]=rsq2;

    basis[0][2][3]=rsq2;
    basis[2][0][3]=rsq2;

    basis[0][1][4]=rsq2;
    basis[1][0][4]=rsq2;

    basis[0][0][5]=rsq3;
    basis[1][1][5]=rsq3;
    basis[2][2][5]=rsq3;

}