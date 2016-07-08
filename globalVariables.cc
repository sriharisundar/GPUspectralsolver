#include "globalVariables.h"

int n1,n2,n3;

double cmat[6][6],cmat33[3][3][3][3];

double identityR2[3][3]={{1,0,0},{0,1,0},{0,0,1}};

double identityR4[3][3][3][3];

double euler[N3][N2][N1][3];

int grainID[N3][N2][N1],phaseID[N3][N2][N1];

void initglobal(void){

	//Initialize 4th order kronecker delta
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				for(int l=0;l<3;l++)
					identityR4[i][j][k][l]=(identityR2[i][k]*identityR2[j][l]+identityR2[i][l]*identityR2[j][k])/2.0;


}