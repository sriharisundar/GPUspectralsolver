#include "globalVariables.h"
#include "printFunctions.h"
#include <cmath>
#include <string>

int n1,n2,n3;
double C0_66[6][6];
fourthOrderTensor C0_3333,cmat3333;
double identityR2[3][3];
double identityR4[3][3][3][3];
double basis[3][3][6];
double strainbar[6],stressbar[6];
double stressref,strainref,errstress,errstrain,error;
double RVEdim[3];
double velgrad33[3][3];
double straingradrate33[3][3],straingradrate6[6];
double rotationrate33[3][3];
double IDstraingradrate[6];
int nsteps,itermax;
int ictrl,ictrl1,ictrl2;
char *outputFile;

double *euler;
int *grainID,*phaseID;
double *ddefgrad,*ddefgradim;
double *straintilde,*stress,*delta;
double *work,*workim;
double *cloc,*fsloc;
fourthOrderTensor *gammaHat;

void initglobal(){

	int i,j,k,l;
	double rsq2=1.0/sqrt(2);		
	double rsq3=1.0/sqrt(3);		
	double rsq6=1.0/sqrt(6);		
	outputFile=new char[100];

	//Initialize 4th order kronecker delta
	for(i=0;i<3;i++)
		for(j=0;j<3;j++){
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
					identityR4[i][j][k][l]=(identityR2[i][k]*identityR2[j][l]+identityR2[i][l]*identityR2[j][k])/2.0;
			for(k=0;k<3;k++)
				basis[i][j][k]=0.0;
		}

	for(i=0;i<6;i++)
		for(j=0;j<6;j++)
			C0_66[i][j]=0;

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

	euler=new double[n3*n2*n1*3];
	grainID=new int[n3*n2*n1];
	phaseID=new int[n3*n2*n1];
	ddefgrad=new double[n3*n2*n1*9];
	ddefgradim=new double[n3*n2*n1*9];
	straintilde=new double[n3*n2*n1*6];
	stress=new double[n3*n2*n1*6];
	work=new double[n3*n2*n1*6];
	workim=new double[n3*n2*n1*6];
	cloc=new double[n3*n2*n1*36];
	fsloc=new double[n3*n2*n1*36];

}