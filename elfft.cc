#include "readInput.h"
#include "globalVariables.h"
#include "solverFunctions.h"
#include <iostream>
#include <fstream>
#include <fftw3.h>
using namespace std;

int main(int argc, char *argv[])
{	
	int i,j,k,n,m;
	double stressbar33[3][3],aux66[6][6],aux3333[3][3][3][3];
	int iteration;

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

	//Initialize stress and strain fields
	for(k=0;k<N3;k++)
		for(j=0;j<N2;j++)
			for(i=0;i<N1;i++){

				for(n=0;n<6;n++)
					straintilde[k][j][i][n]=0;

				for(n=0;n<6;n++){
					if(phaseID[k][j][i]==2)
						stress[k][j][i][n]=0;
					else{

						for(m=0;m<6;m++)
							stress[k][j][i][n]=cloc[k][j][i][n][m]*straintilde[k][j][i][m];
					}
					stressbar[n]+=stress[k][j][i][n]*wgt;
				}
			}

	change_basis(stressbar,stressbar33,aux66,aux3333,1);
	stressref=stressbar33[ictrl1][ictrl2];

	while(step<=nsteps){
		iter
	}

	augmentLagrangian();
	return 0;
}