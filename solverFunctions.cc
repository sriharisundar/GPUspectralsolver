#include "printFunctions.h"
#include "globalVariables.h"
#include "matrixOperations.h"

void augmentLagrangian(void){

	double x[6],dg[6],edot[6],dsg[6],ddg[6];
	double errd=0.0,errs=0.0;

	int i,j,k,n,m;

	for(k=0;k<N3;k++)
		for(j=0;j<N2;j++)
			for(i=0;i<N1;i++){
				
				for(n=0;n<6;n++)
					dg[n]=dbar[n]+dtilde[k][j][i][n];
				
				for(n=0;n<6;n++){
					x[n]=sg[k][j][i][n];
					for(m=0;m<6;m++)
						x[m]+=xlsec66[n][m]*dg[m];
				}

				for(n=0;n<6;n++){
					edot[n]=0.0;
					for(m=0;m<6;m++)
						edot[n]+=fsloc[k][j][i][n][m]*x[m];
				}

				for(n=0;n<6;n++){
					ddg[n]=dg[n]-edot[n];				
					dsg[n]=0.0;
					for(m=0;m<6;m++)
						dsg[n]+=xlsec66[n][m]*(dg[n]-edot[n]);
					sg[k][j][i][n]+=dsg[n];
				}

				errd+=tnorm((double *)ddg,6,1)*wgt;
				errs+=tnorm((double *)dsg,6,1)*wgt;

		}		
}