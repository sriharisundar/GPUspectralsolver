#include "printFunctions.h"
#include "globalVariables.h"
#include "matrixOperations.h"

void augmentLagrangian(void){

	double x[6],dg[6],edot[6],dsg[6],ddg[6];
	errstrain=errstress=0.0;

	int i,j,k,n,m;

	for(k=0;k<N3;k++)
		for(j=0;j<N2;j++)
			for(i=0;i<N1;i++){
				
				for(n=0;n<6;n++)
					dg[n]=strainbar[n]+straintilde[k][j][i][n];
				
				for(n=0;n<6;n++){
					x[n]=stress[k][j][i][n];
					for(m=0;m<6;m++)
						x[m]+=C0_66[n][m]*dg[m];
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
						dsg[n]+=C0_66[n][m]*(dg[n]-edot[n]);
					stress[k][j][i][n]+=dsg[n];
				}

				errstrain+=tnorm((double *)ddg,6,1)*wgt;
				errstress+=tnorm((double *)dsg,6,1)*wgt;

		}		
}

void findGammaHat(void){

	int i,j,k;
	double fourierPoint[3];
	gammaHat=new fourthOrderTensor[N3*N2*(N1/2+1)];

	for(k=0;k<N3;k++){
		fourierPoint[2]=k;
		if(k>N3/2) fourierPoint[2]=k-N3;
		for(j=0;j<N2;j++){
			fourierPoint[1]=j;
			if(j>N2/2) fourierPoint[2]=j-N2;
			for(i=0;i<N1/2+1;i++){
				fourierPoint[0]=i;

			}
		}
	}

}