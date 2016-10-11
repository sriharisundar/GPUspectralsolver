#include "printFunctions.h"
#include "globalVariables.h"
#include "matrixOperations.h"
#include <math.h>
#include <iostream>

using namespace std;

void augmentLagrangian(void){
	
	double prodDim=n1*n2*n3;
	double volumeVoxel=1.0/prodDim;
	double x[6],dg[6],edot[6],dsg[6],ddg[6];

    errstrain=errstress=0.0;

    int i,j,k,n,m;

    for(k=0;k<N3;k++)
        for(j=0;j<N2;j++)
            for(i=0;i<N1;i++){
                
                for(n=0;n<6;n++)
                    dg[n]=strainbar[n]+straintilde[k][j][i][n];
				//print1darray((double *)dg,6);
                
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

                errstrain+=tnorm((double *)ddg,6,1)*volumeVoxel;
                errstress+=tnorm((double *)dsg,6,1)*volumeVoxel;

        }        
}

void findGammaHat(fourthOrderTensor Cref){

	int i,j,k,p,q,r,s,l,m;
    double fourierPoint[3],fourierPointNorm[3],fourierTensor[3][3],normPoint,det=0;
    double G[3][3];

    gammaHat=new fourthOrderTensor[n3*n2*(n1)];

    for(k=0;k<n3;k++){
        fourierPoint[2]=k/(n3*RVEdim[2]);
        if(k>N3/2) fourierPoint[2]=(k-n3)/(n3*RVEdim[2]);
        for(j=0;j<n2;j++){
            fourierPoint[1]=j/(n2*RVEdim[1]);
            if(j>n2/2) fourierPoint[1]=(j-n2)/(n2*RVEdim[1]);
            for(i=0;i<N1;i++){
                fourierPoint[0]=i/(n1*RVEdim[0]);
                if(i>n1/2) fourierPoint[0]=(i-n1)/(n1*RVEdim[0]);
                
                normPoint=tnorm((double *)fourierPoint,3,1);
                        
                fourierPointNorm[0]=fourierPoint[0]/normPoint;
                fourierPointNorm[1]=fourierPoint[1]/normPoint;
                fourierPointNorm[2]=fourierPoint[2]/normPoint;
                
                for(l=0;l<3;l++)
                    for(m=0;m<3;m++)
                    	fourierTensor[l][m]=fourierPointNorm[l]*fourierPointNorm[m];
			
				multiply3333x33(G,Cref,fourierTensor,2,4);
	
				findInverse((double *)G,det,3);

				for(p=0;p<3;p++)
					for(q=0;q<3;q++)
						for(r=0;r<3;r++)
							for(s=0;s<3;s++)
								gammaHat[k*n2*(n1)+j*(n1)+i].tensor[p][q][r][s]=-1*G[p][r]*fourierTensor[q][r];                
            }
        }
    }

			    for(p=0;p<3;p++)
					for(q=0;q<3;q++)
						for(r=0;r<3;r++)
							for(s=0;s<3;s++)
								gammaHat[0].tensor[p][q][r][s]=0;


}