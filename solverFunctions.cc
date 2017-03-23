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

    for(k=0;k<n3;k++)
        for(j=0;j<n2;j++)
            for(i=0;i<n1;i++){
                
                for(n=0;n<6;n++)
                    dg[n]=strainbar[n]+straintilde[(k*n2*(n1)+j*(n1)+i)*6+n];
                
                for(n=0;n<6;n++){
                    x[n]=stress[(k*n2*(n1)+j*(n1)+i)*6+n];
                    for(m=0;m<6;m++)
                        x[n]+=C0_66[n][m]*dg[m];
                }


                for(m=0;m<6;m++){
                    edot[m]=0.0;
                    for(n=0;n<6;n++)
                        edot[m]+=fsloc[(k*n2*(n1)+j*(n1)+i)*36+6*m+n]*x[n];
                }


                for(m=0;m<6;m++){
                    ddg[m]=dg[m]-edot[m];                
                    dsg[m]=0.0;
                    for(n=0;n<6;n++)
                        dsg[m]+=C0_66[m][n]*(dg[m]-edot[m]);
                    stress[(k*n2*(n1)+j*(n1)+i)*6+m]+=dsg[m];
                }

                errstrain+=tnorm((double *)ddg,6,1)*volumeVoxel;
                errstress+=tnorm((double *)dsg,6,1)*volumeVoxel;

//                print1darray(stress[k*n2*(n1)+j*(n1)+i],6);
        }        
}

void findGammaHat(fourthOrderTensor Cref){

	int i,j,k,p,q,r,s,l,m;
    double fourierPoint[3],fourierPointNorm[3],fourierTensor[3][3],normPoint,det=0;
    double G[3][3];

    gammaHat=new fourthOrderTensor[n3*n2*(n1/2+1)];

    for(k=0;k<n3;k++){
        fourierPoint[2]=k;
        if(k>n3/2-1) fourierPoint[2]=(k-n3);
        for(j=0;j<n2;j++){
            fourierPoint[1]=j;
            if(j>n2/2-1) fourierPoint[1]=(j-n2);
            for(i=0;i<n1/2+1;i++){
                fourierPoint[0]=i;
                
                normPoint=tnorm((double *)fourierPoint,3,1);
                        
                fourierPointNorm[0]=fourierPoint[0];
                fourierPointNorm[1]=fourierPoint[1];
                fourierPointNorm[2]=fourierPoint[2];
                
                for(l=0;l<3;l++)
                    for(m=0;m<3;m++)
                    	fourierTensor[l][m]=fourierPointNorm[l]*fourierPointNorm[m];
			
				multiply3333x33((double *) G,Cref,fourierTensor,1,3);
	
				findInverse((double *)G,det,3);

				for(p=0;p<3;p++)
					for(q=0;q<3;q++)
						for(r=0;r<3;r++)
							for(s=0;s<3;s++)
								gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i].tensor[p][q][r][s]=G[p][r]*fourierTensor[q][s]; 

//                            std::cout<<i<<" "<<j<<" "<<k<<" "<<std::endl;
//                            print4darray(gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i].tensor);

            }
        }
    }

			    for(p=0;p<3;p++)
					for(q=0;q<3;q++)
						for(r=0;r<3;r++)
							for(s=0;s<3;s++)
								gammaHat[0].tensor[p][q][r][s]=0;


}