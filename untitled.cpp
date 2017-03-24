#define __CL_ENABLE_EXCEPTIONS
#include "readInput.h"
#include "globalVariables.h"
#include "solverFunctions.h"
#include "matrixOperations.h"
#include "printFunctions.h"
#include <iostream>
#include <fstream>
#include <fftw3.h>
using namespace std;

int main(int argc, char *argv[])
{	
	int i,j,k,n,m;
	double stressbar33[3][3],strainbar33[3][3],work33[3][3],work33im[3][3];
	double aux33[3][3],aux66[6][6],aux3333[3][3][3][3];
	double strain[6],stress6[6];
	double strainout[3][3],stressout[3][3];
	double err2mod;
	double prodDim;
	double volumeVoxel;
	int iteration,step;
	double *delta;

	fftw_complex *out;
	fftw_plan plan_backward;
	fftw_plan plan_forward;


	if (argc<2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	readinput(argv[1]);

	prodDim=n1*n2*n3;
	volumeVoxel=1.0/prodDim;

	delta=new double[n3*n2*n1];
	out=(fftw_complex *) *fftw_alloc_complex(n3*n2*(n1/2+1));

	plan_forward=fftw_plan_dft_r2c_3d(N3,N2,N1,delta,out,FFTW_ESTIMATE);
	plan_backward=fftw_plan_dft_c2r_3d ( N3, N2, N1, out, delta, FFTW_ESTIMATE );	

	cout<<"Output file:"<<outputFile<<endl;
	fstream fieldsOut;
	fieldsOut.open(outputFile,ios::out);

	fstream errorOut;
	errorOut.open("err.out",ios::out);

	// Initialize stress and strain fields
	// sg - stress
	// dtilde-straingradient - straintilde
	
	for(k=0;k<n3;k++)
		for(j=0;j<n2;j++)
			for(i=0;i<n1;i++){

				for(n=0;n<6;n++)
					straintilde[(k*n2*(n1)+j*(n1)+i)*6+n]=0;

				for(n=0;n<6;n++){
					if(phaseID[k*n2*(n1)+j*(n1)+i]==2)
						stress[(k*n2*(n1)+j*(n1)+i)*6+n]=0;
					else{
						stress[(k*n2*(n1)+j*(n1)+i)*6+n]=0;
						for(m=0;m<6;m++){
							stress[(k*n2*(n1)+j*(n1)+i)*6+n]+=cloc[(k*n2*(n1)+j*(n1)+i)*36+n*6+m]*(strainbar[m]);
						}
					}
					stressbar[n]+=stress[(k*n2*(n1)+j*(n1)+i)*6+n]*volumeVoxel;
				}
	}

	change_basis(stressbar,stressbar33,aux66,aux3333,1);
	stressref=stressbar33[ictrl1][ictrl2];

	for(step=1;step<=nsteps;step++){
		iteration=0;

		err2mod=2*error;

		findGammaHat(C0_3333);

		while(iteration<itermax && err2mod > error){
			iteration++;
			cout<<"--------------------------------------------------------------"<<endl;
			cout<<"ITERATION:"<<iteration<<endl;

			// Arrange data for in
			// Perform forward FFT
			cout<<"Forward FFT of polarization field"<<endl<<endl;

			for(n=0;n<6;n++){
				
				for(k=0;k<n3;k++)
					for(j=0;j<n2;j++)
						for(i=0;i<n1;i++){
							delta[(k*n2*(n1)+j*(n1)+i)]=stress[(k*n2*(n1)+j*(n1)+i)*6+n];
							for(m=0;m<6;m++)
								delta[(k*n2*(n1)+j*(n1)+i)]-=C0_66[n][m]*straintilde[(k*n2*(n1)+j*(n1)+i)*6+m];
				}

				fftw_execute(plan_forward);


				for(k=0;k<n3;k++)
					for(j=0;j<n2;j++)
						for(i=0;i<(n1/2+1);i++){
							work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][0];
							workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][1];							
				}

			}

			// Convert stress to tensorial form
			// Multiply with gamma operator
			cout<<"Gamma convolution"<<endl<<endl;
			for(k=0;k<n3;k++)
				for(j=0;j<n2;j++)
					for(i=0;i<(n1/2+1);i++){
						change_basis(&work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6],work33,aux66,aux3333,1);
						change_basis(&workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6],work33im,aux66,aux3333,1);

						multiply3333x33(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i],work33,3,4);
						multiply3333x33(&ddefgradim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i],work33im,3,4);
			}
			
			change_basis(strainbar,strainbar33,aux66,aux3333,1);

			// Arrange data for out
			cout<<"Inverse FFT to get deformation gradient"<<endl<<endl;
			for(m=0;m<3;m++)
				for(n=0;n<3;n++){
					for(k=0;k<n3;k++)
						for(j=0;j<n2;j++)
							for(i=0;i<(n1/2+1);i++){
								out[k*n2*(n1/2+1)+j*(n1/2+1)+i][0]=ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9+3*m+n];
								out[k*n2*(n1/2+1)+j*(n1/2+1)+i][1]=ddefgradim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9+3*m+n];
					}

					out[0][0]=strainbar33[m][n];

					fftw_execute(plan_backward);

					for(k=0;k<n3;k++)
						for(j=0;j<n2;j++)
							for(i=0;i<n1;i++){
								ddefgrad[(k*n2*(n1)+j*(n1)+i)*9+3*m+n]=delta[(k*n2*(n1)+j*(n1)+i)]/prodDim;
							}

			}

			// Get symmetric part of defgrad		
			for(k=0;k<n3;k++)
				for(j=0;j<n2;j++)
					for(i=0;i<n1;i++){
						symmetric(&ddefgrad[(k*n2*(n1)+j*(n1)+i)*9],(double *) aux33);
						change_basis(&straintilde[(k*n2*(n1)+j*(n1)+i)*6],aux33,aux66,aux3333,2);
			}


			cout<<"Augmented Lagrangian method for stress update"<<endl<<endl;

			augmentLagrangian();
			
			for(n=0;n<6;n++)
				stressbar[n]=0;


			for(k=0;k<n3;k++)
				for(j=0;j<n2;j++)
					for(i=0;i<n1;i++){
						for(n=0;n<6;n++)
							stressbar[n]+=stress[(k*n2*(n1)+j*(n1)+i)*6+n]*volumeVoxel;
			}

			// Write stressbar and strainbar for each iteration to output file

			change_basis(stressbar,stressbar33,aux66,aux3333,1);
			stressref=stressbar33[ictrl1][ictrl2];
			cout<<"STRESS FIELD ERROR:"<<errstress/stressref<<endl;
			cout<<"STRAIN FIELD ERROR:"<<errstrain/strainref<<endl;

			errorOut<<iteration<<" ";
			errorOut<<errstress/stressref<<" ";
			errorOut<<errstrain/strainref<<" ";
			errorOut<<straingradrate33[0][0]<<" ";
			errorOut<<straingradrate33[1][1]<<" ";
			errorOut<<straingradrate33[2][2]<<" ";
			errorOut<<straingradrate33[1][2]<<" ";
			errorOut<<straingradrate33[0][2]<<" ";
			errorOut<<straingradrate33[0][1]<<" ";
			errorOut<<stressbar33[0][0]<<" ";
			errorOut<<stressbar33[1][1]<<" ";
			errorOut<<stressbar33[2][2]<<" ";
			errorOut<<stressbar33[1][2]<<" ";
			errorOut<<stressbar33[0][2]<<" ";
			errorOut<<stressbar33[0][1]<<" ";

			errorOut<<endl;

		} // Iteration loop

	} // Step loop

	// Write field output
	for(k=0;k<n3;k++)
		for(j=0;j<n2;j++)
			for(i=0;i<n1;i++){
				
				for(m=0;m<6;m++){
					strain[m]=strainbar[m]+straintilde[(k*n2*(n1)+j*(n1)+i)*6+m];
					stress6[m]=stress[(k*n2*(n1)+j*(n1)+i)*6+m];
				}

				change_basis(strain,strainout,aux66,aux3333,1);
				change_basis(stress6,stressout,aux66,aux3333,1);

				fieldsOut<<i+1<<" ";
				fieldsOut<<j+1<<" ";
				fieldsOut<<k+1<<" ";
				fieldsOut<<grainID[k*n2*(n1)+j*(n1)+i]<<" ";

				fieldsOut<<strainout[0][0]<<" ";
				fieldsOut<<strainout[1][1]<<" ";
				fieldsOut<<strainout[2][2]<<" ";
				fieldsOut<<strainout[1][2]<<" ";
				fieldsOut<<strainout[0][2]<<" ";
				fieldsOut<<strainout[0][1]<<" ";

				fieldsOut<<stressout[0][0]<<" ";
				fieldsOut<<stressout[1][1]<<" ";
				fieldsOut<<stressout[2][2]<<" ";
				fieldsOut<<stressout[1][2]<<" ";
				fieldsOut<<stressout[0][2]<<" ";
				fieldsOut<<stressout[0][1]<<" ";

				fieldsOut<<endl;
	}

	fftw_free(out);
	return 0;
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
			
				multiply3333x33((double *) G,Cref,fourierTensor,2,4);
	
				findInverse((double *)G,det,3);

				for(p=0;p<3;p++)
					for(q=0;q<3;q++)
						for(r=0;r<3;r++)
							for(s=0;s<3;s++)
								gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i].tensor[p][q][r][s]=-1*G[p][r]*fourierTensor[q][s]; 
            }
        }
    }

			    for(p=0;p<3;p++)
					for(q=0;q<3;q++)
						for(r=0;r<3;r++)
							for(s=0;s<3;s++)
								gammaHat[0].tensor[p][q][r][s]=0;


}