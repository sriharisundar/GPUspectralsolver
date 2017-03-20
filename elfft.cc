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
	double stressbar33[3][3],work33[3][3],work33im[3][3];
	double aux33[3][3],aux66[6][6],aux3333[3][3][3][3];
	double strain[6],stress6[6];
	double strainout[3][3],stressout[3][3];
	double err2mod;
	double prodDim;
	double volumeVoxel;
	int iteration,step;
	double *delta;
	double *deltaout;

	int N1,N2,N3; 
	fftw_complex *out;
	fftw_plan plan_backward;
	fftw_plan plan_forward;


	if (argc<2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	readinput(argv[1]);

	N1=n1;
	N2=n2;
	N3=n3;

	prodDim=n1*n2*n3;
	volumeVoxel=1.0/prodDim;

	delta=new double[n3*n2*n1];
	deltaout=new double[n3*n2*n1];
	out=(fftw_complex *) *fftw_alloc_complex(n3*n2*(n1/2+1));

	plan_forward=fftw_plan_dft_r2c_3d(N3,N2,N1,delta,out,FFTW_ESTIMATE);
	plan_backward=fftw_plan_dft_c2r_3d ( N3, N2, N1, out, deltaout, FFTW_ESTIMATE );	

	cout<<"Output file:"<<outputFile<<endl;
	fstream fieldsOut;
	fieldsOut.open(outputFile,ios::out);

	fstream errorOut;
	errorOut.open("err.out",ios::out);

	// Initialize stress and strain fields
	// sg - stress
	// dtilde-straingradient - straintilde
	//for(k=0;k<n3;k++)
		//for(j=0;j<n2;j++)
			//for(i=0;i<n1;i++){
				//cout<<k<<" "<<j<<" "<<i<<endl;
				//print2darray(&cloc[(k*n2*(n1)+j*(n1)+i)*36],6);
			//}
		
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
			//debughermitkimaakacout<<"--------------------------------------------------------------"<<endl;
			//debughermitkimaakacout<<"ITERATION:"<<iteration<<endl;

			// Arrange data for in
			// Perform forward FFT
			//debughermitkimaakacout<<"Forward FFT of polarization field"<<endl<<endl;

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
							//if(iteration == 1)
								//cout<<i<<" "<<j<<" "<<k<<" "<<out[k*n2*(n1/2+1)+j*(n1/2+1)+i][0]<<endl;
							work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][0];
							workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][1];							
				}

			}

			// Convert stress to tensorial form
			// Multiply with gamma operator
			//debughermitkimaakacout<<"Gamma convolution"<<endl<<endl;
			for(k=0;k<n3;k++)
				for(j=0;j<n2;j++)
					for(i=0;i<(n1/2+1);i++){
						change_basis(&work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6],work33,aux66,aux3333,1);
						change_basis(&workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6],work33im,aux66,aux3333,1);

						//if(iteration == 1){
							//cout<<i<<" "<<j<<" "<<k<<" "<<endl;
							//print2darray((double*)work33,3);}

						multiply3333x33(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i],work33,3,4);
						multiply3333x33(&ddefgradim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i],work33im,3,4);

						//if(iteration == 1){
							//print2darray(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],3);
							//print4darray(gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i].tensor);}


			}
			
			// Arrange data for out
			//debughermitkimaakacout<<"Inverse FFT to get deformation gradient"<<endl<<endl;
			for(m=0;m<3;m++)
				for(n=0;n<3;n++){

								if(iteration == 1)
									cout<<m<<" "<<n<<endl;
					for(k=0;k<N3;k++)
						for(j=0;j<N2;j++)
							for(i=0;i<(N1/2+1);i++){
								out[k*N2*(N1/2+1)+j*(N1/2+1)+i][0]=ddefgrad[(k*n2*(n1)+j*(n1)+i)*9+3*m+n];
								out[k*N2*(N1/2+1)+j*(N1/2+1)+i][1]=ddefgradim[(k*n2*(n1)+j*(n1)+i)*9+3*m+n];
					}

					fftw_execute(plan_backward);

					for(k=0;k<N3;k++)
						for(j=0;j<N2;j++)
							for(i=0;i<N1;i++){
								if(iteration == 1)
									cout<<i<<" "<<j<<" "<<k<<" "<<deltaout[(k*n2*(n1)+j*(n1)+i)]<<endl;
								ddefgrad[(k*n2*(n1)+j*(n1)+i)*9+3*m+n]=deltaout[(k*n2*(n1)+j*(n1)+i)]/prodDim;
							}

			}

			// Get symmetric part of defgrad
		
			for(k=0;k<n3;k++)
				for(j=0;j<n2;j++)
					for(i=0;i<n1;i++){
						symmetric(&ddefgrad[(k*n2*(n1)+j*(n1)+i)*9],(double *) aux33);
						change_basis(&straintilde[(k*n2*(n1)+j*(n1)+i)*6],aux33,aux66,aux3333,2);
			}


			//debughermitkimaakacout<<"Augmented Lagrangian method for stress update"<<endl<<endl;
//                                        if(iteration==1)
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
			//debughermitkimaakacout<<"STRESS FIELD ERROR:"<<errstress/stressref<<endl;
			//debughermitkimaakacout<<"STRAIN FIELD ERROR:"<<errstrain/strainref<<endl;

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