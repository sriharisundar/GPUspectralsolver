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
	double prodDim=n1*n2*n3;
	double volumeVoxel=1.0/prodDim;
	int iteration,step;

	fftw_complex *in;
	fftw_complex *out;
	fftw_plan plan_backward;
	fftw_plan plan_forward;

	if (argc!=2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	initglobal();
	readinput(argv[1]);

	in=(fftw_complex *) *fftw_alloc_complex(N3*N2*(N1));
	out=(fftw_complex *) *fftw_alloc_complex(N3*N2*(N1));

	plan_forward=fftw_plan_dft_3d(N3,N2,N1,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	plan_backward=fftw_plan_dft_3d ( N3, N2, N1, out, in, FFTW_BACKWARD,FFTW_ESTIMATE );	

	cout<<"Output file:"<<outputFile<<endl;
	fstream fieldsOut;
	fieldsOut.open(outputFile,ios::out);

	fstream errorOut;
	errorOut.open("err.out",ios::out);

	//Initialize stress and strain fields
	// sg - stress
	// dtilde-straingradient - straintilde
	
	for(k=0;k<N3;k++)
		for(j=0;j<N2;j++)
			for(i=0;i<N1;i++){

				for(n=0;n<6;n++)
					straintilde[k][j][i][n]=0;

				for(n=0;n<6;n++){
					if(phaseID[k][j][i]==2)
						stress[k][j][i][n]=0;
					else{
						stress[k][j][i][n]=0;
						for(m=0;m<6;m++)
							stress[k][j][i][n]+=cloc[k][j][i][n][m]*(strainbar[m]);
					}
					stressbar[n]+=stress[k][j][i][n]*volumeVoxel;
				}
	}
//
//			for(i=0;i<N1;i++)
//		for(j=0;j<N2;j++)
//	for(k=0;k<N3;k++){
//					cout<<"GrainID:"<<grainID[k][j][i]<<" "<<endl;
////					print2darray(cloc[k][j][i]);
//					print1darray(stress[k][j][i],6);}
//
//	print1darray((double *) strainbar,6);

	change_basis(stressbar,stressbar33,aux66,aux3333,1);
	stressref=stressbar33[ictrl1][ictrl2];

	for(step=1;step<=nsteps;step++){
		iteration=0;

		err2mod=2*error;

		findGammaHat(C0_3333);

		while(iteration<itermax && err2mod > error){
			iteration++;
			//debugssh cout<<"--------------------------------------------------------------"<<endl;
			//debugssh cout<<"ITERATION:"<<iteration<<endl;

			//arrange data for in
			//perform forward FFT
			//debugssh cout<<"Forward FFT of polarization field"<<endl<<endl;
			for(n=0;n<6;n++){
				
				for(k=0;k<N3;k++)
					for(j=0;j<N2;j++)
						for(i=0;i<N1;i++){
							in[k*N2*N1+j*N1+i][0]=stress[k][j][i][n];
							in[k*N2*N1+j*N1+i][1]=0;
							for(m=0;m<6;m++)
								in[k*N2*N1+j*N1+i][0]-=C0_66[n][m]*straintilde[k][j][i][m];
						}

				fftw_execute(plan_forward);
				

				for(k=0;k<N3;k++)
					for(j=0;j<N2;j++)
						for(i=0;i<(N1);i++){
							work[k][j][i][n]=out[k*N2*(N1)+j*(N1)+i][0];
							workim[k][j][i][n]=out[k*N2*(N1)+j*(N1)+i][1];
							//print1darray((double *)work[k][j][i],6);
				}
			}

			//convert stress to tensorial form
			//multiply with gamma operator
			//debugssh cout<<"Gamma convolution"<<endl<<endl;
			for(k=0;k<N3;k++)
				for(j=0;j<N2;j++)
					for(i=0;i<(N1);i++){

						change_basis(work[k][j][i],work33,aux66,aux3333,1);
						change_basis(workim[k][j][i],work33im,aux66,aux3333,1);

////						cout<<"Fourier point:"<<i<<" "<<j<<" "<<k<<" "<<endl;
//						print2darray(work33);
//						print2darray(work33im);

						multiply3333x33(ddefgrad[k][j][i],gammaHat[k*n2*(n1)+j*(n1)+i],work33,3,4);
						multiply3333x33(ddefgradim[k][j][i],gammaHat[k*n2*(n1)+j*(n1)+i],work33im,3,4);
//						print2darray(ddefgrad[k][j][i]);
//						print2darray(ddefgradim[k][j][i]);
			}
			
			//arrange data for out
			//debugssh cout<<"Inverse FFT to get deformation gradient"<<endl<<endl;
			for(m=0;m<3;m++)
				for(n=0;n<3;n++){

					for(k=0;k<N3;k++)
						for(j=0;j<N2;j++)
							for(i=0;i<(N1);i++){
								out[k*N2*(N1)+j*(N1)+i][0]=ddefgrad[k][j][i][m][n];
								out[k*N2*(N1)+j*(N1)+i][1]=ddefgradim[k][j][i][m][n];
					}

					fftw_execute(plan_backward);

					for(k=0;k<N3;k++)
						for(j=0;j<N2;j++)
							for(i=0;i<N1;i++){
////								cout<<delta[k][j][i]<<endl;
								ddefgrad[k][j][i][m][n]=in[k*N2*(N1)+j*(N1)+i][0]/prodDim;
							}

			}

			//get symmetric part of defgrad
		
			for(k=0;k<N3;k++)
				for(j=0;j<N2;j++)
					for(i=0;i<N1;i++){
//						print2darray(ddefgrad[k][j][i]);
						symmetric(ddefgrad[k][j][i],aux33);
						change_basis(straintilde[k][j][i],aux33,aux66,aux3333,2);
			}


			//debugssh cout<<"Augmented Lagrangian method for stress update"<<endl<<endl;
			augmentLagrangian();
			for(n=0;n<6;n++)
				stressbar[n]=0;


			for(k=0;k<N3;k++)
				for(j=0;j<N2;j++)
					for(i=0;i<N1;i++)
						for(n=0;n<6;n++)
							stressbar[n]+=stress[k][j][i][n]*volumeVoxel;

			change_basis(stressbar,stressbar33,aux66,aux3333,1);
			stressref=stressbar33[ictrl1][ictrl2];
			//print1darray((double *)stressbar,6);
			//print2darray(C0_66);

            //debugssh cout<<"STRESS FIELD ERROR:"<<errstress/stressref<<endl;
			//debugssh cout<<"STRAIN FIELD ERROR:"<<errstrain/strainref<<endl;
		}
	}

	for(i=0;i<N1;i++)
		for(j=0;j<N2;j++)
			for(k=0;k<N3;k++){
				
				for(m=0;m<6;m++){
					strain[m]=strainbar[m]+straintilde[k][j][i][m];
					stress6[m]=stress[k][j][i][m];
				}

				change_basis(strain,strainout,aux66,aux3333,1);
				change_basis(stress6,stressout,aux66,aux3333,1);

				fieldsOut<<strainout[0][0]<<" ";
				fieldsOut<<strainout[1][1]<<" ";
				fieldsOut<<strainout[2][2]<<" ";
				fieldsOut<<strainout[2][3]<<" ";
				fieldsOut<<strainout[1][3]<<" ";
				fieldsOut<<strainout[1][2]<<" ";

				fieldsOut<<stressout[0][0]<<" ";
				fieldsOut<<stressout[1][1]<<" ";
				fieldsOut<<stressout[2][2]<<" ";
				fieldsOut<<stressout[2][3]<<" ";
				fieldsOut<<stressout[1][3]<<" ";
				fieldsOut<<stressout[1][2]<<" ";

				fieldsOut<<endl;
	}

	fftw_free(out);
	return 0;
}