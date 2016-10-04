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
	double stressbar33[3][3],aux66[6][6],aux3333[3][3][3][3],work33[3][3],work33im[3][3];
	double err2mod;
	int iteration,step;

	fftw_complex *out;
	fftw_plan plan_backward;
	fftw_plan plan_forward;

	if (argc!=2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	initglobal();
	readinput(argv[1]);

	//out = new fftw_complex [N1*N2*(N3/2+1)];
	out=(fftw_complex *) *fftw_alloc_complex(N3*N2*(N1/2+1));

	plan_forward=fftw_plan_dft_r2c_3d(N3,N2,N1,(double*) delta,out,FFTW_ESTIMATE);
	plan_backward=fftw_plan_dft_c2r_3d ( N3, N2, N1, out, (double*) delta, FFTW_ESTIMATE );	

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

						for(m=0;m<6;m++)
							stress[k][j][i][n]=cloc[k][j][i][n][m]*straintilde[k][j][i][m];
					}
					stressbar[n]+=stress[k][j][i][n]*wgt;

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
			//cout<<"ITER:"<<iteration<<endl;

			for(n=0;n<6;n++){
				delta[k][j][i]=stress[k][j][i][n];
				for(k=0;k<N3;k++)
					for(j=0;j<N2;j++)
						for(i=0;i<N1;i++)
							for(m=0;m<6;m++)
								delta[k][j][i]-=C0_66[n][m]*straintilde[k][j][i][m];
				
				fftw_execute(plan_forward);
				
				for(k=0;k<N3;k++)
					for(j=0;j<N2;j++)
						for(i=0;i<(N1/2+1);i++){
							work[k][j][i][n]=out[k*N2*(N1/2+1)+j*(N1/2+1)+i][0];
							workim[k][j][i][n]=out[k*N2*(N1/2+1)+j*(N1/2+1)+i][1];
						}
			}

			for(k=0;k<N3;k++)
				for(j=0;j<N2;j++)
					for(i=0;i<(N1/2+1);i++){

						change_basis(work[k][j][i],work33,aux66,aux3333,1);
						change_basis(workim[k][j][i],work33im,aux66,aux3333,1);

						

					}
		}
	}

	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	fftw_free(out);
	augmentLagrangian();
	return 0;
}