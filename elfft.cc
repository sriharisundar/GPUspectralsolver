#define __CL_ENABLE_EXCEPTIONS
#include "readInput.h"
#include "globalVariables.h"
#include "solverFunctions.h"
#include "matrixOperations.h"
#include "printFunctions.h"
#include <iostream>
#include <fstream>
#include <fftw3.h>
#include "clFFT.h"
#include <CL/cl.hpp>
using namespace std;

#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_DEFAULT
#endif

// From Fabian Dournac
// https://dournac.org/info/fft_gpu
int FFT_3D_OpenCL(float *tab[], const char* direction, int sizex, int sizey, int sizez)
{
	int i;

	/* OpenCL variables. */
	cl_int err;
	cl_platform_id platform = 0;
	cl_device_id device = 0;
	cl_context ctx = 0;
	cl_command_queue queue = 0;

	/* Input and Output  buffer. */
	cl_mem buffersIn[2]  = {0, 0};
	cl_mem buffersOut[2] = {0, 0};

	/* Temporary buffer. */
	cl_mem tmpBuffer = 0;

	/* Size of temp buffer. */
	size_t tmpBufferSize = 0;
	int status = 0;
	int ret = 0;

	/* Total size of FFT. */
	size_t N = sizex*sizey*sizez;

	/* FFT library realted declarations. */
	clfftPlanHandle planHandle;
	clfftDim dim = CLFFT_3D;
	size_t clLengths[3] = {sizex, sizey, sizez};

	/* Setup OpenCL environment. */
	err = clGetPlatformIDs(1, &platform, NULL);
	err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);

	/* Create an OpenCL context. */
	ctx = clCreateContext(NULL, 1, &device, NULL, NULL, &err);

	/* Create a command queue. */
	queue = clCreateCommandQueue(ctx, device, 0, &err);

	/* Setup clFFT. */
	clfftSetupData fftSetup;
	err = clfftInitSetupData(&fftSetup);
	err = clfftSetup(&fftSetup);

	/* Create a default plan for a complex FFT. */
	err = clfftCreateDefaultPlan(&planHandle, ctx, dim, clLengths);

	/* Set plan parameters. */
	err = clfftSetPlanPrecision(planHandle, CLFFT_SINGLE);
	err = clfftSetLayout(planHandle, CLFFT_COMPLEX_PLANAR, CLFFT_COMPLEX_PLANAR);
	err = clfftSetResultLocation(planHandle, CLFFT_OUTOFPLACE);

	/* Bake the plan. */
	err = clfftBakePlan(planHandle, 1, &queue, NULL, NULL);

	/* Real and Imaginary arrays. */
	cl_float* inReal  = (cl_float*) malloc (N * sizeof (cl_float));
	cl_float* inImag  = (cl_float*) malloc (N * sizeof (cl_float));
	cl_float* outReal = (cl_float*) malloc (N * sizeof (cl_float));
	cl_float* outImag = (cl_float*) malloc (N * sizeof (cl_float));

	/* Initialization of inReal, inImag, outReal and outImag. */
        for(i=0; i<N; i++)
	{
		inReal[i]  = tab[0][i];
		inImag[i]  = 0.0f;
		outReal[i] = 0.0f;
		outImag[i] = 0.0f;
	}

	/* Create temporary buffer. */
	status = clfftGetTmpBufSize(planHandle, &tmpBufferSize);

	if ((status == 0) && (tmpBufferSize > 0)) {
		tmpBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, tmpBufferSize, 0, &err);
		if (err != CL_SUCCESS)
			printf("Error with tmpBuffer clCreateBuffer\n");
	}

	/* Prepare OpenCL memory objects : create buffer for input. */
	buffersIn[0] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
			N * sizeof(cl_float), inReal, &err);
	if (err != CL_SUCCESS)
		printf("Error with buffersIn[0] clCreateBuffer\n");

	/* Enqueue write tab array into buffersIn[0]. */
	err = clEnqueueWriteBuffer(queue, buffersIn[0], CL_TRUE, 0, N *
			sizeof(float),
			inReal, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		printf("Error with buffersIn[0] clEnqueueWriteBuffer\n");

	/* Prepare OpenCL memory objects : create buffer for input. */
	buffersIn[1] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
			N * sizeof(cl_float), inImag, &err);
	if (err != CL_SUCCESS)
		printf("Error with buffersIn[1] clCreateBuffer\n");

	/* Enqueue write tab array into buffersIn[0]. */
	err = clEnqueueWriteBuffer(queue, buffersIn[1], CL_TRUE, 0, N * sizeof(float),
			inImag, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		printf("Error with buffersIn[1] clEnqueueWriteBuffer\n");

	/* Prepare OpenCL memory objects : create buffer for output. */
	buffersOut[0] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N *
			sizeof(cl_float), outReal, &err);
	if (err != CL_SUCCESS)
		printf("Error with buffersOut[0] clCreateBuffer\n");

	buffersOut[1] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N *
			sizeof(cl_float), outImag, &err);
	if (err != CL_SUCCESS)
		printf("Error with buffersOut[1] clCreateBuffer\n");

	/* Execute Forward or Backward FFT. */
	if(strcmp(direction,"forward") == 0)
	{
		/* Execute the plan. */
		err = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &queue, 0, NULL, NULL,
				buffersIn, buffersOut, tmpBuffer);
	}
	else if(strcmp(direction,"backward") == 0)
	{
		/* Execute the plan. */
		err = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &queue, 0, NULL, NULL,
				buffersIn, buffersOut, tmpBuffer);
	}

	/* Wait for calculations to be finished. */
	err = clFinish(queue);

	/* Fetch results of calculations : Real and Imaginary. */
	err = clEnqueueReadBuffer(queue, buffersOut[0], CL_TRUE, 0, N * sizeof(float), tab[0],
			0, NULL, NULL);
	err = clEnqueueReadBuffer(queue, buffersOut[1], CL_TRUE, 0, N * sizeof(float), tab[1],
			0, NULL, NULL);

	/* Release OpenCL memory objects. */
	clReleaseMemObject(buffersIn[0]);
	clReleaseMemObject(buffersIn[1]);
	clReleaseMemObject(buffersOut[0]);
	clReleaseMemObject(buffersOut[1]);
	clReleaseMemObject(tmpBuffer);

	/* Release the plan. */
	err = clfftDestroyPlan(&planHandle);

	/* Release clFFT library. */
	clfftTeardown();

	/* Release OpenCL working objects. */
	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	return ret;
}

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
	int fftchoice;

	fftw_complex *in;
	fftw_complex *out;
	fftw_plan plan_backward;
	fftw_plan plan_forward;



	if (argc!=2){
		cout<<"Pass input file name as argument"<<endl;
		return 0;
	}

	readinput(argv[1]);

	if (argc!=3){
		cout<<"Enter 1 for FFTW or 2 for clFFT";
		cin>>fftchoice;
	}
	else
		fftchoice=atoi(argv[2]);

	if (fftchoice==1){
		in=(fftw_complex *) *fftw_alloc_complex(N3*N2*(N1));
		out=(fftw_complex *) *fftw_alloc_complex(N3*N2*(N1));

		plan_forward=fftw_plan_dft_3d(N3,N2,N1,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
		plan_backward=fftw_plan_dft_3d ( N3, N2, N1, out, in, FFTW_BACKWARD,FFTW_ESTIMATE );
	}

	else { 
		
	;}

	cout<<"Output file:"<<outputFile<<endl;
	fstream fieldsOut;
	fieldsOut.open(outputFile,ios::out);

	fstream errorOut;
	errorOut.open("err.out",ios::out);

	// Initialize stress and strain fields
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
						for(i=0;i<N1;i++){
							work[k][j][i][n]=out[k*N2*N1+j*N1+i][0];
							workim[k][j][i][n]=out[k*N2*N1+j*N1+i][1];
				}

			}

			// Convert stress to tensorial form
			// Multiply with gamma operator
			cout<<"Gamma convolution"<<endl<<endl;
			for(k=0;k<N3;k++)
				for(j=0;j<N2;j++)
					for(i=0;i<(N1);i++){
						change_basis(work[k][j][i],work33,aux66,aux3333,1);
						change_basis(workim[k][j][i],work33im,aux66,aux3333,1);
						
						multiply3333x33(ddefgrad[k][j][i],gammaHat[k*n2*(n1)+j*(n1)+i],work33,3,4);
						multiply3333x33(ddefgradim[k][j][i],gammaHat[k*n2*(n1)+j*(n1)+i],work33im,3,4);
			}
			
			// Arrange data for out
			cout<<"Inverse FFT to get deformation gradient"<<endl<<endl;
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
							for(i=0;i<N1;i++)
								ddefgrad[k][j][i][m][n]=in[k*N2*(N1)+j*(N1)+i][0]/prodDim;

			}

			// Get symmetric part of defgrad
		
			for(k=0;k<N3;k++)
				for(j=0;j<N2;j++)
					for(i=0;i<N1;i++){
						symmetric(ddefgrad[k][j][i],aux33);
						change_basis(straintilde[k][j][i],aux33,aux66,aux3333,2);
			}


			cout<<"Augmented Lagrangian method for stress update"<<endl<<endl;
			augmentLagrangian();
			
			for(n=0;n<6;n++)
				stressbar[n]=0;


			for(k=0;k<N3;k++)
				for(j=0;j<N2;j++)
					for(i=0;i<N1;i++){
						for(n=0;n<6;n++)
							stressbar[n]+=stress[k][j][i][n]*volumeVoxel;
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
	for(k=0;k<N3;k++)
		for(j=0;j<N2;j++)
			for(i=0;i<N1;i++){
				
				for(m=0;m<6;m++){
					strain[m]=strainbar[m]+straintilde[k][j][i][m];
					stress6[m]=stress[k][j][i][m];
				}

				change_basis(strain,strainout,aux66,aux3333,1);
				change_basis(stress6,stressout,aux66,aux3333,1);

				fieldsOut<<i+1<<" ";
				fieldsOut<<j+1<<" ";
				fieldsOut<<k+1<<" ";
				fieldsOut<<grainID[k][j][i]<<" ";

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