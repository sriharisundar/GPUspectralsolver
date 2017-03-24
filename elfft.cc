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

namespace util {

inline std::string loadProgram(std::string input)
{
    std::ifstream stream(input.c_str());
    if (!stream.is_open()) {
        std::cout << "Cannot open file: " << input << std::endl;
        exit(1);
    }

     return std::string(
        std::istreambuf_iterator<char>(stream), 
        (std::istreambuf_iterator<char>()));
}

}

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
    int fftchoice;

//fftw variables
    double *delta;
    fftw_complex *out;
    fftw_plan plan_backward, plan_forward;

//clFFT variables
    cl_int err;
    cl::Buffer d_stress;
    cl::Buffer d_straintilde;
    cl::Buffer d_gammaHat;
    cl::Buffer d_C0_66;

    cl::Context context(DEVICE);        

    cl::CommandQueue queue(context);

    cl::Program program(context,util::loadProgram("kernelFunctions.cl"),true);

    cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,int,int> findAuxiliaryStress(program,"findAuxiliaryStress");

    clfftPlanHandle planHandleFWD,planHandleBWD;
    clfftDim dim = CLFFT_3D;
    clfftSetupData fftSetup;
    size_t clLengths[3];
    size_t clStridesFWD[3],clStridesBWD[3];

    if (argc<2){
        cout<<"Pass input file name as argument"<<endl;
        return 0;
    }

    //if (argc<3){
    //    cout<<"Enter 1 for FFTW or 2 for clFFT:";
    //    cin>>fftchoice;
    //}
    //else
    //    fftchoice=atoi(argv[2]);

    readinput(argv[1]);

    prodDim=n1*n2*n3;
    volumeVoxel=1.0/prodDim;

    delta=new double[n3*n2*n1];
    out=(fftw_complex *) *fftw_alloc_complex(n3*n2*(n1/2+1));

    plan_forward=fftw_plan_dft_r2c_3d(n3, n2, n1, delta, out, FFTW_ESTIMATE);
    plan_backward=fftw_plan_dft_c2r_3d(n3, n2, n1, out, delta, FFTW_ESTIMATE);

    clLengths[0]=n1;
    clLengths[1]=n2;
    clLengths[2]=n3;

    clStridesFWD[0]=6;
    clStridesFWD[1]=6*n1;
    clStridesFWD[2]=6*n1*n2;

    clStridesBWD[0]=9;
    clStridesBWD[1]=9*n1;
    clStridesBWD[2]=9*n1*n2;

    d_stress=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double)*6*prodDim);
    d_straintilde=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double)*6*prodDim);
    d_gammaHat=cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(fourthOrderTensor)*prodDim);
    d_C0_66=cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double)*36);

    // Setup clFFT
    err=clfftInitSetupData(&fftSetup);
    err=clfftSetup(&fftSetup);

    // Create a default plan for a complex FFT. This will be modified further to satisfy our conditions
    err=clfftCreateDefaultPlan(&planHandleFWD, context(), dim, clLengths);

    // Set plan parameters
    err=clfftSetPlanPrecision(planHandleFWD, CLFFT_DOUBLE);
    err=clfftSetResultLocation(planHandleFWD, CLFFT_OUTOFPLACE);
    
    // Batching is required as we will be executing the plan for all the six components of the auxiliary stress tensor together
    err=clfftSetPlanBatchSize(planHandleFWD, 6);
    
    // Distance specifies the separation between the starting elements of the arrays corresponding to each batch for input & output arrays
    err=clfftSetPlanDistance(planHandleFWD, 1, 1);
    
    /* Due to batching of the transform, the separation between consecutive elements for one particular batch, along each dimension
       is not the default values anymore. The x stride becomes 6, y stride becomes 6*lenX, z stride becomes 6*lenX*lenY. This is 
       specified by the clStrides variable defined previously. */
    err=clfftSetPlanInStride(planHandleFWD, dim, clStridesFWD);
    err=clfftSetPlanOutStride(planHandleFWD, dim, clStridesFWD);

    /* Since we will be using the Hermitian symmetry which is applicable when the input data for the forward transform is all real,
       we choose the plan input and output accordingly. For real transforms clFFT selects the direction of transform based on type
       of the input and output buffers. Thus we create two plan handles one for forward transform and one for reverse transform. */
    err=clfftCopyPlan(&planHandleBWD,context(), planHandleFWD);

    // Specify the type of input and output buffers for forward and reverse transforms
    err=clfftSetLayout(planHandleFWD, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
    err=clfftSetLayout(planHandleBWD, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL);
    
    err=clfftSetPlanBatchSize(planHandleFWD, 9);

    err=clfftSetPlanInStride(planHandleBWD, dim, clStridesBWD);
    err=clfftSetPlanOutStride(planHandleBWD, dim, clStridesBWD);


    // Bake the plan
    err=clfftBakePlan(planHandleFWD, 1, &queue(), NULL, NULL);
    err=clfftBakePlan(planHandleBWD, 1, &queue(), NULL, NULL);

    cout<<"Output file:"<<outputFile<<endl;
    fstream fieldsOut;
    fieldsOut.open(outputFile, ios::out);

    fstream errorOut;
    errorOut.open("err.out", ios::out);

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

    change_basis(stressbar, stressbar33, aux66, aux3333, 1);
    stressref=stressbar33[ictrl1][ictrl2];

    for(step=1;step<=nsteps;step++){
        iteration=0;

        err2mod=2*error;

        findGammaHat(C0_3333);

        err = queue.enqueueWriteBuffer(d_C0_66, CL_TRUE, 0, 36*sizeof(double), C0_66);

        err = queue.enqueueWriteBuffer(d_gammaHat, CL_TRUE, 0, prodDim*sizeof(fourthOrderTensor), 
                                    gammaHat);

        err = queue.enqueueWriteBuffer(d_stress, CL_TRUE, 0, prodDim*6*sizeof(double), 
                                    stress);

        err = queue.enqueueWriteBuffer(d_straintilde, CL_TRUE, 0, prodDim*6*sizeof(double), 
                                    straintilde);

        while(iteration<itermax && err2mod > error){
            iteration++;
            cout<<"--------------------------------------------------------------"<<endl;
            cout<<"ITERATION:"<<iteration<<endl;

            // Arrange data for in
            // Perform forward FFT
            cout<<"Forward FFT of polarization field"<<endl<<endl;

            for(n=0;n<6;n++){
                
            err=clFinish(queue());

//                for(k=0;k<n3;k++)
//                    for(j=0;j<n2;j++)
//                        for(i=0;i<n1;i++){
//                            delta[(k*n2*(n1)+j*(n1)+i)]=stress[(k*n2*(n1)+j*(n1)+i)*6+n];
//                            for(m=0;m<6;m++)
//                                delta[(k*n2*(n1)+j*(n1)+i)]-=C0_66[n][m]*straintilde[(k*n2*(n1)+j*(n1)+i)*6+m];
//                }
//
//                fftw_execute(plan_forward);
//
//                for(k=0;k<n3;k++)
//                    for(j=0;j<n2;j++)
//                        for(i=0;i<(n1/2+1);i++){
//                            work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][0];
//                            workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][1];                            
                //}

            }

            // Convert stress to tensorial form
            // Multiply with gamma operator
            cout<<"Gamma convolution"<<endl<<endl;
            for(k=0;k<n3;k++)
                for(j=0;j<n2;j++)
                    for(i=0;i<(n1/2+1);i++){
                        change_basis(&work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6], work33, aux66, aux3333, 1);
                        change_basis(&workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6], work33im, aux66, aux3333, 1);

                        multiply3333x33(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9], gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i], work33, 3, 4);
                        multiply3333x33(&ddefgradim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9], gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i], work33im, 3, 4);
            }
            
            //change_basis(strainbar, strainbar33, aux66, aux3333, 1);

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

                    //out[0][0]=strainbar33[m][n];

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
                        symmetric(&ddefgrad[(k*n2*(n1)+j*(n1)+i)*9], (double *) aux33);
                        change_basis(&straintilde[(k*n2*(n1)+j*(n1)+i)*6], aux33, aux66, aux3333, 2);
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
            change_basis(stressbar, stressbar33, aux66, aux3333, 1);
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

    clfftTeardown();

    // Write field output
    for(k=0;k<n3;k++)
        for(j=0;j<n2;j++)
            for(i=0;i<n1;i++){
                
                for(m=0;m<6;m++){
                    strain[m]=strainbar[m]+straintilde[(k*n2*(n1)+j*(n1)+i)*6+m];
                    stress6[m]=stress[(k*n2*(n1)+j*(n1)+i)*6+m];
                }

                change_basis(strain, strainout, aux66, aux3333, 1);
                change_basis(stress6, stressout, aux66, aux3333, 1);

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

//    err=clfftCreateDefaultPlan(&planHandleBWD, context(), dim, clLengths);
//    err=clfftSetPlanPrecision(planHandleBWD, CLFFT_DOUBLE);
//    err=clfftSetResultLocation(planHandleBWD, CLFFT_OUTOFPLACE);
//    err=clfftSetPlanBatchSize(planHandleBWD, 6);
//    err=clfftSetPlanDistance(planHandleBWD, 1, 1);
//    err=clfftSetPlanInStride(planHandleBWD, dim, clStrides);
//    err=clfftSetPlanOutStride(planHandleBWD, dim, clStrides);