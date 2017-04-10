#define __CL_ENABLE_EXCEPTIONS
#include "readInput.h"
#include "globalVariables.h"
#include "solverFunctions.h"
#include "matrixOperations.h"
#include "printFunctions.h"
#include <iostream>
#include <iomanip>
#include <fstream>
//#include <fftw3.h>
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

unsigned getDeviceList(std::vector<cl::Device>& devices)
{
  cl_int err;

  // Get list of platforms
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);

  // Enumerate devices
  for (int i = 0; i < platforms.size(); i++)
  {
    cl_uint num = 0;
    std::vector<cl::Device> plat_devices;
    platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &plat_devices);
    devices.insert(devices.end(), plat_devices.begin(), plat_devices.end());
  }

  return devices.size();
}

int main(int argc, char *argv[])
{    
    int i,j,k,n,m;
    double stressbar33[3][3],work33[3][3],work33im[3][3];
    double aux33[3][3],aux66[6][6],aux3333[3][3][3][3];
    double strain[6],stress6[6];
    double strainout[3][3],stressout[3][3];
    double err2mod;
    long prodDim,prodDimHermitian;
    double volumeVoxel;
    int iteration,step;
    int fftchoice;

    double* stressaux;
    double* errors;

    vector6_complex* stressFourier;
    tensor33_complex* ddefgradFourier;

////fftw variables
    double *delta;
//    fftw_complex *out;
//    fftw_plan plan_backward, plan_forward;

//clFFT variables
    cl_int err;
    cl::Buffer d_stress;
    cl::Buffer d_straintilde;
    cl::Buffer d_C0_66;
    cl::Buffer d_gammaHat;
    cl::Buffer d_strainbar;
    cl::Buffer d_fsloc;
    cl::Buffer d_stressFourier;
    cl::Buffer d_ddefgradFourier;
    cl::Buffer d_ddefgrad;
    cl::Buffer d_errors;

    cl_uint deviceIndex = 0;

    // Get list of devices
    std::vector<cl::Device> devices;
    unsigned numDevices = getDeviceList(devices);

    // Check device index in range
    if (deviceIndex >= numDevices)
    {
      std::cout << "Invalid device index (try '--list')\n";
      return EXIT_FAILURE;
    }

    cl::Device device = devices[0];

    std::vector<cl::Device> chosen_device;
    chosen_device.push_back(device);
    cl::Context context(chosen_device);
    cl::CommandQueue queue(context, device);

    //cl::Context context(DEVICE);        

    //cl::CommandQueue queue(context);

    //cl::Program program(context,util::loadProgram("/home/hpc/srihari_hpc/GPUspectralsolver/kernelFunctions.cl"));
    cl::Program program(context,util::loadProgram("/data2/srihari/DDP/GPUspectralsolver/kernelFunctions.cl"));
    
    try{program.build({device});}
    catch(cl::Error error){
        std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)<<"\n";
        exit(1);
    }

    cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,long> findAuxiliaryStress(program,"findAuxiliaryStress");
    cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,long> convolute(program,"convolute");
    cl::make_kernel<cl::Buffer,cl::Buffer,long> getStrainTilde(program,"getStrainTilde");
    cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,long> augmentLagrangianCL(program,"augmentLagrangian");

    clfftPlanHandle planHandleFWD,planHandleBWD;
    clfftDim dim = CLFFT_3D;
    clfftSetupData fftSetup;
    size_t clLengths[3];
    size_t clStridesFWDin[3], clStridesFWDout[3];
    size_t clStridesBWDin[3], clStridesBWDout[3];

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
    prodDimHermitian=(n1/2+1)*n2*n3;
    volumeVoxel=1.0/prodDim;

    stressaux=new double[prodDim*6];
    errors=new double[prodDim*2];

    stressFourier = new vector6_complex[prodDimHermitian];
    ddefgradFourier = new tensor33_complex[prodDimHermitian];

//    delta=new double[n3*n2*n1];
//    out=(fftw_complex *) *fftw_alloc_complex(n3*n2*(n1/2+1));

//    plan_forward=fftw_plan_dft_r2c_3d(n3, n2, n1, delta, out, FFTW_ESTIMATE);
//    plan_backward=fftw_plan_dft_c2r_3d(n3, n2, n1, out, delta, FFTW_ESTIMATE);

    clLengths[0]=n1;
    clLengths[1]=n2;
    clLengths[2]=n3;

    clStridesFWDin[0]=6;
    clStridesFWDin[1]=6*n1;
    clStridesFWDin[2]=6*n1*n2;

    clStridesFWDout[0]=6;
    clStridesFWDout[1]=6*(n1/2+1);
    clStridesFWDout[2]=6*(n1/2+1)*n2;
    
    clStridesBWDin[0]=9;
    clStridesBWDin[1]=9*(n1/2+1);
    clStridesBWDin[2]=9*(n1/2+1)*n2;

    clStridesBWDout[0]=9;
    clStridesBWDout[1]=9*n1;
    clStridesBWDout[2]=9*n1*n2;

    d_stress=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(vector6)*prodDim);
    d_straintilde=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(vector6)*prodDim);
    d_C0_66=cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(tensor66));
    d_strainbar=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(vector6));
    d_fsloc=cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(tensor66)*prodDim);
    d_gammaHat=cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(fourthOrderTensor)*prodDimHermitian);
    d_stressFourier=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(vector6_complex)*prodDimHermitian);
    d_ddefgradFourier=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(tensor33_complex)*prodDimHermitian);
    d_ddefgrad=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(tensor33)*prodDim);
    d_errors=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double)*2*prodDim);

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
    err=clfftSetPlanInStride(planHandleFWD, dim, clStridesFWDin);
    err=clfftSetPlanOutStride(planHandleFWD, dim, clStridesFWDout);

    /* Since we will be using the Hermitian symmetry which is applicable when the input data for the forward transform is all real,
       we choose the plan input and output accordingly. For real transforms clFFT selects the direction of transform based on type
       of the input and output buffers. Thus we create two plan handles one for forward transform and one for reverse transform. */
    err=clfftCopyPlan(&planHandleBWD,context(), planHandleFWD);

    // Specify the type of input and output buffers for forward and reverse transforms
    err=clfftSetLayout(planHandleFWD, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
    err=clfftSetLayout(planHandleBWD, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL);
    
    err=clfftSetPlanBatchSize(planHandleBWD, 9);
    err=clfftSetPlanDistance(planHandleBWD, 1, 1);

    err=clfftSetPlanInStride(planHandleBWD, dim, clStridesBWDin);
    err=clfftSetPlanOutStride(planHandleBWD, dim, clStridesBWDout);

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

        err = queue.enqueueWriteBuffer(d_C0_66, CL_TRUE, 0, sizeof(tensor66), C0_66);

        err = queue.enqueueWriteBuffer(d_gammaHat, CL_TRUE, 0, prodDimHermitian*sizeof(fourthOrderTensor), 
                                    gammaHat);

        err = queue.enqueueWriteBuffer(d_strainbar, CL_TRUE, 0, sizeof(vector6), strainbar);

        err = queue.enqueueWriteBuffer(d_fsloc, CL_TRUE, 0, prodDim*sizeof(tensor66), 
                                       fsloc);

        queue.finish();

        while(iteration<itermax && err2mod > error){
            
            iteration++;

            errstrain=errstress=0.0;

            err = queue.enqueueWriteBuffer(d_stress, CL_TRUE, 0, sizeof(vector6)*prodDim, 
                                        stress);

            err = queue.enqueueWriteBuffer(d_straintilde, CL_TRUE, 0, sizeof(vector6)*prodDim, 
                                        straintilde);

            queue.finish();

            cout<<"--------------------------------------------------------------"<<endl;
            cout<<"ITERATION:"<<iteration<<endl;

            // Arrange data for in
            // Perform forward FFT

            try{findAuxiliaryStress(
                cl::EnqueueArgs(queue,cl::NDRange(prodDim*6),cl::NDRange(6)),
                d_stress, d_straintilde, d_C0_66, prodDim);}
            catch(cl::Error error){
                cout<<error.what()<<"("<<error.err()<<")"<<endl;
            }

            //err = queue.enqueueReadBuffer(d_stress, CL_TRUE, 0, sizeof(vector6)*prodDim, 
                                          //stressaux);

            cout<<"Forward FFT of polarization field"<<endl<<endl;

            err = clfftEnqueueTransform(planHandleFWD, CLFFT_FORWARD, 1, &queue(), 0, NULL, NULL,
                    &d_stress(), &d_stressFourier(), NULL);

            queue.finish();
//            err = queue.enqueueReadBuffer(d_stressFourier, CL_TRUE, 0, sizeof(vector6_complex)*prodDimHermitian, 
//                                          stressFourier);

            cout<<"Gamma convolution"<<endl<<endl;
            try{convolute(
                cl::EnqueueArgs(queue,cl::NDRange(prodDimHermitian)),
                d_stressFourier, d_ddefgradFourier, d_gammaHat, prodDimHermitian);}
            catch(cl::Error error){
                cout<<error.what()<<"("<<error.err()<<")"<<endl;
            }            

            queue.finish();
//            err = queue.enqueueReadBuffer(d_ddefgradFourier, CL_TRUE, 0, sizeof(tensor33_complex)*prodDimHermitian, 
//                                          ddefgradFourier);

            cout<<"Inverse FFT to get deformation gradient"<<endl<<endl;
            err = clfftEnqueueTransform(planHandleBWD, CLFFT_BACKWARD, 1, &queue(), 0, NULL, NULL,
                    &d_ddefgradFourier(), &d_ddefgrad(), NULL);

            try{queue.finish();}
            catch(cl::Error error){
                cout<<error.what()<<"("<<error.err()<<")"<<endl;
            }

//            err = queue.enqueueReadBuffer(d_ddefgrad, CL_TRUE, 0, sizeof(tensor33)*prodDim, 
//                                          ddefgrad);

            try{getStrainTilde(
                cl::EnqueueArgs(queue,cl::NDRange(prodDim)),
                d_ddefgrad, d_straintilde, prodDim);}
            catch(cl::Error error){
                cout<<error.what()<<"("<<error.err()<<")"<<endl;
            }

//            if(iteration==2){
//                for(k=0;k<n3;k++)
//                    for(j=0;j<n2;j++)
//                        for(i=0;i<n1;i++){
//                            cout<<i<<" "<<j<<" "<<k<<endl;
//                            //print1darray(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],9);
//                            print1darray(&straintilde[(k*n2*(n1)+j*(n1)+i)*6],6);
//            }}

            if(iteration==1)
                for(k=0;k<n3;k++)
                    for(j=0;j<n2;j++)
                        for(i=0;i<(n1/2+1);i++){
                            //cout <<i<<" "<<j<<" "<<k<<endl;
                            //debugreversefft for(int count=0;count<6;count++)
                                //debugreversefft cout<<setw(10)<<stressFourier[(k*n2*(n1/2+1)+j*(n1/2+1)+i)].vector[count][0]<<" ";
                            //debugreversefft cout<<endl;
                            //debugreversefft for(int count=0;count<6;count++)
                                //debugreversefft cout<<setw(10)<<stressFourier[(k*n2*(n1/2+1)+j*(n1/2+1)+i)].vector[count][1]<<" ";
                            //debugreversefft cout<<endl<<endl;
                            //for(int count1=0;count1<3;count1++){
                                //for(int count2=0;count2<3;count2++)
                                    //cout<<setw(10)<<ddefgradFourier[(k*n2*(n1/2+1)+j*(n1/2+1)+i)].tensor[count1][count2][0]<<" ";
                                ////cout<<endl;
                            //;}    
                            //cout<<endl<<endl;    
                            //for(int count1=0;count1<3;count1++){
                                //for(int count2=0;count2<3;count2++)
                                    //cout<<setw(10)<<ddefgradFourier[(k*n2*(n1/2+1)+j*(n1/2+1)+i)].tensor[count1][count2][1]<<" ";
                                ////cout<<endl;
                            //;}    
                            //cout<<endl<<endl;    
                        }

            try{augmentLagrangianCL(
                cl::EnqueueArgs(queue,cl::NDRange(6*prodDim),cl::NDRange(6)),
                d_stress, d_straintilde, d_fsloc, d_strainbar, d_C0_66, d_errors, prodDim);}
            catch(cl::Error error){
                cout<<error.what()<<"("<<error.err()<<")"<<endl;
            }

            err = queue.enqueueReadBuffer(d_straintilde, CL_TRUE, 0, sizeof(vector6)*prodDim, 
                                          straintilde);

            err = queue.enqueueReadBuffer(d_stress, CL_TRUE, 0, sizeof(vector6)*prodDim, 
                                          stress);
            
            err = queue.enqueueReadBuffer(d_errors, CL_TRUE, 0, sizeof(double)*2*prodDim, 
                                          errors);

    	    queue.finish();

            cout<<"Augmented Lagrangian method for stress update"<<endl<<endl;

//            augmentLagrangian();

            for(k=0;k<n3;k++)
                for(j=0;j<n2;j++)
                    for(i=0;i<n1;i++){
                        errstrain+=errors[(k*n2*(n1)+j*(n1)+i)*2+0]*volumeVoxel;
                        errstress+=errors[(k*n2*(n1)+j*(n1)+i)*2+1]*volumeVoxel;
            }                        

            for(k=0;k<n3;k++)
                for(j=0;j<n2;j++)
                    for(i=0;i<n1;i++){
                            //print1darray(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],9);
                        //if(iteration==1)
                            //print1darray(&straintilde[(k*n2*(n1)+j*(n1)+i)*6],6);
                        //if(iteration==1)
                            //cout<<i<<" "<<j<<" "<<k<<" actual"<<endl;
                            //print1darray(&stress[(k*n2*(n1)+j*(n1)+i)*6],6);}
                        //else if(iteration==2)
                            //cout<<i<<" "<<j<<" "<<k<<" aux"<<endl;
                            //print1darray(&stressaux[(k*n2*(n1)+j*(n1)+i)*6],6);}
                        //print2darray((double *) C0_66,6);
            }



//            cout<<"Forward FFT of polarization field"<<endl<<endl;
//
//            for(n=0;n<6;n++){
//
//                for(k=0;k<n3;k++)
//                    for(j=0;j<n2;j++)
//                        for(i=0;i<n1;i++){
//                            delta[(k*n2*(n1)+j*(n1)+i)]=stress[(k*n2*(n1)+j*(n1)+i)*6+n];
//                            for(m=0;m<6;m++)
//                                delta[(k*n2*(n1)+j*(n1)+i)]-=C0_66[n][m]*straintilde[(k*n2*(n1)+j*(n1)+i)*6+m];
//                }
//
////                fftw_execute(plan_forward);
//
//                for(k=0;k<n3;k++)
//                    for(j=0;j<n2;j++)
//                        for(i=0;i<(n1/2+1);i++){
//                            work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][0];
//                            workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6+n]=out[k*n2*(n1/2+1)+j*(n1/2+1)+i][1];                            
//                }
//
//            }
//
//            // Convert stress to tensorial form
//            // Multiply with gamma operator
//            cout<<"Gamma convolution"<<endl<<endl;
//            for(k=0;k<n3;k++)
//                for(j=0;j<n2;j++)
//                    for(i=0;i<(n1/2+1);i++){
//                        change_basis(&work[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6], work33, aux66, aux3333, 1);
//                        change_basis(&workim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*6], work33im, aux66, aux3333, 1);
//
//                        multiply3333x33(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9], gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i], work33, 3, 4);
//                        multiply3333x33(&ddefgradim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9], gammaHat[k*n2*(n1/2+1)+j*(n1/2+1)+i], work33im, 3, 4);
//            }
            //
//            //change_basis(strainbar, strainbar33, aux66, aux3333, 1);
//
//            // Arrange data for out
//            cout<<"Inverse FFT to get deformation gradient"<<endl<<endl;
//            for(m=0;m<3;m++)
//                for(n=0;n<3;n++){
//                    for(k=0;k<n3;k++)
//                        for(j=0;j<n2;j++)
//                            for(i=0;i<(n1/2+1);i++){
//                                out[k*n2*(n1/2+1)+j*(n1/2+1)+i][0]=ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9+3*m+n];
//                                out[k*n2*(n1/2+1)+j*(n1/2+1)+i][1]=ddefgradim[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9+3*m+n];
//                    }
//
//                    //out[0][0]=strainbar33[m][n];
//
////                    fftw_execute(plan_backward);
//
//                    for(k=0;k<n3;k++)
//                        for(j=0;j<n2;j++)
//                            for(i=0;i<n1;i++){
//                                ddefgrad[(k*n2*(n1)+j*(n1)+i)*9+3*m+n]=delta[(k*n2*(n1)+j*(n1)+i)]/prodDim;
//                            }
//
//            }
//            // Get symmetric part of defgrad
        //
//            for(k=0;k<n3;k++)
//                for(j=0;j<n2;j++)
//                    for(i=0;i<n1;i++){
//                        symmetric(&ddefgrad[(k*n2*(n1)+j*(n1)+i)*9], (double *) aux33);
//                        change_basis(&straintilde[(k*n2*(n1)+j*(n1)+i)*6], aux33, aux66, aux3333, 2);
//            }
//

            
//            if(iteration==2){
//                for(k=0;k<n3;k++)
//                    for(j=0;j<n2;j++)
//                        for(i=0;i<n1;i++){
//                            cout<<i<<" "<<j<<" "<<k<<endl;
//                            //print1darray(&ddefgrad[(k*n2*(n1/2+1)+j*(n1/2+1)+i)*9],9);
//                            print1darray(&straintilde[(k*n2*(n1)+j*(n1)+i)*6],6);
//                            print1darray(&stress[(k*n2*(n1)+j*(n1)+i)*6],6);
//                }
//            }

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

    return 0;
}
