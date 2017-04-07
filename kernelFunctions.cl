#include "/home/hpc/srihari_hpc/GPUspectralsolver/kernelUtilFunctions.h"

__kernel void findAuxiliaryStress(
    __global vector6* d_stress, 
    __global vector6* d_straintilde, 
    __global const double* d_C0_66, 
    const unsigned int prodDim)
{
    int i,j,m;

    i=get_group_id(0);
    j=get_local_id(0);

    if(i<prodDim){
        for(m=0;m<6;m++)
            d_stress[i].vector[j]-=d_C0_66[j*6+m]*d_straintilde[i].vector[m];
    }
}

__kernel void convolute(
    __global vector6_complex* d_stressFourier, 
    __global tensor33_complex* d_ddefgradFourier, 
    __global fourthOrderTensor* d_gammaHat, 
    const unsigned int prodDimHermitian)
{
    int i;

    i=get_global_id(0);
     
    if(i<prodDimHermitian){
        change_basis_fourier(d_stressFourier,d_ddefgradFourier,1,i);
        multiply3333x33(d_ddefgradFourier,d_gammaHat,3,4,i);
    }

}

__kernel void getStrainTilde(
    __global tensor33* d_ddefgrad,
    __global vector6* d_straintilde, 
    const unsigned int prodDim)
{
    int i,k,l;

    tensor33 strain33;

    i=get_global_id(0);

    if(i<prodDim){
        
        for(k=0; k<3;k++)
            for(l=0; l<3;l++)
                strain33.tensor[k][l]=0.5*(d_ddefgrad[i].tensor[k][l]+d_ddefgrad[i].tensor[l][k]);

        change_basis(d_straintilde,&strain33,2,i);
    }
}

//
//__kernel void augmentLagrangian(){
//
//	cl_double prodDim=n1*n2*n3;
//	cl_double volumeVoxel=1.0/prodDim;
//	cl_double x[6],dg[6],edot[6],dsg[6],ddg[6];
//
//    errstrain=errstress=0.0;
//
//    cl_int i,j,k,n,m;
//
//    for(k=0;k<n3;k++)
//        for(j=0;j<n2;j++)
//            for(i=0;i<n1;i++){
                //
//                for(n=0;n<6;n++)
//                    dg[n]=strainbar[n]+straintilde[(k*n2*(n1)+j*(n1)+i)*6+n];
                //
//                for(n=0;n<6;n++){
//                    x[n]=stress[(k*n2*(n1)+j*(n1)+i)*6+n];
//                    for(m=0;m<6;m++)
//                        x[n]+=C0_66[n][m]*dg[m];
//                }
//
//
//                for(m=0;m<6;m++){
//                    edot[m]=0.0;
//                    for(n=0;n<6;n++)
//                        edot[m]+=fsloc[(k*n2*(n1)+j*(n1)+i)*36+6*m+n]*x[n];
//                }
//
//
//                for(m=0;m<6;m++){
//                    ddg[m]=dg[m]-edot[m];                
//                    dsg[m]=0.0;
//                    for(n=0;n<6;n++)
//                        dsg[m]+=C0_66[m][n]*(dg[m]-edot[m]);
//                    stress[(k*n2*(n1)+j*(n1)+i)*6+m]+=dsg[m];
//                }
//
//                errstrain+=tnorm((double *)ddg,6,1)*volumeVoxel;
//                errstress+=tnorm((double *)dsg,6,1)*volumeVoxel;
//
////                print1darray(stress[k*n2*(n1)+j*(n1)+i],6);
//        }        
//}
/////////
