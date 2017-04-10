#include "/home/hpc/srihari_hpc/GPUspectralsolver/kernelUtilFunctions.h"
//#include "/data2/srihari/DDP/GPUspectralsolver/kernelUtilFunctions.h"

__kernel void findAuxiliaryStress(
    __global vector6* stress, 
    __global vector6* stressaux, 
    __global vector6* straintilde, 
    __global const tensor66* C0_66, 
    const unsigned long prodDim)
{
    int i,j,m;

    i=get_group_id(0);
    j=get_local_id(0);

    if(i<prodDim){
	stressaux[i].vector[j]=stress[i].vector[j];
        for(m=0;m<6;m++)
            stressaux[i].vector[j]-=C0_66->tensor[j][m]*straintilde[i].vector[m];
    }
}

__kernel void convolute(
    __global vector6_complex* stressFourier, 
    __global tensor33_complex* ddefgradFourier, 
    __global fourthOrderTensor* gammaHat, 
    const unsigned long prodDimHermitian)
{
    int i;

    i=get_global_id(0);
     
    if(i<prodDimHermitian){
        change_basis_fourier(stressFourier,ddefgradFourier,1,i);
        multiply3333x33(ddefgradFourier,gammaHat,3,4,i);
    }

}

__kernel void getStrainTilde(
    __global tensor33* ddefgrad,
    __global vector6* straintilde, 
    const unsigned long prodDim)
{
    int i,k,l;

    tensor33 strain33;

    i=get_global_id(0);

    if(i<prodDim){
        
        for(k=0; k<3;k++)
            for(l=0; l<3;l++)
                strain33.tensor[k][l]=0.5*(ddefgrad[i].tensor[k][l]+ddefgrad[i].tensor[l][k]);

        change_basis(straintilde,&strain33,2,i);
    }
}

__kernel void augmentLagrangian(
    __global vector6* stress, 
    __global vector6* straintilde, 
    __global tensor66* fsloc, 
    __global const vector6* strainbar,
    __global const tensor66* C0_66,
    __global double* errors,
    const unsigned long prodDim)
{
    int i,j,m;

    i=get_group_id(0);
    j=get_local_id(0);

    double volumeVoxel=1.0/prodDim;
    __local double x[6],dg[6],edot[6],dsg[6],ddg[6];

    double errorstress=0,errorstrain=0;    

    if(i<prodDim && j<6){
        dg[j]=strainbar->vector[j]+straintilde[i].vector[j];

        x[j]=stress[i].vector[j];
        for(m=0;m<6;m++)
            x[j]+=C0_66->tensor[j][m]*dg[m];

        edot[j]=0.0;
        for(m=0;m<6;m++)
            edot[j]+=fsloc[i].tensor[j][m]*x[m];


        ddg[j]=dg[j]-edot[j];                
        dsg[j]=0.0;
        for(m=0;m<6;m++)
            dsg[j]+=C0_66->tensor[j][m]*(dg[m]-edot[m]);
        stress[i].vector[j]+=dsg[j];

        if(j==0){
            for(m=0;m<6;m++)
                errorstrain += ddg[m]*ddg[m];
            errorstrain=sqrt(errorstrain);
            errors[i*2+0]=errorstrain;
        }

        if(j==1){
            for(m=0;m<6;m++)
                errorstress += dsg[m]*dsg[m];
            errorstress=sqrt(errorstress);
            errors[i*2+1]=errorstress;
        }

    }
}
