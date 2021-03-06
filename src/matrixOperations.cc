#include "matrixOperations.h"
#include "globalVariables.h"
#include "printFunctions.h"
#define _USE_MATH_DEFINES 
#define TINY 1e-20
#include <cmath>
#include <iostream>
#include <assert.h>
#include <eigen/Eigen/Dense>

using namespace Eigen;

void change_basis(double CE2[6],double C2[3][3],double CE4[6][6],double C4[3][3][3][3],int iopt){

	int i,j,k,l,m,n;

	switch(iopt){
		case 1:
			for(i=0;i<3;i++)
				for(j=0;j<3;j++){
					C2[i][j]=0;
					for(k=0;k<6;k++)
						C2[i][j]+=CE2[k]*basis[i][j][k];
				}
			break;

		case 2:
			for(k=0;k<6;k++){
				CE2[k]=0.0;
				for(i=0;i<3;i++)
					for(j=0;j<3;j++)
						CE2[k]+=C2[i][j]*basis[i][j][k];
			}
			break;

		case 3:
			for(i=0;i<3;i++)
				for(k=0;k<3;k++)
					for(j=0;j<3;j++)
						for(l=0;l<3;l++){
							C4[i][j][k][l]=0.0;
							for(n=0;n<6;n++)
								for(m=0;m<6;m++)
									C4[i][j][k][l]+=CE4[n][m]*basis[i][j][n]*basis[k][l][m];
						}
			break;

		case 4:
			for(n=0;n<6;n++)
				for(m=0;m<6;m++){
					CE4[n][m]=0.0;
					for(i=0;i<3;i++)
						for(k=0;k<3;k++)
							for(j=0;j<3;j++)
								for(l=0;l<3;l++)
									CE4[n][m]+=C4[i][j][k][l]*basis[i][j][n]*basis[k][l][m];
				}
			break;
	}
}

/* 
	iopt=1, Transforms second order matrix c2 into fourth order tensor c4 
	iopt=2, c4 to c2  
	iopt=3, transforms with inv.fact.
	iopt=4, to go from 6x6 to 3x3x3x3 with aijkl antisymmetry */

void voigt(double C2[6][6], double C4[3][3][3][3], int iopt){
	
	int ijv[6][2]={{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};
	int i,j,i1,i2,j1,j2;
	double f;

	switch(iopt){
		case 1 :
			for(i=0;i<6;i++){
				i1=ijv[i][0];
				i2=ijv[i][1];
				for(j=0;j<6;j++){
					j1=ijv[j][0];
					j2=ijv[j][1];
                    C4[i1][i2][j1][j2]=C2[i][j];
			        C4[i2][i1][j1][j2]=C2[i][j];
			        C4[i1][i2][j2][j1]=C2[i][j];
			        C4[i2][i1][j2][j1]=C2[i][j];
			    }
			}
			break;
		
		case 2:
			for(i=0;i<6;i++){
				i1=ijv[i][0];
				i2=ijv[i][1];
				for(j=0;j<6;j++){
					j1=ijv[j][0];
					j2=ijv[j][1];
                    C2[i][j]=C4[i1][i2][j1][j2];
                }
            }
            
            break;

		case 3:           
			for(i=0;i<6;i++){
				i1=ijv[i][0];
				i2=ijv[i][1];
				for(j=0;j<6;j++){
					j1=ijv[j][0];
					j2=ijv[j][1];
					f=1.0;
					if(i>2) f=0.5;
					if(j>2) f=0.5*f;
                    C4[i1][i2][j1][j2]=f*C2[i][j];
			        C4[i2][i1][j1][j2]=f*C2[i][j];
			        C4[i1][i2][j2][j1]=f*C2[i][j];
			        C4[i2][i1][j2][j1]=f*C2[i][j];
				}
			}
			
			break;

		case 4:
			for(i=0;i<6;i++){
				i1=ijv[i][0];
				i2=ijv[i][1];
				for(j=0;j<6;j++){
					j1=ijv[j][0];
					j2=ijv[j][1];
                    
                    if(i<3){
	                    C4[i1][i2][j1][j2]=C2[i][j];
				        C4[i2][i1][j1][j2]=C2[i][j];
				        C4[i1][i2][j2][j1]=C2[i][j];
				        C4[i2][i1][j2][j1]=C2[i][j];
				    }
				    else{
	                    C4[i1][i2][j1][j2]=C2[i][j];
				        C4[i2][i1][j1][j2]=-1*C2[i][j];
				        C4[i1][i2][j2][j1]=C2[i][j];
				        C4[i2][i1][j2][j1]=-1*C2[i][j];				    	
				    }
			    }
			}
			
			break;
	}
}

/*  a --> transformation matrix
	euler[0]-phi1
	euler[1]-Phi
	euler[2]-phi2
	all angles in degrees
	Calculate the euler angles associated with the transformation
	matrix a(i,j) if iopt=1 and viceversa if iopt=2
	a(i,j) transforms from system sa to system ca.*/

void transformationMatrix(double a[][3],double euler[3],int iopt){
	
    double sphi1,cphi1,sPhi,cPhi,sphi2,cphi2;
	double pi=M_PI;

	switch(iopt){

		case 1:
	        euler[1]=acos(a[3][3]);
	        if(std::abs(a[3][3])>=0.9999){	      
	          euler[2]=0.0;
	          euler[0]=atan2(a[1][2],a[1][1]);
	      	}
	        else{
	          sPhi=sin(euler[1]);
	          euler[0]=atan2(a[3][1]/sPhi,-a[3][2]/sPhi);
	          euler[2]=atan2(a[1][3]/sPhi,a[2][3]/sPhi);
	        }
	        euler[0]=euler[0]*180.0/pi;
	        euler[1]=euler[1]*180.0/pi;
	        euler[2]=euler[2]*180.0/pi;
	        
	        break;

	    case 2:
	        sphi1=sin(euler[0]*pi/180.0);
	        cphi1=cos(euler[0]*pi/180.0);
	        sPhi=sin(euler[1]*pi/180.0);
	        cPhi=cos(euler[1]*pi/180.0);
	        sphi2=sin(euler[2]*pi/180.0);
	        cphi2=cos(euler[2]*pi/180.0);
	        
	        a[0][0]=cphi2*cphi1-sphi1*sphi2*cPhi;
	        a[1][0]=-sphi2*cphi1-sphi1*cphi2*cPhi;
	        a[2][0]=sphi1*sPhi;
	        a[0][1]=cphi2*sphi1+cphi1*sphi2*cPhi;
	        a[1][1]=-sphi1*sphi2+cphi1*cphi2*cPhi;
	        a[2][1]=-sPhi*cphi1;
	        a[0][2]=sPhi*sphi2;
	        a[1][2]=cphi2*sPhi;
	        a[2][2]=cPhi;

	        break;
	}
}

/*
	aIn --> 4th order tensor in reference frame
	aOut --> transformed 4th order tensor
	q --> Transformation matrix for frame a to frame b transform
	Transforms 4th order tensor using transformation matrix q()
	iopt=1 transform from frame a to frame b
	iopt=2 transform from frame b to frame a
*/

void transformFourthOrderTensor(double aIn[3][3][3][3], double aOut[3][3][3][3], double q[3][3], int iopt){

	switch(iopt){
		case 1:
			for(int i1=0;i1<3;i1++)
			for(int j1=0;j1<3;j1++)
			for(int k1=0;k1<3;k1++)
			for(int l1=0;l1<3;l1++){
				aOut[i1][j1][k1][l1]=0;		
				for(int i2=0;i2<3;i2++)
				for(int j2=0;j2<3;j2++)
				for(int k2=0;k2<3;k2++)
				for(int l2=0;l2<3;l2++)
					aOut[i1][j1][k1][l1]=aOut[i1][j1][k1][l1]+q[i2][i1]*q[j2][j1]*q[k2][k1]*q[l2][l1]*aIn[i2][j2][k2][l2];
			}
			break;

		case 2:
			for(int i1=0;i1<3;i1++)
			for(int j1=0;j1<3;j1++)
			for(int k1=0;k1<3;k1++)
			for(int l1=0;l1<3;l1++){
				aOut[i1][j1][k1][l1]=0;		
				for(int i2=0;i2<3;i2++)
				for(int j2=0;j2<3;j2++)
				for(int k2=0;k2<3;k2++)
				for(int l2=0;l2<3;l2++)
					aOut[i1][j1][k1][l1]=aOut[i1][j1][k1][l1]+q[i1][i2]*q[j1][j2]*q[k1][k2]*q[l1][l2]*aIn[i2][j2][k2][l2];
			}
	}
}



/*
	aIn --> 2nd order tensor in reference frame
	aOut --> transformed 2nd order tensor
	q --> Transformation matrix for frame a to frame b transform
	Transforms 4th order tensor using transformation matrix q()
	iopt=1 transform from frame a to frame b
	iopt=2 transform from frame b to frame a
*/
void transformSecondOrderTensor(double aIn[3][3], double aOut[3][3], double q[3][3], int iopt){

	switch(iopt){
		case 1:
			for(int i1=0;i1<3;i1++)
			for(int j1=0;j1<3;j1++){
				aOut[i1][j1]=0;
				for(int i2=0;i2<3;i2++)
				for(int j2=0;j2<3;j2++)
					aOut[i1][j1]=aOut[i1][j1]+q[i2][i1]*q[j2][j1]*aIn[i2][j2];
			}
			break;

		case 2:
			for(int i1=0;i1<3;i1++)
			for(int j1=0;j1<3;j1++){
				aOut[i1][j1]=0;
				for(int i2=0;i2<3;i2++)
				for(int j2=0;j2<3;j2++)
					aOut[i1][j1]=aOut[i1][j1]+q[i1][i2]*q[j1][j2]*aIn[i2][j2];
			}
	}
}

/*
	Find norm of matrix (linearized)
*/
double tnorm(double *a,int nrows, int ncols){

	double output=0.0;

	for(int i=0;i<nrows*ncols;i++)
		output += a[i]*a[i];

	output=sqrt(output);

	return output;
}

/*
	c --> Output tensor
	A --> input tensor 4th order
	b --> input tensor 2nd order
	m,n --> non repeated index numbers in right hand side
	ex: If you want to do c_ij=A_ijkl:b_jl, then m=2, n=4
*/
void multiply3333x33(double *c, fourthOrderTensor A, double b[3][3], int m, int n){
	int i,j,k,l;

	//print2darray(b);
	
	if(m==3 && n==4)
		for(i=0;i<3;i++)
		    for(j=0;j<3;j++){
		    	c[i*3+j]=0;
	    		for(k=0;k<3;k++)
				    for(l=0;l<3;l++)
				    	c[i*3+j]+=A.tensor[i][j][k][l]*b[k][l];
			}

	else if(m==2 && n==4)
		for(i=0;i<3;i++)
		    for(k=0;k<3;k++){
		    	c[i*3+k]=0;
	    		for(j=0;j<3;j++)
				    for(l=0;l<3;l++)
				    	c[i*3+k]+=A.tensor[i][j][k][l]*b[j][l];
				}
}

int findInverse(double* in, double det, int order){

	int i,j;

	typedef Matrix <double, Dynamic, Dynamic> MatrixXd;
	MatrixXd matin(order,order);
	MatrixXd matout(order,order);

	for(i=0;i<order;i++)
		for(j=0;j<order;j++)
			matin(i,j)=in[j*order+i];

	matout=matin.inverse();
	//det=matin.determinant();

	for(i=0;i<order;i++)
		for(j=0;j<order;j++)
			in[j*order+i]=matout(i,j);


	return 1;
}

/*
	Find symmetric part of a matrix
*/
void symmetric(double *aIn, double *aOut){

	for(int i=0; i<3;i++)
		for(int j=0; j<3;j++)
			aOut[i*3+j]=0.5*(aIn[i*3+j]+aIn[j*3+i]);

}

/*
	Find antisymmetric part of a matrix
*/
void antisymmetric(double aIn[3][3], double aOut[3][3]){

	for(int i=0; i<3;i++)
		for(int j=0; j<3;j++)
			aOut[i][j]=0.5*(aIn[i][j]-aIn[j][i]);

}