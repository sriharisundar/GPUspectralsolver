#include "matrixOperations.h"
#include "globalVariables.h"
#define _USE_MATH_DEFINES 
#define TINY 1e-20
#include <cmath>
#include <iostream>
#include <assert.h>

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

double tnorm(double *a,int nrows, int ncols){

	double output=0.0;

	for(int i=0;i<nrows*ncols;i++)
		output += a[i]*a[i];

	output=sqrt(output);

	return output;
}

void forwardsub(double* a, int order, double b[]){
	
	int i,j;
	double sum;

	for(i=0;i<order;i++){
		sum=0;
		for(j=0;j<i;j++)
			sum+=*((a+i*order)+j)*b[j];
		b[i]=(b[i]-sum)*(*((a+i*order)+j));
	}
}

void backwardsub(double* a, int order, double b[]){
	
	int i,j;
	double sum;
	
	for(i=order-1;i>=0;i--){
		sum=0;
		for(j=order-1;j>i;j--)
			sum+=*((a+i*order)+j)*b[j];
		b[i]=(b[i]-sum)/(*((a+i*order)+j));
	}
}


//!  Find pivot element
/*!
*   The function pivot finds the largest element for a pivot in "jcol"
*   of Matrix "a", performs interchanges of the appropriate
*   rows in "a", and also interchanges the corresponding elements in
*   the order vector.
*
*
*  \param     a      -  n by n Matrix of coefficients
*  \param   order  - integer vector to hold row ordering
*  \param    jcol   - column of "a" being searched for pivot element
*
* 	modified from code at
*	http://www.johnloomis.org/ece538/notes/Matrix/ludcmp.html
*/
int pivot(double *a, int order, int jcol)
{
	//int i, ipvt,n;
	//double big, anext;
	//n = a.rows();
//
	///*
	//*  Find biggest element on or below diagonal.
	//*  This will be the pivot row.
	//*/
//
	//ipvt = jcol;
	//big = fabs(*((a+ipvt*order)+ipvt));
	//for (i = ipvt+1; i<n; i++) {
		//anext = fabs(*((a+i*order)+jcol));
		//if (anext>big) {
			//big = anext;
			//ipvt = i;
		//}
	//}
	//assert(fabs(big)>TINY); // otherwise Matrix is singular
//	
	///*
	//*   Interchange pivot row (ipvt) with current row (jcol).
	//*/
//	
	//if (ipvt==jcol) return 0;
	//a.swaprows(jcol,ipvt);
	//i = order[jcol];
	//order[jcol] = order[ipvt];
	//order[ipvt] = i;
	//return 1;
}

//! finds LU decomposition of Matrix
/*!
*   The function ludcmp computes the lower L and upper U triangular
*   matrices equivalent to the A Matrix, such that L U = A.  These
*   matrices are returned in the space of A, in compact form.
*   The U Matrix has ones on its diagonal.  Partial pivoting is used
*   to give maximum valued elements on the diagonal of L.  The order of
*   the rows after pivoting is returned in the integer vector "order".
*   This should be used to reorder the right-hand-side vectors before
*   solving the system A x = b.
*
*
* \param       a     - n by n Matrix of coefficients
* \param      order - integer vector holding row order after pivoting.
*
* 	modified from code at
* 	http://www.johnloomis.org/ece538/notes/Matrix/ludcmp.html
*/

int ludcmp(double* a, int order)
{	return 1;
	//int i, j, k, n, nm1;
	//int flag = 1;    /* changes sign with each row interchange */
	//double sum, diag;
//
	///* establish initial ordering in order vector */
//	
	//for (i=0; i<n; i++) order[i] = i;
//
	///* do pivoting for first column and check for singularity */
//
	//if (pivot(a,order,0)) flag = -flag;
	//diag = 1.0/a[0][0];
	//for (i=1; i<n; i++) a[0][i] *= diag;
//	
	///* 
	//*  Now complete the computing of L and U elements.
	//*  The general plan is to compute a column of L's, then
	//*  call pivot to interchange rows, and then compute
	//*  a row of U's.
	//*/
//	
	//nm1 = n - 1;
	//for (j=1; j<nm1; j++) {
		///* column of L's */
		//for (i=j; i<n; i++) {
			//sum = 0.0;
			//for (k=0; k<j; k++) sum += a[i][k]*a[k][j];
			//a[i][j] -= sum;
		//}
		///* pivot, and check for singularity */
		//if (pivot(a,order,j)) flag = -flag;
		///* row of U's */
		//diag = 1.0/a[j][j];
		//for (k=j+1; k<n; k++) {
			//sum = 0.0;
			//for (i=0; i<j; i++) sum += a[j][i]*a[i][k];
			//a[j][k] = (a[j][k]-sum)*diag;
		//}
	//}
//
	///* still need to get last element in L Matrix */
//
	//sum = 0.0;
	//for (k=0; k<nm1; k++) sum += a[nm1][k]*a[k][nm1];
	//a[nm1][nm1] -= sum;
	//return flag;
}

int lusolve(double* a, int order, double b[]){

	forwardsub(a,order,b);
	backwardsub(a,order,b);
	return 1;
}

int findInverse(double* in, int order){

	int i,j,isingular,bksubstat;
	double A[6][6],B[6];

	for(i=0;i<order;i++){
		for(j=0;j<order;j++)
			A[i][j]=0;
		A[i][i]=1;
	}

	isingular=ludcmp(in,order);

	if(isingular==0)
		return 0;

	for(i=0;i<6;i++)
		bksubstat=lusolve(in,order,A[i]);

	for(i=0;i<order;i++)
		for(j=0;j<order;j++)
			*((in+i*order)+j)=A[i][j];

	return 1;

}