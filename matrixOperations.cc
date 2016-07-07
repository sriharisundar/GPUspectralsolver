#include "matrixOperations.h"
#define _USE_MATH_DEFINES 
#include <cmath>
void change_basis(double CE2[6][6],double C2[3][3],double CE4[6][6],double C4[3][3][3][3],int iopt){

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

//euler[0]-phi1
//euler[1]-Phi
//euler[2]-phi2
void transformationMatrix(double a[][3],double euler[],int iopt){
	
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