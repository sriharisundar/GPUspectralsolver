#include "matrixOperations.h"
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
					f=1;
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

void euler(double a[3][3],double euler[3],int iopt){

}