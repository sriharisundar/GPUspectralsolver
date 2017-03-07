#include "printFunctions.h"
#include <iostream>
#include <iomanip>


void print1darray(double *a,int size){
	int i;
	for(i=0;i<size;i++)
		std::cout<<std::setw(10)<<a[i]<<" ";
	std::cout<<std::endl;	

}

void print2darray(double *a, int size){
	int i,j;

	for(i=0;i<size;i++){
		for(j=0;j<size;j++)
			std::cout<<std::setw(10)<<a[i*size+j]<<" ";
		std::cout<<std::endl;
	}	
	std::cout<<std::endl;	

}

void print3darray(double a[][3][6]){
	int i,j,k;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<6;k++)
				std::cout<<a[i][j][k]<<std::setw(8)<<" ";
			std::cout<<"\t";
		}
		std::cout<<std::endl;
	}	
	std::cout<<std::endl;	

}


void print4darray(double a[3][3][3][3]){
	int i,j,k,l;

	for(i=0;i<3;i++){
		for(k=0;k<3;k++){
			for(j=0;j<3;j++){
				for(l=0;l<3;l++)
					std::cout<<std::setprecision(10)<<std::setw(10)<<a[i][j][k][l]<<" ";
				std::cout<<"\t \t";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;	
	}
	std::cout<<std::endl;	

}
