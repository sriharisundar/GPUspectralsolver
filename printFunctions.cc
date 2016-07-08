#include "printFunctions.h"
#include <iostream>
#include <iomanip>

void print2darray(double a[][6]){
	int i,j;

	for(i=0;i<6;i++){
		for(j=0;j<6;j++)
			std::cout<<std::setw(10)<<a[i][j]<<" ";
		std::cout<<std::endl;
	}	
	std::cout<<std::endl;	

}

void print2darray(double a[][3]){
	int i,j;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++)
			std::cout<<std::setw(10)<<a[i][j]<<" ";
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
					std::cout<<std::setw(10)<<a[i][j][k][l]<<" ";
				std::cout<<"\t \t";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;	
	}
	std::cout<<std::endl;	

}
