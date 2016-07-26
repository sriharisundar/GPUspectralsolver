#include <iostream>
#include "matrixOperations.h"
using namespace std;

int main(){

	double **a = new double*[3];
	for(int i=0;i<3;i++)
		a[i]=new double[3];
	
	a[0][0]=5;
	a[1][1]=5;
	a[2][2]=5;

	minv(a,3);

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
			cout<<a[i][j]<<" ";
		cout<<endl;
	}


	return 0;
}