#include <iostream>
#include "matrixOperations.h"
using namespace std;

int main(){
	
	double det;
	double a[6][6]={0};
	
	a[0][0]=5;
	a[1][1]=5;
	a[2][2]=5;
	a[3][3]=5;
	a[4][4]=5;
	a[5][5]=5;

	findInverse((double *)a,det,6);

	for(int i=0;i<6;i++){
		for(int j=0;j<6;j++)
			cout<<a[i][j]<<" ";
		cout<<endl;
	}

	return 0;
}
