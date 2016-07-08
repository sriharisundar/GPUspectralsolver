#ifndef matrix_operations
#define matrix_operations

void change_basis(double CE2[6][6],double C2[3][3],double CE4[6][6],double C4[3][3][3][3],int IOPT);

/* 
TRANSFORMS SECOND ORDER MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=1 
IOPT=2, C4 to C2  
IF IOPT=3,TRANSFORMS WITH INV.FACT.
IOPT=4, TO GO FROM 6x6 TO 3x3x3x3 WITH Aijkl ANTISYMMETRY */
void voigt(double C2[6][6], double C4[3][3][3][3], int iopt);

/*  a --> transformation matrix
	euler[0]-phi1
	euler[1]-Phi
	euler[2]-phi2

CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.*/
void transformationMatrix(double a[][3],double euler[3],int iopt);

void transformFourthOrderTensor(double aIn[3][3][3][3], double aOut[3][3][3][3], double q[3][3]);

void transformSecondOrderTensor(double aIn[3][3], double aOut[3][3], double q[3][3]);

#endif