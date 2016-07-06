#ifndef matrix_operations
#define matrix_operations
void change_basis(double CE2[6][6],double C2[3][3],double CE4[6][6],double C4[3][3][3][3],int IOPT);

/*TRANSFORMS SECOND ORDER MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=1 
IOPT=2, C4 to C2  
IF IOPT=3,TRANSFORMS WITH INV.FACT.
IOPT=4, TO GO FROM 6x6 TO 3x3x3x3 WITH Aijkl ANTISYMMETRY */
void voigt(double C2[6][6], double C4[3][3][3][3], int iopt);

void transformationMatrix(double a[3][3],double euler[3],int iopt);
#endif