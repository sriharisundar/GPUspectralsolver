#ifndef matrix_operations
#define matrix_operations

void change_basis(double CE2[6],double C2[3][3],double CE4[6][6],double C4[3][3][3][3],int iopt);

/* 
	iopt=1, Transforms second order matrix c2 into fourth order tensor c4 
	iopt=2, c4 to c2  
	iopt=3, transforms with inv.fact.
	iopt=4, to go from 6x6 to 3x3x3x3 with aijkl antisymmetry */
void voigt(double C2[6][6], double C4[3][3][3][3], int iopt);

/*  a --> transformation matrix
	euler[0]-phi1
	euler[1]-Phi
	euler[2]-phi2
	all angles in degrees
	Calculate the euler angles associated with the transformation
	matrix a(i,j) if iopt=1 and viceversa if iopt=2
	a(i,j) transforms from system sa to system ca.*/
void transformationMatrix(double a[][3],double euler[3],int iopt);

/*
	aIn --> 4th order tensor in reference frame
	aOut --> transformed 4th order tensor
	q --> Transformation matrix for frame a to frame b transform
	Transforms 4th order tensor using transformation matrix q()
	iopt=1 transform from frame a to frame b
	iopt=2 transform from frame b to frame a
*/
void transformFourthOrderTensor(double aIn[3][3][3][3], double aOut[3][3][3][3], double q[3][3], int iopt);

/*
	aIn --> 2nd order tensor in reference frame
	aOut --> transformed 2nd order tensor
	q --> Transformation matrix for frame a to frame b transform
	Transforms 4th order tensor using transformation matrix q()
	iopt=1 transform from frame a to frame b
	iopt=2 transform from frame b to frame a
*/
void transformSecondOrderTensor(double aIn[3][3], double aOut[3][3], double q[3][3], int iopt);

/*
	Find norm of matrix (linearized).
*/
double tnorm(double *a,int nrows, int ncols);

/*
	Use LU decomposition to find the inverse of a matrix.
	Currently implemented for 6x6 matrices only
*/
int findInverse(double *in, int order);


#endif