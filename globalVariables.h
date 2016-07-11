#ifndef global_variables
#define global_variables

#define N1 32
#define N2 32
#define N3 32

extern int n1,n2,n3;

extern double cmat[6][6],cmat33[3][3][3][3];
extern double identityR2[3][3];
extern double identityR4[3][3][3][3];
extern double basis[3][3][6];
extern double euler[N3][N2][N1][3];
extern int grainID[N3][N2][N1],phaseID[N3][N2][N1];

void initglobal(void);

#endif