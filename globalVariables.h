#ifndef global_variables
#define global_variables

#define N1 32
#define N2 32
#define N3 32

//extern int n1,n2,n3;

extern double cmat3333[3][3][3][3];
extern double xlsec66[6][6],xlsec3333[3][3][3][3];
extern double identityR2[3][3];
extern double identityR4[3][3][3][3];
extern double basis[3][3][6];
extern double euler[N3][N2][N1][3];
extern int grainID[N3][N2][N1],phaseID[N3][N2][N1];
extern double dbar[6];
extern double velgrad[N3][N2][N1][3][3],velgradim[N3][N2][N1][3][3];
extern double straintilde[N3][N2][N1][6],stress[N3][N2][N1][6];
extern double cloc[N3][N2][N1][6][6],fsloc[N3][N2][N1][6][6];
extern double wgt;
extern double stressref,strainref,errstress,errstrain,error;
extern double RVEdim[3];
extern double velgrad33[3][3];
extern double straingradrate33[3][3],straingradrate6[6];
extern double rotationrate33[3][3];
extern double IDstraingradrate[6];
extern int nsteps,itermax;
extern char *outputFile;

void initglobal(void);

#endif