#ifndef structures
#define structures
typedef struct {
	double tensor[3][3][3][3];
}fourthOrderTensor;
#endif

#ifndef global_variables
#define global_variables

// Common data for all voxels

extern int n1,n2,n3;
extern double C0_66[6][6];
extern fourthOrderTensor C0_3333,cmat3333;
extern double identityR2[3][3];
extern double identityR4[3][3][3][3];
extern double basis[3][3][6];
extern double strainbar[6],stressbar[6];
extern double stressref,strainref,errstress,errstrain,error;
extern double RVEdim[3];
extern double velgrad33[3][3];
extern double straingradrate33[3][3],straingradrate6[6];
extern double rotationrate33[3][3];
extern double IDstraingradrate[6];
extern int nsteps,itermax;
extern int ictrl,ictrl1,ictrl2;
extern char *outputFile;

//Voxel specific data, needs to be dynamically allocated
extern double *euler;
extern int *grainID,*phaseID;
extern double *ddefgrad,*ddefgradim;
extern double *straintilde,*stress;
extern double *work,*workim;
extern double *cloc,*fsloc;
extern fourthOrderTensor *gammaHat;

void initglobal();

#endif