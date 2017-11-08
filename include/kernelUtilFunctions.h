typedef struct {
    double vector[6];
}vector6;

typedef struct {
    double vector[6][2];
}vector6_complex;

typedef struct {
    double tensor[3][3];
}tensor33;

typedef struct {
    double tensor[3][3][2];
}tensor33_complex;

typedef struct {
    double tensor[6][6];
}tensor66;

typedef struct {
    double tensor[3][3][3][3];
}fourthOrderTensor;

void change_basis_fourier(__global vector6_complex* a,__global tensor33_complex* b, int iopt, int id)
{
    int i,j,k,l,m,n;

    double rsq2=rsqrt(2.0);        
    double rsq3=rsqrt(3.0);        
    double rsq6=rsqrt(6.0);        

    double basis[3][3][6];

    basis[0][0][1]=-rsq6;
    basis[1][1][1]=-rsq6;
    basis[2][2][1]= 2.0*rsq6;

    basis[0][0][0]=-rsq2;
    basis[1][1][0]=rsq2;

    basis[1][2][2]=rsq2;
    basis[2][1][2]=rsq2;

    basis[0][2][3]=rsq2;
    basis[2][0][3]=rsq2;

    basis[0][1][4]=rsq2;
    basis[1][0][4]=rsq2;

    basis[0][0][5]=rsq3;
    basis[1][1][5]=rsq3;
    basis[2][2][5]=rsq3;

    switch(iopt){
        case 1:
            for(i=0;i<3;i++)
                for(j=0;j<3;j++){
                    b[id].tensor[i][j][0]=0;
                    b[id].tensor[i][j][1]=0;
                    for(k=0;k<6;k++){
                        b[id].tensor[i][j][0]+=a[id].vector[k][0]*basis[i][j][k];
                        b[id].tensor[i][j][1]+=a[id].vector[k][1]*basis[i][j][k];
                    }
                }
            break;

        case 2:
            for(k=0;k<6;k++){
                a[id].vector[k][0]=0.0;
                a[id].vector[k][1]=0.0;
                for(i=0;i<3;i++)
                    for(j=0;j<3;j++){
                        a[id].vector[k][0]+=b[id].tensor[i][j][0]*basis[i][j][k];
                        a[id].vector[k][1]+=b[id].tensor[i][j][1]*basis[i][j][k];
                }
            }
            break;
    }
}

void change_basis(__global vector6* a, tensor33* b, int iopt, int id)
{
    int i,j,k,l,m,n;

    double rsq2=rsqrt(2.0);        
    double rsq3=rsqrt(3.0);        
    double rsq6=rsqrt(6.0);        

    double basis[3][3][6];

    basis[0][0][1]=-rsq6;
    basis[1][1][1]=-rsq6;
    basis[2][2][1]= 2.0*rsq6;

    basis[0][0][0]=-rsq2;
    basis[1][1][0]=rsq2;

    basis[1][2][2]=rsq2;
    basis[2][1][2]=rsq2;

    basis[0][2][3]=rsq2;
    basis[2][0][3]=rsq2;

    basis[0][1][4]=rsq2;
    basis[1][0][4]=rsq2;

    basis[0][0][5]=rsq3;
    basis[1][1][5]=rsq3;
    basis[2][2][5]=rsq3;

    switch(iopt){
        case 1:
            for(i=0;i<3;i++)
                for(j=0;j<3;j++){
                    b->tensor[i][j]=0;
                    for(k=0;k<6;k++)
                        b->tensor[i][j]+=a[id].vector[k]*basis[i][j][k];
                }
            break;

        case 2:
            for(k=0;k<6;k++){
                a[id].vector[k]=0.0;
                for(i=0;i<3;i++)
                    for(j=0;j<3;j++)
                        a[id].vector[k]+=b->tensor[i][j]*basis[i][j][k];
            }
            break;
    }
}

void multiply3333x33(__global tensor33_complex* b, __global fourthOrderTensor* A, int m, int n, int id)
{
    int i,j,k,l,realorim;
    
    double c[3][3];

    if(m==3 && n==4)
        for(realorim = 0; realorim < 2; ++realorim){
            for(i=0;i<3;i++)
                for(j=0;j<3;j++){
                    c[i][j]=0;
                    for(k=0;k<3;k++)
                        for(l=0;l<3;l++)
                            c[i][j]+=A[id].tensor[i][j][k][l]*b[id].tensor[k][l][realorim];
                }

            for(i=0;i<3;i++)
                for(j=0;j<3;j++)
                    b[id].tensor[i][j][realorim]=c[i][j];
        }
}