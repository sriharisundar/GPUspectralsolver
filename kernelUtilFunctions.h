typedef struct {
	cl_double tensor[3][3][3][3];
}fourthOrderTensor;

//void change_basis(double CE2[6],double C2[3][3],double CE4[6][6],double C4[3][3][3][3],int iopt){
void change_basis(double* a,double* b,int iopt){

    int i,j,k,l,m,n;

    double rsq2=1.0/sqrt(2);        
    double rsq3=1.0/sqrt(3);        
    double rsq6=1.0/sqrt(6);        

    double basis[3][3][6];

    basis[0][0][1]=-rsq6;
    basis[1][1][1]=-rsq6;
    basis[2][2][1]= 2.0*rsq6;

    basis[0][0][0]=-rsq2;
    basis[1][1][0]= rsq2;

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
                    b[i*3+j]=0;
                    for(k=0;k<6;k++)
                        b[i*3+j]+=a[k]*basis[i][j][k];
                }
            break;

        case 2:
            for(k=0;k<6;k++){
                a[k]=0.0;
                for(i=0;i<3;i++)
                    for(j=0;j<3;j++)
                        a[k]+=b[i*3+j]*basis[i][j][k];
            }
            break;

        case 3:
            for(i=0;i<3;i++)
                for(k=0;k<3;k++)
                    for(j=0;j<3;j++)
                        for(l=0;l<3;l++){
                            b[i*3*3*3+j*3*3+k*3+l]=0.0;
                            for(n=0;n<6;n++)
                                for(m=0;m<6;m++)
                                    b[i*3*3*3+j*3*3+k*3+l]+=a[n*3+m]*basis[i][j][n]*basis[k][l][m];
                        }
            break;

        case 4:
            for(n=0;n<6;n++)
                for(m=0;m<6;m++){
                    a[n*3+m]=0.0;
                    for(i=0;i<3;i++)
                        for(k=0;k<3;k++)
                            for(j=0;j<3;j++)
                                for(l=0;l<3;l++)
                                    a[n*3+m]+=b[i*3*3*3+j*3*3+k*3+l]*basis[i][j][n]*basis[k][l][m];
                }
            break;
    }
}

void multiply3333x33(double *c, fourthOrderTensor A, double b[3][3], int m, int n){
    int i,j,k,l;

    //print2darray(b);
    
    if(m==3 && n==4)
        for(i=0;i<3;i++)
            for(j=0;j<3;j++){
                c[i*3+j]=0;
                for(k=0;k<3;k++)
                    for(l=0;l<3;l++)
                        c[i*3+j]+=A.tensor[i][j][k][l]*b[k][l];
            }

    else if(m==2 && n==4)
        for(i=0;i<3;i++)
            for(k=0;k<3;k++){
                c[i*3+k]=0;
                for(j=0;j<3;j++)
                    for(l=0;l<3;l++)
                        c[i*3+k]+=A.tensor[i][j][k][l]*b[j][l];
                }
}