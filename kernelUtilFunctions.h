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
                        b[i*3+j]+=CE2[k]*basis[i][j][k];
                }
            break;

        case 2:
            for(k=0;k<6;k++){
                CE2[k]=0.0;
                for(i=0;i<3;i++)
                    for(j=0;j<3;j++)
                        CE2[k]+=C2[i][j]*basis[i][j][k];
            }
            break;

        case 3:
            for(i=0;i<3;i++)
                for(k=0;k<3;k++)
                    for(j=0;j<3;j++)
                        for(l=0;l<3;l++){
                            C4[i][j][k][l]=0.0;
                            for(n=0;n<6;n++)
                                for(m=0;m<6;m++)
                                    C4[i][j][k][l]+=CE4[n][m]*basis[i][j][n]*basis[k][l][m];
                        }
            break;

        case 4:
            for(n=0;n<6;n++)
                for(m=0;m<6;m++){
                    CE4[n][m]=0.0;
                    for(i=0;i<3;i++)
                        for(k=0;k<3;k++)
                            for(j=0;j<3;j++)
                                for(l=0;l<3;l++)
                                    CE4[n][m]+=C4[i][j][k][l]*basis[i][j][n]*basis[k][l][m];
                }
            break;
    }
}
